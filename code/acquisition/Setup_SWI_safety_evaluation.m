%% Description
%Copy of SetUp_SWI_Widebeam where only non-steered beams are included for
%safety measurements. This script repeats a single transmit event,
%determined by TestedMode and TestedFreq

clear all;close all force;clc;

%% System parameters
filename = ('S5_1_SWI_Luuk');
ReprocessRcvData = false;

%Everything is derived from the following two fields. 
TestedMode = "WB"; %WB, FC, DW, SW_push, SW_tracking
TestedFreq = "Positive"; %Fundamental, Positive (harmonic), Negative (harmonic), Coded
OrientationMode = "Focused"; %Focused or Widebeam. MUST be set to match the clinical SetUp_SWI_Widebeam.m run (controls the widebeam na/txFocus/FPS).

% Specify Trans structure array.
Trans = iusInitTrans_Max('S5-1', [], []);
Trans.frequency = 3.125;
%Trans.Bandwidth = [1,5];
Trans.Bandwidth = [1,6];
lambda_mm = 1540/(Trans.frequency*1e6)*1e3;

%General parameters
PulsesPerBarker = 1; %Each element in the Barker code corresponds to a number of cycles (default: 1)
UpsampleFactor = 1; %Include zero cycles in excitation pulse to clearly separate each Barker element
CodedExcitationType = "Chirp"; %"Barker" or "Chirp"
NumCyclesChirp = 10;
FundamentalTransmitFrequency = 3.25;
HarmonicTransmitFrequency = 1.9531; %Match clinical SetUp_SWI_Widebeam.m (1.9531)

%All B-mode imaging modes share the same start and end depth
Bmode.startDepth = 0;
Bmode.endDepth = 250;
Bmode.PRI = 270; %Match clinical SetUp_SWI_Widebeam.m base PRI

GeneralWindow = "rectangular"; %rectangular or hanning or kaiser

%Tune-able parameters widebeam B-mode imaging
%Note that the amount of focusing and angle of the PData regions should
%depend on the number of widebeams that are used.

switch(TestedFreq)
    case "Fundamental"
        Bmode_WB.HarmonicImaging = false;
        Bmode_WB.CodedExcitation = false;
    case "Positive"
        Bmode_WB.HarmonicImaging = true;
        Bmode_WB.CodedExcitation = false;
    case "Negative"
        Bmode_WB.HarmonicImaging = true;
        Bmode_WB.CodedExcitation = false;
    case "Coded"
        Bmode_WB.HarmonicImaging = false;
        Bmode_WB.CodedExcitation = true;
end

Bmode_WB.theta = (pi/180)*80;  % 80 degree sector
%Widebeam parameters matched to clinical SetUp_SWI_Widebeam.m (governed by OrientationMode)
if(OrientationMode=="Focused")
    Bmode_WB.DesiredFPS = 90;
    Bmode_WB.na = 21; % 21
    Bmode_WB.rayDelta = Bmode_WB.theta/(Bmode_WB.na-1); % angle increment between rays
    Bmode_WB.RegionAngle = Bmode_WB.rayDelta*3; %Angle that each Pdata region spans.
    Bmode_WB.txFocus = -250; %-700 or 600 (closer to 0 is weaker focus and wider beam)
else
    Bmode_WB.DesiredFPS = 30;
    Bmode_WB.na = 73; % 73
    Bmode_WB.rayDelta = Bmode_WB.theta/(Bmode_WB.na-1); % angle increment between rays
    Bmode_WB.RegionAngle = Bmode_WB.rayDelta*6; %Angle that each Pdata region spans.
    Bmode_WB.txFocus = 600; %-700 or 600 (closer to 0 is weaker focus and wider beam)
end
Bmode_WB.PRI_us = Bmode.PRI; %Match clinical (270)
Bmode_WB.S_curve = 0; %0 for no S-curve, 1 almost linear, 5 most S-curve, "Custom" for manually saved S-curve
Bmode_WB.TxWindowOption = "hanning"; %Match clinical (hanning). Also used for shear wave tracking.
%Bmode_WB.Nframes = 2;

switch(Bmode_WB.TxWindowOption)
    case "rectangular"
        Bmode_WB.TxWindow = ones(1,Trans.numelements);
    case "hanning"
        Bmode_WB.TxWindow = hanning(Trans.numelements)';
    case "kaiser"
        Bmode_WB.TxWindow = kaiser(Trans.numelements,1)';
end

%Fixed/derived
Bmode_WB.startDepth = Bmode.startDepth;
Bmode_WB.endDepth = Bmode.endDepth;
Bmode_WB.aperture = 79*Trans.spacing; % aperture based on 80 elements
Bmode_WB.dApex = (Bmode_WB.aperture/4)/tan(Bmode_WB.theta/2); % dist. to virt. apex
Bmode_WB.maxAcqLength = ceil(Bmode_WB.endDepth/cos(Bmode_WB.theta/2));
if(Bmode_WB.HarmonicImaging)
    Bmode_WB.ttnFrame = max([Bmode_WB.PRI_us, 1e6/Bmode_WB.DesiredFPS-2*Bmode_WB.na*Bmode_WB.PRI_us]);
    Bmode_WB.ActualFPS = 1e6/(2*Bmode_WB.na*Bmode_WB.PRI_us+Bmode_WB.ttnFrame);
    Bmode_WB.MaxNa = floor(1e6/(2*Bmode_WB.PRI_us*Bmode_WB.DesiredFPS));
else
    Bmode_WB.ttnFrame = max([Bmode_WB.PRI_us, 1e6/Bmode_WB.DesiredFPS-Bmode_WB.na*Bmode_WB.PRI_us]);
    Bmode_WB.ActualFPS = 1e6/(Bmode_WB.na*Bmode_WB.PRI_us+Bmode_WB.ttnFrame);
    Bmode_WB.MaxNa = floor(1e6/(Bmode_WB.PRI_us*Bmode_WB.DesiredFPS));
end
if(~isfield(Bmode_WB,'Nframes')); Bmode_WB.Nframes = 2*ceil(Bmode_WB.ActualFPS/2); end %Set number of frames to frames in 1 second

fprintf("Widebeams currently use %i beams with a PRI of %i to realize a frame rate of %.1f FPS. The number of angles can be at most %i to reach %i FPS.\n",Bmode_WB.na,Bmode_WB.PRI_us,Bmode_WB.ActualFPS,Bmode_WB.MaxNa,Bmode_WB.DesiredFPS)

%Tune-able parameters focused B-mode
switch(TestedFreq)
    case "Fundamental"
        Bmode_FC.HarmonicImaging = false;
        Bmode_FC.CodedExcitation = false;
    case "Positive"
        Bmode_FC.HarmonicImaging = true;
        Bmode_FC.CodedExcitation = false;
    case "Negative"
        Bmode_FC.HarmonicImaging = true;
        Bmode_FC.CodedExcitation = false;
    case "Coded"
        Bmode_FC.HarmonicImaging = false;
        Bmode_FC.CodedExcitation = true;
end

Bmode_FC.na = 73; %73
Bmode_FC.PRI_us = Bmode.PRI; %Match clinical (270)
Bmode_FC.txFocus = 160; % initial focal depth
Bmode_FC.DesiredFPS = 30; % Preferred frame rate
Bmode_FC.txFNum = 3.5;  % set to desired f-number value for transmit (range: 1.0 - 20)
Bmode_FC.TxWindowOption = "kaiser"; %Match clinical (kaiser)
Bmode_FC.Nframes = 2;

%Fixed/derived
Bmode_FC.startDepth = Bmode.startDepth;
Bmode_FC.endDepth = Bmode.endDepth;
Bmode_FC.theta = -Bmode_WB.theta/2; %Angle to span with focused transmits
Bmode_FC.rayDelta = 2*(-Bmode_FC.theta)/(Bmode_FC.na-1); %Angle step between rays
Bmode_FC.aperture = Trans.numelements*Trans.spacing;
Bmode_FC.radius = (Bmode_FC.aperture/2)/tan(-Bmode_FC.theta); % dist. to virt. apex
Bmode_FC.Ntx_half=round((Bmode_FC.txFocus/Bmode_FC.txFNum)/Trans.spacing/2); % no. of elements in 1/2 aperture.
if Bmode_FC.Ntx_half > (Trans.numelements/2 - 1), Bmode_FC.Ntx_half = floor(Trans.numelements/2 - 1); end
Bmode_FC.Ntx= 2*Bmode_FC.Ntx_half+1;
switch(Bmode_FC.TxWindowOption)
    case "rectangular"
        Bmode_FC.TxWindow = ones(1,Bmode_FC.Ntx);
    case "hanning"
        Bmode_FC.TxWindow = hanning(Bmode_FC.Ntx);
    case "kaiser"
        Bmode_FC.TxWindow = kaiser(Bmode_FC.Ntx,1);
end
Bmode_FC.maxAcqLength = ceil(Bmode_FC.endDepth/cos(Bmode_FC.theta/2));
if(Bmode_FC.HarmonicImaging)
    Bmode_FC.ttnFrame = max([Bmode_FC.PRI_us, 1e6/Bmode_FC.DesiredFPS-2*Bmode_FC.na*Bmode_FC.PRI_us]);
    Bmode_FC.ActualFPS = 1e6/(2*Bmode_FC.na*Bmode_FC.PRI_us+Bmode_FC.ttnFrame);
    Bmode_FC.MaxNa = floor(1e6/(2*Bmode_FC.PRI_us*Bmode_FC.DesiredFPS));
else
    Bmode_FC.ttnFrame = max([Bmode_FC.PRI_us, 1e6/Bmode_FC.DesiredFPS-Bmode_FC.na*Bmode_FC.PRI_us]);
    Bmode_FC.ActualFPS = 1e6/(Bmode_FC.na*Bmode_FC.PRI_us+Bmode_FC.ttnFrame);
    Bmode_FC.MaxNa = floor(1e6/(Bmode_FC.PRI_us*Bmode_FC.DesiredFPS));
end
if(~isfield(Bmode_FC,'Nframes')); Bmode_FC.Nframes = 2*ceil(Bmode_FC.ActualFPS/2); end %Set number of frames to frames in 1 second

fprintf("Focused beams currently use %i beams with a PRI of %i to realize a frame rate of %.1f FPS. The number of angles can be at most %i to reach %i FPS.\n",Bmode_FC.na,Bmode_FC.PRI_us,Bmode_FC.ActualFPS,Bmode_FC.MaxNa,Bmode_FC.DesiredFPS)

%Tune-able parameters ultrafast B-mode
switch(TestedFreq)
    case "Fundamental"
        Bmode_DW.HarmonicImaging = false;
        Bmode_DW.CodedExcitation = false;
    case "Positive"
        Bmode_DW.HarmonicImaging = true;
        Bmode_DW.CodedExcitation = false;
    case "Negative"
        Bmode_DW.HarmonicImaging = true;
        Bmode_DW.CodedExcitation = false;
    case "Coded"
        Bmode_DW.HarmonicImaging = false;
        Bmode_DW.CodedExcitation = true;
end
Bmode_DW.HarmonicOverlap = false; %If true, we use a 3 transmit angles where the center one has opposite polarity and the other two have 0.5 weight.
Bmode_DW.na = 1;

if(Bmode_DW.HarmonicOverlap)
    Bmode_DW.CodedExcitation = false;
    Bmode_DW.HarmonicImaging = false;
    assert(Bmode_DW.na==3,"When using harmonic overlap in diverging waves, the number of angles should be 3.")
end

if(Bmode_DW.na==1)
    Bmode_DW.StartAngle = 0; %Degree
else
    Bmode_DW.StartAngle = 6; %Degree
end
Bmode_DW.DesiredFPS = 5000;
Bmode_DW.PRI_us = Bmode.PRI; %Match clinical (270)
Bmode_DW.TxWindowOption = "kaiser"; %Match clinical (kaiser)
%Bmode_DW.Nframes = 2;

switch(Bmode_DW.TxWindowOption)
    case "rectangular"
        Bmode_DW.TxWindow = ones(1,Trans.numelements);
    case "hanning"
        Bmode_DW.TxWindow = hanning(Trans.numelements)';
    case "kaiser"
        Bmode_DW.TxWindow = kaiser(Trans.numelements,1)';
end

%Fixed/derived
Bmode_DW.Angles = deg2rad(linspace(Bmode_DW.StartAngle,-Bmode_DW.StartAngle,Bmode_DW.na));
Bmode_DW.theta = -Bmode_WB.theta;
Bmode_DW.fullAngle = -Bmode_WB.theta;
Bmode_DW.aperture = Trans.numelements*Trans.spacing;
Bmode_DW.radius = (Bmode_DW.aperture/2)/tan(-Bmode_DW.theta);
Bmode_DW.txFocus = -25; 
%Bmode_DW.txFocus = -4*Bmode_DW.radius; 
Bmode_DW.startDepth = 0;
Bmode_DW.endDepth = Bmode.endDepth;
%Bmode_DW.endDepth = Bmode.endDepth;
Bmode_DW.maxAcqLength = ceil(Bmode_DW.endDepth/cos(Bmode_DW.theta/2));
if(Bmode_DW.HarmonicImaging)
    Bmode_DW.ttnFrame = max([Bmode_DW.PRI_us, 1e6/Bmode_DW.DesiredFPS-2*Bmode_DW.na*Bmode_DW.PRI_us]);
    Bmode_DW.ActualFPS = 1e6/(2*Bmode_DW.na*Bmode_DW.PRI_us+Bmode_DW.ttnFrame);
    Bmode_DW.MaxNa = floor(1e6/(2*Bmode_DW.PRI_us*Bmode_DW.DesiredFPS));
else
    Bmode_DW.ttnFrame = max([Bmode_DW.PRI_us, 1e6/Bmode_DW.DesiredFPS-Bmode_DW.na*Bmode_DW.PRI_us]);
    Bmode_DW.ActualFPS = 1e6/(Bmode_DW.na*Bmode_DW.PRI_us+Bmode_DW.ttnFrame);
    Bmode_DW.MaxNa = floor(1e6/(Bmode_DW.PRI_us*Bmode_DW.DesiredFPS));
end
if(~isfield(Bmode_DW,'Nframes')); Bmode_DW.Nframes = 2*ceil(Bmode_DW.ActualFPS/2); end %Set number of frames to frames in 1 second

fprintf("Diverging waves currently use %i angles with a PRI of %i to realize a frame rate of %.1f FPS. The number of angles can be at most %i to reach %i FPS.\n",Bmode_DW.na,Bmode_DW.PRI_us,Bmode_DW.ActualFPS,Bmode_DW.MaxNa,Bmode_DW.DesiredFPS)

%Tune-able parameters shear wave imaging
switch(TestedFreq)
    case "Fundamental"
        SW.HarmonicImaging = false;
        SW.CodedExcitation = false;
    case "Positive"
        SW.HarmonicImaging = true;
        SW.CodedExcitation = false;
    case "Negative"
        SW.HarmonicImaging = true;
        SW.CodedExcitation = false;
    case "Coded"
        SW.HarmonicImaging = false;
        SW.CodedExcitation = true;
end

SW.startDepth = 0;
SW.endDepth = Bmode.endDepth;
SW.DesiredFPS = 20; %Match clinical
SW.PushApodization = "none"; %kaiser none
SW.PushApodBeta = 2;
SW.WaitForRpeak = false;
SW.DoAutomaticSWS = false;
SW.ObserveLeftOfPush = true;
SW.Ndetect          = 60;      % Set SW.Ndetect = number of detect acquisitions.
SW.PushFrequency = 2.25;
SW.pushCycle  = 1500; %Default 1500
SW.nb_push_elmts = 41; %All? 41 or 79
SW.PRI_us = Bmode.PRI; %Match clinical (270)
SW.ClutterFilter = "None";
SW.ClutterFilterN = 3;
SW.Nframes   = 6;

SW.txFocus = 600; %Tracking only. Closer to 0 = more diverging beam
SW.aperture = Bmode_WB.aperture;
SW.dApex = Bmode_WB.dApex;

%Fixed/derived
SW.powermax   = 250;      % scaling of the display function
assert(mod(SW.nb_push_elmts,2)==1,"Number of push elements must be odd")
SW.maxAcqLength = ceil(SW.endDepth/cos(Bmode_WB.theta/2));
if(SW.HarmonicImaging)
    SW.ttnFrame = max([SW.PRI_us, 1e6/SW.DesiredFPS-2*SW.Ndetect*SW.PRI_us]);
    SW.ActualFPS = 1e6/(2*SW.Ndetect*SW.PRI_us+SW.ttnFrame);
else
    SW.ttnFrame = max([SW.PRI_us, 1e6/SW.DesiredFPS-SW.Ndetect*SW.PRI_us]);
    SW.ActualFPS = 1e6/(SW.Ndetect*SW.PRI_us+SW.ttnFrame);
end
if(~isfield(SW,'Nframes')); SW.Nframes = 2*ceil(SW.ActualFPS/2); end %Set number of frames to frames in 1 second

% Define ROI for IQ processing, can be changed in GUI
SW.FocusX_wavelength = 0; %10
SW.FocusX_cm = SW.FocusX_wavelength*1540/3.125e6*100;
SW.FocusZ_wavelength = round(40/lambda_mm);
SW.FocusZ_cm = SW.FocusZ_wavelength*1540/3.125e6*100;
SW.ROI_width = 80;
SW.ROI_height = 60;

if(Bmode_WB.HarmonicImaging)
    pgn = 2; %increased processing gain for harmonic imaging
else
    pgn = 1;
end

if(SW.HarmonicImaging)
    SW.Ndetect = SW.Ndetect/2; %decreased number of acquisitions for SW tracking (since each acquisition now requires two transmits)
end

%Tune-able parameters for strain estimation
%Bmode_strain.Nframes = 3*Bmode_WB.ActualFPS;
Bmode_strain.Nframes = 2;

%% Define system parameters.
Resource.Parameters.numTransmit = 80;      % number of transmit channels.
Resource.Parameters.numRcvChannels = 80;    % number of receive channels.
Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.VDAS.dmaTimeout = 15000; %Increase DMA timeout
Resource.Parameters.simulateMode = 0;
%Resource.Parameters.simulateMode = 1; %forces simulate mode, even if hardware is present.
%Resource.Parameters.simulateMode = 2 %stops sequence and processes RcvData continuously.

%% Imaging parameters
[PData,DisplayWindow] = CalculatePData(Bmode_WB,Bmode_FC,Bmode_DW,SW,"init");

% Specify Media object. 'pt1.m' script defines array of point targets.
Media.MP(1,:) = [-45,0,30,1.0];
Media.MP(2,:) = [-15,0,30,1.0];
Media.MP(3,:) = [15,0,30,1.0];
Media.MP(4,:) = [45,0,30,1.0];
Media.MP(5,:) = [-15,0,60,1.0];
Media.MP(6,:) = [-15,0,90,1.0];
Media.MP(7,:) = [-15,0,120,1.0];
Media.MP(8,:) = [-15,0,150,1.0];
Media.MP(9,:) = [-45,0,120,1.0];
Media.MP(10,:) = [15,0,120,1.0];
Media.MP(11,:) = [45,0,120,1.0];
Media.MP(12,:) = [-10,0,69,1.0];
Media.MP(13,:) = [-5,0,75,1.0];
Media.MP(14,:) = [0,0,78,1.0];
Media.MP(15,:) = [5,0,80,1.0];
Media.MP(16,:) = [10,0,81,1.0];
Media.MP(17,:) = [-75,0,120,1.0];
Media.MP(18,:) = [75,0,120,1.0];
Media.MP(19,:) = [-15,0,180,1.0];
Media.numPoints = 19;
Media.attenuation = -0.5;
Media.function = 'movePoints';

%Single, stationary point
% Media.MP(1,:) = [0,0,100,1.0];
% Media.numPoints = 1;
% Media.attenuation = -0.5;

%% Specify Resources.
% RcvBuffer for widebeam B-mode data
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 4096; %4096
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 10;

% InterBuffer for widebeam Bmode
Resource.InterBuffer(1).numFrames = 1;
Resource.InterBuffer(1).pagesPerFrame = 1;

% ImageBuffer for reference Bmode image
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 10;

%Displaywindow
Resource.DisplayWindow = DisplayWindow;


%% Transmit parameters
% Specify Transmit waveform structure.
% - detect waveform positive

%Define multiple TW structures
%1: Normal B-mode pulse
TW(1).type = 'parametric';
TW(1).Parameters = [FundamentalTransmitFrequency,0.67,2,1];

%2: Coded excitation
switch(CodedExcitationType)
    case "Barker"
        %Length 13 Barker
        TW(2).type = 'envelope';
        TW(2).envNumCycles = 13*PulsesPerBarker*UpsampleFactor;
        TW(2).envFrequency = ones(1,TW(2).envNumCycles,1)*FundamentalTransmitFrequency;
        TW(2).envPulseWidth = 0.67.*[+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1];
        TW(2).envPulseWidth = upsample(TW(2).envPulseWidth, UpsampleFactor);
        TW(2).envPulseWidth = repelem(TW(2).envPulseWidth,PulsesPerBarker);
        %TW(2).envPulseWidth = hanning(TW(2).envNumCycles)'.*TW(2).envPulseWidth;
    case "Chirp"
        TW(2).type = 'envelope';
        TW(2).envNumCycles = NumCyclesChirp;
        TW(2).envFrequency = linspace(FundamentalTransmitFrequency-0.75,FundamentalTransmitFrequency+0.75,TW(2).envNumCycles);
        %TW(2).envFrequency = flip(TW(2).envFrequency);
        TW(2).envPulseWidth = 0.67.*ones(1,TW(2).envNumCycles);
        %TW(2).envPulseWidth = hanning(TW(2).envNumCycles)'.*TW(2).envPulseWidth;
end

%3: Push pulse
TW(3).type = 'parametric';
TW(3).Parameters = [SW.PushFrequency,0.67,SW.pushCycle*2,1]; %Match clinical (duty 0.67)

%4: Harmonic B-mode pulse positive polarity
TW(4).type = 'parametric';
TW(4).Parameters = [HarmonicTransmitFrequency,0.67,2,1];

%5: Harmonic B-mode pulse negative polarity
TW(5).type = 'parametric';
TW(5).Parameters = [HarmonicTransmitFrequency,0.67,2,-1];

%6: Harmonic coded excitation positive polarity
switch(CodedExcitationType)
    case "Barker"
        %Length 13 Barker
        TW(6).type = 'envelope';
        TW(6).envNumCycles = 13*PulsesPerBarker*UpsampleFactor;
        TW(6).envFrequency = ones(1,TW(6).envNumCycles,1)*HarmonicTransmitFrequency;
        TW(6).envPulseWidth = 0.67.*[+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1];
        TW(6).envPulseWidth = upsample(TW(6).envPulseWidth, UpsampleFactor);
        TW(6).envPulseWidth = repelem(TW(6).envPulseWidth,PulsesPerBarker);
        %TW(6).envPulseWidth = hanning(TW(6).envNumCycles)'.*TW(6).envPulseWidth;
    case "Chirp"
        TW(6).type = 'envelope';
        TW(6).envNumCycles = NumCyclesChirp;
        TW(6).envFrequency = linspace(HarmonicTransmitFrequency-0.5,HarmonicTransmitFrequency+0.5,TW(6).envNumCycles);
        TW(6).envPulseWidth = 0.67.*ones(1,TW(6).envNumCycles);
        %TW(6).envPulseWidth = hanning(TW(6).envNumCycles)'.*TW(6).envPulseWidth;
end

%7: Harmonic coded excitation negative polarity
TW(7) = TW(6);
TW(7).envPulseWidth = -TW(7).envPulseWidth;

%Compute all TW waveforms
[~, ~, ~, ~, Computed_TW(2)] = computeTWWaveform(TW(2));
[~, ~, ~, ~, Computed_TW(6)] = computeTWWaveform(TW(6));
[~, ~, ~, ~, Computed_TW(7)] = computeTWWaveform(TW(7));
Computed_TW(2).equalize = []; %Otherwise an error occurs that this field missing
[~, ~, ~, ~, Computed_TW(1)] = computeTWWaveform(TW(1));
[~, ~, ~, ~, Computed_TW(3)] = computeTWWaveform(TW(3));
[~, ~, ~, ~, Computed_TW(4)] = computeTWWaveform(TW(4));
[~, ~, ~, ~, Computed_TW(5)] = computeTWWaveform(TW(5));
TW = Computed_TW;

% Set TPC profile 5 high voltage limit.
TPC(1).hv = 25;

TPC(5).maxHighVoltage = Trans.maxHighVoltage;
TPC(5).hv = 25; %35

% Specify TX structure array.
TX = repmat(struct('waveform', 1, ...
    'Origin', [0.0,0.0,0.0], ...
    'focus', Bmode_WB.txFocus, ...
    'Steer', [0.0,0.0], ...
    'Apod', ones(1,80), ...
    'Delay', zeros(1,80), ...
    'focusX',[],...
    'pushElements',[],...
    'FocalPt',[],...
    'TXPD', [], ...
    'peakCutOff', 1, ...
    'peakBLMax', 500), 1, 2*Bmode_WB.na+2*Bmode_FC.na+2*Bmode_DW.na+2+1);


% - Set event specific TX attributes for B-mode imaging
TXorgs = Bmode_WB.dApex*tan(Bmode_WB.Angles);
h = waitbar(0,'Program TX parameters, please wait!');
[Bmode_WB,TXwaveform1,TXwaveform2] = SelectTXwaveforms(Bmode_WB);
for n = 1:Bmode_WB.na
    TX(n).waveform = TXwaveform1;
    TX(n).Origin = [TXorgs(n), 0.0, 0.0];
    TX(n).Apod = Bmode_WB.TxWindow;
    TX(n).Steer = [Bmode_WB.Angles(n),0.0];
    TX(n).Delay = computeTXDelays(TX(n));
    TX(n).TXPD = computeTXPD(TX(n),PData);

    TX(n+Bmode_WB.na).waveform = TXwaveform2;
    TX(n+Bmode_WB.na).Origin = TX(n).Origin;
    TX(n+Bmode_WB.na).Apod = TX(n).Apod;
    TX(n+Bmode_WB.na).Steer = TX(n).Steer;
    TX(n+Bmode_WB.na).Delay = TX(n).Delay;
    TX(n+Bmode_WB.na).TXPD = TX(n).TXPD;
    waitbar(n/Bmode_WB.na)
end
close(h)

%Focused B-mode
Bmode_FC.Angles = Bmode_FC.theta:Bmode_FC.rayDelta:(Bmode_FC.theta + (Bmode_FC.na-1)*Bmode_FC.rayDelta);
Bmode_FC.TXorgs = Bmode_FC.radius*tan(Bmode_FC.Angles);

k = 2*Bmode_WB.na;
[Bmode_FC,TXwaveform1,TXwaveform2] = SelectTXwaveforms(Bmode_FC);
for n = k+1:k+Bmode_FC.na   % P.numRays transmit events
    TX(n).waveform = TXwaveform1;
    TX(n).Origin = [Bmode_FC.TXorgs(n-k),0.0,0.0];
    TX(n).focus = Bmode_FC.txFocus;
    TX(n).Steer = [Bmode_FC.Angles(n-k),0.0];
    TX(n).ApertureCenter = max([min([find(Trans.ElementPos/lambda_mm>=Bmode_FC.TXorgs(n-k),1),Trans.numelements-Bmode_FC.Ntx_half]),Bmode_FC.Ntx_half+1]);
    TX(n).Apod = zeros(1,Trans.numelements);
    TX(n).Apod(TX(n).ApertureCenter-Bmode_FC.Ntx_half:TX(n).ApertureCenter+Bmode_FC.Ntx_half) = Bmode_FC.TxWindow;
    TX(n).Delay = computeTXDelays(TX(n));

    TX(n+Bmode_FC.na) = TX(n);
    TX(n+Bmode_FC.na).waveform = TXwaveform2;
end

%Diverging wave
k = 2*Bmode_WB.na+2*Bmode_FC.na;
[Bmode_DW,TXwaveform1,TXwaveform2] = SelectTXwaveforms(Bmode_DW);
for n = k+1: k+Bmode_DW.na
    TX(n).waveform = TXwaveform1;
    TX(n).Apod = Bmode_DW.TxWindow;
    TX(n).focus = Bmode_DW.txFocus;
    TX(n).Steer(1) = -Bmode_DW.Angles(n-k);
    TX(n).Delay = computeTXDelays(TX(n));

    TX(n+Bmode_DW.na) = TX(n);
    TX(n+Bmode_DW.na).waveform = TXwaveform2;
end

%Set event specific TX attributes for SW tracking and push
TX = CalculatePushTx(TX);

%% Receive parameters
wlsPer128 = 128/(4*2); % wavelengths in 128 samples for 4 samplesPerWave

%Bandpass filter
BandpassFilterOption = "fir1"; %S.InputBandpassFilter
%BandpassFilterOption = "Default";

Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
    'startDepth', 0, ...
    'endDepth', wlsPer128*ceil(Bmode_WB.maxAcqLength/wlsPer128), ... %wlsPer128*ceil(maxAcqLength/wlsPer128) or maxAcqLength
    'TGC', 1, ...
    'bufnum', 1, ...
    'framenum', 1, ...
    'acqNum', 1, ...
    'sampleMode', 'NS200BW', ... %NS200BW BS100BW BS67BW BS50BW
    'InputFilter',[], ...
    'demodFrequency',0, ...
    'mode', 0, ...
    'callMediaFunc', 0), 1, Resource.RcvBuffer(1).numFrames);

for i = 1:Resource.RcvBuffer(1).numFrames %loop over B-mode frames
        Receive(i).callMediaFunc = 1;
        Receive(i).framenum = i;
        Receive(i).acqNum = 1;
        Receive(i).startDepth = 0;
        Receive(i).endDepth = Bmode_WB.maxAcqLength;
        Receive(i).bufnum = 1;
        Receive(i).demodFrequency = FundamentalTransmitFrequency;
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [200,400,590,709,769,830,890,950];
TGC.rangeMax = Bmode_WB.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

%RcvProfile
RcvProfile.AntiAliasCutoff = 10;
RcvProfile.DCsubtract = 'on';
RcvProfile.PgaGain = 30; % overall gain in dB of the preamp output buffer prior to the A/D. 24 or 30 dB (default 24)
RcvProfile.LnaGain = 24; % overall gain in dB of the fixed-gain low noise input amplifier stage. 15, 18 or 24 dB (default 18).
RcvProfile.LnaZinSel = 31; % A value of 0 gives the lowest input impedance setting. Input impedance increases with increasing values, up to the highest impedance with a value of 31. Default 0

Recon = repmat(struct(...
    'senscutoff', 0.0, ...
    'pdatanum', 1, ...
    'rcvBufFrame',-1, ...
    'IntBufDest', [1,-1], ...
    'ImgBufDest', [1,-1], ...
    'RINums',1),1,1); 


%% Process parameters
% Specify Process structure array. (1) is used for B-mode imaging
% Specify Process structure array.
pers = 20;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
    'framenum',-1,...   % (-1 => lastFrame)
    'pdatanum',1,...    % number of PData structure to use
    'pgain',pgn,...            % pgain is image processing gain
    'reject',2,...      % reject level
    'persistMethod','simple',...
    'persistLevel',pers,...
    'interpMethod','4pt',...
    'grainRemoval','none',... %none medium
    'processMethod','none',... %none reduceSpeckle2
    'averageMethod','none',...
    'compressMethod','power',...
    'compressFactor',40,...
    'display',1,...      % display image after processing
    'displayWindow',1};

% EF1 is external function for UI control
Process(2).classname = 'External';
Process(2).method = 'UIControl';
Process(2).Parameters = {'srcbuffer','none'};

%% SeqControl and Events for shearwave generation
% - Change to Profile 1 (widebeams)
SeqControl(1).command = 'setTPCProfile';
SeqControl(1).condition = 'immediate';
SeqControl(1).argument = 1;

% - Change to Profile 5 (shear wave push)
SeqControl(2).command = 'setTPCProfile';
SeqControl(2).condition = 'immediate';
SeqControl(2).argument = 5;

% General PRF imaging
SeqControl(3).command = 'timeToNextAcq';
SeqControl(3).argument = 500; %2000 per second

%PRF shear wave push
SeqControl(4).command = 'timeToNextAcq';
SeqControl(4).argument = 5e5; %2 per second

% Return to Matlab
SeqControl(5).command = 'returnToMatlab';

nsc = length(SeqControl)+1;

% Specify Event structure arrays.
n = 1;

Event(n).info = 'ext func for UI control';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 2;
Event(n).seqControl = 0;
n = n+1;

%% Regular flash imaging starts from event(nStartFlash)
nStartFlash = n;
TTNA_SeqControl = 3;

%Get indices of transmit events at 0 degree steering. 
T = [TX(:).Steer];
T = T(1:2:end);
TX0 = find(T==0); %[widebeam pos, widebeam neg, focused pos, focused neg, diverging pos, diverging neg, tracking pos, tracking neg, push]
assert(length(TX0)==9,"Length of TX0 is not correct. Check if each imaging mode (widebeam, focused, diverging, tracking, push) has a transmit at 0 degree steer. ")

switch(TestedMode)
    case "WB"
        TPCprofileSeqControl = 1; 

        switch(TestedFreq)
            case "Fundamental"
                TX_index = TX0(1); 
                disp("Fundamental widebeam B-mode")
            case "Positive"
                TX_index = TX0(1); 
                disp("Harmonic widebeam B-mode positive polarity")
            case "Negative"
                TX_index = TX0(2);
                disp("Harmonic widebeam B-mode negative polarity")
            case "Coded"
                TX_index = TX0(1); 
                disp("Fundamental widebeam B-mode with coded excitation")
        end
    case "FC" 
        TPCprofileSeqControl = 1;
        switch(TestedFreq)
            case "Fundamental"
                TX_index = TX0(3); 
                disp("Fundamental focused B-mode")
            case "Positive"
                TX_index = TX0(3); 
                disp("Harmonic focused B-mode positive polarity")
            case "Negative"
                TX_index = TX0(4);
                disp("Harmonic focused B-mode negative polarity")
            case "Coded"
                TX_index = TX0(3); 
                disp("Fundamental focused B-mode with coded excitation")
        end
    case "DW"
        TPCprofileSeqControl = 1;
        switch(TestedFreq)
            case "Fundamental"
                TX_index = TX0(5); 
                disp("Fundamental diverging B-mode")
            case "Positive"
                TX_index = TX0(5); 
                disp("Harmonic diverging B-mode positive polarity")
            case "Negative"
                TX_index = TX0(6);
                disp("Harmonic diverging B-mode negative polarity")
            case "Coded"
                TX_index = TX0(5); 
                disp("Fundamental diverging B-mode with coded excitation")
        end
    case "SW_tracking"
        TPCprofileSeqControl = 1;
        switch(TestedFreq)
            case "Fundamental"
                TX_index = TX0(7);
                disp("Fundamental SW tracking")
            case "Positive"
                TX_index = TX0(7);
                disp("Harmonic SW tracking positive polarity")
            case "Negative"
                TX_index = TX0(8);
                disp("Harmonic SW tracking negative polarity")
            case "Coded"
                TX_index = TX0(7);
                disp("Fundamental SW tracking with coded excitation")
        end
    case "SW_push"
        TPCprofileSeqControl = 2;
        TTNA_SeqControl = 4;
        switch(TestedFreq)
            case "Fundamental"
                TX_index = TX0(9);
                disp("Fundamental SW push")
            case "Positive"
                error("Shear wave push should be fundamental")
            case "Negative"
                error("Shear wave push should be fundamental")
            case "Coded"
                error("Shear wave push should be fundamental")
        end
end

%Plot TXPD of selected transmit
figure
TXPD_event = computeTXPD(TX(TX_index),PData);
imagesc(TXPD_event(:,:,1))
axis image
axis off
title('TXPD of selected TX event')

MaxMI = 1.9; 
PlotMIvoltage(MaxMI,TX(TX_index),TW)

ReconInfo = repmat(struct('mode', 'replaceIntensity', ... %Accumulate IQ: accumIQ_normalize or accumIQ
    'Pre',[],...
    'Post',[],...
    'txnum', TX_index, ... %Number of transmit object used for acquisition
    'rcvnum', 1, ... %Number of receive object used for acquisition
    'scaleFactor', 1, ...
    'pagenum',1, ... %Indicates page no. for multi-page InterBuffer
    'threadSync', 1,...
    'regionnum', 1), 1, 1);


% Switch to TPC profile 1 (low power) for widebeam imaging
Event(n).info = 'Switch TPC profile.';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = TPCprofileSeqControl;
n = n+1;

nStartAcqFlash = n; 

for j = 1:10
    Event(n).info = 'Transmit pulse';
    Event(n).tx = TX_index;
    Event(n).rcv = j;
    Event(n).recon = 1;
    Event(n).process = 1;
    Event(n).seqControl = [TTNA_SeqControl,nsc];
    SeqControl(nsc).command = 'transferToHost';
    nsc = nsc + 1;
    n = n+1;
end

Event(n-1).seqControl = [Event(n-1).seqControl,5]; %Return to Matlab before jump

Event(n).info = 'Jump back';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = nsc;
SeqControl(nsc).command = 'jump';
SeqControl(nsc).argument = nStartAcqFlash;

%% User specified UI Control Elements
%Create UI
import vsv.seq.uicontrol.VsSliderControl;
import vsv.seq.uicontrol.VsButtonControl
[UIPos,SG,UI] = CreateUI(Resource);

EF(1).Function = vsv.seq.function.ExFunctionDef('UIControl', @UIControl);

% ROI position
ROIpos = [0-SW.ROI_width/2,TX(end).focus-SW.ROI_height/2,SW.ROI_width,SW.ROI_height];

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 1;

% Save all the structures to a .mat file.
save(['MatFiles/',filename]);

% Uncomment the following line to automatically run VSX every time you run
% this SetUp script (note that if VSX finds the variable 'filename' in the
% Matlab workspace, it will load and run that file instead of prompting the
% user for the file to be used):

%VSX

return

%% **** Callback routines used by UIControls (UI) ****

function PgainCallback(~,~,UIValue)
Process = evalin('base', 'Process');
Process(1).Parameters{find(strcmp(Process(1).Parameters,'pgain'))+1} = UIValue;
assignin('base','Process',Process);

Control.Parameters = {'Process'};
Control.Command = 'update&Run';
assignin('base','Control', Control);

end

function ObservePointCallback(~,~,UIState)

SW = evalin('base','SW');

%Observe left or right from focus point
if UIState == 1
    SW.ObserveLeftOfPush = true;
else
    SW.ObserveLeftOfPush = false;
end

assignin('base', 'SW', SW);

TX = evalin('base','TX');
TX = CalculatePushTx(TX);
assignin('base','TX',TX);

Control.Parameters = {'TX'};
Control.Command = 'update&Run';
assignin('base','Control', Control);

end

function SensCutoffCallback(~,~,UIValue)
%Sensitivity cutoff change
ReconL = evalin('base', 'Recon');
for i = 1:size(ReconL,2)
    ReconL(i).senscutoff = UIValue;
end
assignin('base','Recon',ReconL);
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'Recon'};
assignin('base','Control', Control);
end

function RangeChangeCallback(hObject,~,UIValue)

Trans = evalin('base','Trans'); %Needed for computeTGCWaveform (evalin('caller'))
TX = evalin('base','TX');

%Range change
simMode = evalin('base','Resource.Parameters.simulateMode');
% No range change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',evalin('base','Bmode_WB.endDepth'));
    return
end

%Update PData and Resource.DisplayWindow
Bmode_WB = evalin('base','Bmode_WB');
Bmode_FC = evalin('base','Bmode_FC');
Bmode_DW = evalin('base','Bmode_DW');
Bmode_WB.endDepth = UIValue;
assignin('base','Bmode_WB',Bmode_WB);

SW = evalin('base','SW');
[PData,~] = CalculatePData(Bmode_WB,Bmode_FC,Bmode_DW,SW,"update");
assignin('base','PData',PData);

%Update TXPD
for k = 1:length(TX)
    TX(k).TXPD = [];
    TX(k).TXPD = computeTXPD(TX(k),PData);
end

assignin('base','TX',TX);

%Update Receive
Receive = evalin('base', 'Receive');
maxAcqLength = ceil(Bmode_WB.endDepth/cos(Bmode_WB.theta/2));
for i = 1:size(Receive,2)
    Receive(i).endDepth = maxAcqLength;
end
assignin('base','Receive',Receive);

Bmode_WB.maxAcqLength = maxAcqLength;
assignin('base','Bmode_WB',Bmode_WB);

%Update TGC
TGC = evalin('base', 'TGC');
TGC.rangeMax = Bmode_WB.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);
assignin('base','TGC',TGC);

%Update Control to update&run
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','Receive','TGC','Recon','TX'};
assignin('base','Control', Control);
end

function NumPushElemCallback(~,~,UIValue)
%PushElements adjustment
TX = evalin('base','TX');
UI = evalin('base','UI');
freeze = evalin('base','freeze');
SW = evalin('base','SW');

pushElements = round(UIValue);

% If focusAdj is off, or at freeze, or Shearwave fig is closed, change value back to the old one
if freeze == 1
    set(UI(4).handle(2),'Value',TX(end).pushElements);
    set(UI(4).handle(3),'String',TX(end).pushElements);
else
    TX(end).oldElements = TX(end).pushElements;
    SW.nb_push_elmts = pushElements;
    assignin('base', 'SW', SW);

    TX = CalculatePushTx(TX);
    OverlaySW_UI();
    assignin('base','TX',TX);

    msg = ['Change push elements to ',num2str(pushElements),' elements...']

    Control.Parameters = {'TX'};
    Control.Command = 'update&Run';
    assignin('base','Control', Control);
end
end

function IQCaxisCallback(~,~,UIValue)
SW = evalin('base','SW');
SW.powermax = UIValue;
assignin('base', 'SW', SW);
end

function PushCyclesCallback(~,~,UIValue)
%SW.pushCycle change
TW = evalin('base', 'TW');
SW = evalin('base','SW');

% TW(3) is used for push
pushCycle = round(UIValue);
TW(3).oldCycles = evalin('base','SW.pushCycle');
TW(3).Parameters(3) = pushCycle*2;

SW.pushCycle = pushCycle;
assignin('base','SW',SW)

[~,~,~,~,TW(3)] = computeTWWaveform(TW(3));

msg = ['Change push cycles to ',num2str(pushCycle,'%4.0f'),' cycles...']

assignin('base','TW',TW);
Control.Parameters = {'TW'};
Control.Command = 'update&Run';
assignin('base','Control', Control);

end

function IQLoopFpsCallback(~,~,UIValue)
%fps adjustment for replay shearwave imaging
UI = evalin('base','UI');
loopfps = evalin('base','loopfps');
freeze = evalin('base','freeze');

if freeze
    assignin('base','loopfps', UIValue);
else
    set(UI(6).handle(2),'Value',loopfps);
    set(UI(6).handle(3),'String',loopfps);
end
end

function replayIQ(varargin)
freeze = evalin('base','freeze');
% figClose = evalin('base','figClose');

if freeze
    replay = 'on';
    assignin('base','replay','on');

    SW = evalin('base','SW');

    % need copyBuffers to access IQData
    Control.Command = 'copyBuffers';
    runAcq(Control); % NOTE:  If runacq() has an error, it reports it then exits MATLAB.

    UI = evalin('base','UI');
    TX = evalin('base','TX');
    MovieData = evalin('base','MovieData');
    SWIfigHandle = evalin('base','SWIfigHandle');

    set(UI(8).handle,'Visible','off');
    set(UI(9).handle,'Visible','on');

    disp('replay shearwave imaging....');

    IQMovie(SW.Ndetect-1) = struct('cdata',[],'colormap',[]);
    SWIimageHandle = evalin('base','SWIimageHandle');

    evalin('base','clear IQMovie');

    while strcmp(replay,'on')

        for i = 1:SW.Ndetect-1
            tic
            replay = evalin('base','replay');
            freeze = evalin('base','freeze');
            loopfps = evalin('base','loopfps');

            % if GUI is closed, return
            if evalin('base','vsExit == 1')
                return
            end

            % unfreeze will stop replay
            if(~freeze || strcmp(replay,'off') == 1)
                set(UI(8).handle,'Visible','on');
                set(UI(9).handle,'Visible','off');
                assignin('base','replay','off');
                return
            end

            set(SWIimageHandle,'CData',squeeze(MovieData(:,:,i)));
            drawnow
            %caxis(get(SWIfigHandle,'currentAxes'),[50,SW.powermax]);
            caxis(get(SWIfigHandle,'currentAxes'),[SW.PowerMin,SW.PowerMax]);

            frameData = getframe(SWIfigHandle);
            % padding to a multiple of two for MPEG-4 format
            [row,col,~]=size(frameData.cdata);
            if ~isequal(mod(row,2),0)
                frameData.cdata(end,:,:) = [];
            end
            if ~isequal(mod(col,2),0)
                frameData.cdata(:,end,:) = [];
            end
            IQMovie(i) = frameData;
            pause(1/loopfps-toc);
        end

        assignin('base','IQMovie',IQMovie);
    end
end

end

function saveIQ(varargin)
UI = evalin('base','UI');
freeze = evalin('base','freeze');

if evalin('base','exist(''IQMovie'',''var'')')
    if freeze
        assignin('base','replay','off');
        set(UI(8).handle,'Visible','on');
        set(UI(9).handle,'Visible','off');

        IQMovie = evalin('base','IQMovie');
        loopfps = evalin('base','loopfps');
        [fn,pn,filterindex] = uiputfile('*.mp4','Save Shearwave movie as');
        if ~isequal(fn,0) % fn will be zero if user hits cancel
            fn = strrep(fullfile(pn,fn), '''', '''''');

            vidObj = VideoWriter(fn,'MPEG-4');
            vidObj.Quality = 100;
            vidObj.FrameRate = loopfps;
            open(vidObj);
            writeVideo(vidObj,IQMovie);
            close(vidObj);

            fprintf('The shearwave movie has been saved at %s \n',fn);
        else
            disp('The shearwave movie is not saved.');
        end
    end
else
    msgbox('replay is not finished!');
end
end

function ROIWidthCallback(~,~,UIValue)
%ROI Width change
PData = evalin('base','PData');
ROIHandle = evalin('base','ROIHandle');
SWIfigHandle = evalin('base','SWIfigHandle');
SW = evalin('base','SW');

ROIpos = get(ROIHandle,'Position');

newX = ROIpos(1) - (UIValue-ROIpos(3))/2;

%TODO: Define edge cases
ROIpos(1) = newX;
ROIpos(3) = UIValue;
ROI_width = UIValue;
assignin('base','ROIpos',ROIpos);
SW.ROI_width = ROI_width;
assignin('base','SW',SW);

PData(2).Region(1).Shape.width = ROI_width;
PData(2).Region(1) = computeRegions(PData(2));

assignin('base','PData',PData);

Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','Recon'};
assignin('base','Control', Control);
assignin('base','newIQ',1);

pos = get(SWIfigHandle,'Position');
set(SWIfigHandle,'Position',[pos(1:2),[SW.ROI_width+3,SW.ROI_height]*7]);

end

function ROIHeigthCallback(~,~,UIValue)
%ROI Height change
Bmode_WB = evalin('base','Bmode_WB');
PData = evalin('base','PData');
ROIHandle = evalin('base','ROIHandle');
SWIfigHandle = evalin('base','SWIfigHandle');
SW = evalin('base','SW');

ROIpos = get(ROIHandle,'Position');

newY = ROIpos(2) - (UIValue-ROIpos(4))/2;
if newY < Bmode_WB.startDepth, newY = Bmode_WB.startDepth;
elseif newY > Bmode_WB.endDepth - UIValue, newY = Bmode_WB.endDepth - UIValue; end

ROIpos(2) = newY;
ROIpos(4) = UIValue;
ROI_height = UIValue;
assignin('base','ROIpos',ROIpos);

SW.ROI_height = ROI_height;
assignin('base','SW',SW);

PData(2).Region(1).Shape.Position(3) = newY;
PData(2).Region(1).Shape.height = ROI_height;
PData(2).Region(1) = computeRegions(PData(2));

assignin('base','PData',PData);

Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','Recon'};
assignin('base','Control', Control);
assignin('base','newIQ',1);

pos = get(SWIfigHandle,'Position');
set(SWIfigHandle,'Position',[pos(1:2),[SW.ROI_width+3,SW.ROI_height]*7]);
end



%% **** Callback routines used by External function definition (EF) ****

function UIControl(varargin)
% set WindowButtonDownFcn
Resource = evalin('base','Resource');
set(Resource.DisplayWindow(1).figureHandle,'WindowButtonDownFcn',@moveFocus);
set(Resource.DisplayWindow(1).figureHandle,'WindowButtonMotionFcn',[]);
assignin('base','Resource',Resource);

end

function [] = OverlaySW_UI()

% Focus mark, SWIROI,  and Transducer rect on bmode figure
persistent recHandle markHandle recTrans
Trans = evalin('base','Trans');
TX = evalin('base','TX');
SW = evalin('base','SW');
lambda_mm = evalin('base','lambda_mm');
x = TX(end).focusX;
z = TX(end).focus;
bmodeFigHandle = evalin('base','Resource.DisplayWindow(1).figureHandle');
if ishandle(bmodeFigHandle)
    ind_active1 = find(TX(end).Apod,1,"first");
    if  isempty(recHandle) || ~ishandle(recHandle)
        figure(bmodeFigHandle), hold on,
        recHandle = rectangle('Position',[x-SW.ROI_width/2,z-SW.ROI_height/2,SW.ROI_width,SW.ROI_height],'EdgeColor','w','LineWidth',2);
        %recTrans  = rectangle('Position',[x-TX(end).pushElements*Trans.spacing/2,4,TX(end).pushElements*Trans.spacing,4],'EdgeColor','r','FaceColor','r');
        recTrans  = rectangle('Position',[Trans.ElementPos(ind_active1,1)/lambda_mm,-4,TX(end).pushElements*Trans.spacing,4],'EdgeColor','r','FaceColor','r');
        markHandle = plot(x,z,'xr','MarkerFaceColor','r','MarkerSize',8,'Linewidth',2);hold off;
        assignin('base','TransHandle',recTrans);
        assignin('base','ROIHandle',recHandle);
        assignin('base','markHandle',markHandle);
    else
        set(markHandle,'XData',x);
        set(markHandle,'YData',z);
        %set(recTrans,'Position',[x-TX(end).pushElements*Trans.spacing/2,4,TX(end).pushElements*Trans.spacing,4]);
        set(recTrans,'Position',[Trans.ElementPos(ind_active1,1)/lambda_mm,-4,TX(end).pushElements*Trans.spacing,4]);
    end
end
end


function processIQ(IBuffer, QBuffer)
IQBuffer = complex(IBuffer, QBuffer); % VTS-1142 separate I and Q buffers
%processIQFunction: Computes power estimates from IQData
%		Im = I(k) * Q(k+1) - I(k+1) * Q(k)
%		Re = I(k) * I(k+1) + Q(k) * Q(k+1)
%		Power = sqrt(Im*Im + Re*Re);

disp("In ProcessIQ function")

ROIpos = evalin('base','ROIpos');
SWIfigHandle = evalin('base','SWIfigHandle');
SW = evalin('base','SW');
PData = evalin('base','PData');
lambda_mm = evalin('base','lambda_mm');
HarmonicImaging = evalin('base','SW.HarmonicImaging');
SeqControl = evalin('base','SeqControl');
ROILA = (PData(2).Region.PixelsLA+1);

Depth = floor(SW.ROI_height/PData(2).PDelta(3));%size(IQBuffer,1);%PData(2).Size(1);
Width = length(ROILA)/Depth;

if ~isequal(mod(Width,1),0)
    Width = floor(SW.ROI_width/PData(2).PDelta(1));%size(IQBuffer,2);%PData(2).Size(2);
    Depth = length(ROILA)/Width;
end

MovieData = zeros(Depth,Width,SW.Ndetect-1);

axisChannel = linspace(ROIpos(1),ROIpos(1)+ROIpos(3),Width);
axisDepth   = linspace(ROIpos(2),ROIpos(2)+ROIpos(4),Depth);

%Calculate power of IQ 3 and 4
Power = EstimatePower(IQBuffer(:,:,1,3),IQBuffer(:,:,1,4),ROILA,Depth,Width,SW.ObserveLeftOfPush);

%Set color limits
PowerMin = prctile(Power,10,"all");
PowerMax = prctile(Power,90,"all");
SW.PowerMin = PowerMin;
SW.PowerMax = PowerMax;
assignin('base','SW',SW)

%     %Estimate axial velocities
[m_mode,AxisSpatial,AxisTemporal] = EstimateAxialVelocities(IQBuffer,HarmonicImaging,SeqControl,PData,lambda_mm,SW.FocusZ_wavelength);
Q = any(m_mode,2);
ylim1 = AxisSpatial(find(Q,1,"first"))-5;
ylim2 = AxisSpatial(find(Q,1,"last"))+5;

figure(SWIfigHandle)
clf;
subplot(2,1,1)
myHandle = imagesc(axisChannel,axisDepth,Power);
caxis(get(SWIfigHandle,'currentAxes'),[SW.PowerMin,SW.PowerMax]);
title('Shear Wave Visualization');
xlabel('Wavelength');
ylabel('Wavelength');
colormap('gray')
axis image

subplot(2,1,2)
imagesc(AxisTemporal, AxisSpatial, m_mode);
xlabel('Time (ms)');
ylabel('Lateral (mm)');
colormap jet
clim([-5 5])
colorbar;
if(SW.DoAutomaticSWS)
    [Line_time,Line_spatial,SWS] = EstimateSWSradon(m_mode,AxisSpatial,AxisTemporal);
    hold on
    plot(Line_time,Line_spatial,'k','LineWidth', 2);
    title(sprintf('Shear Wave M-mode: Speed = %.2f m/s', SWS))
else
    title('Shear Wave M-mode' )
end
%ylim([ylim1,ylim2])

MovieData(:,:,1) = Power;
assignin('base','SWIimageHandle',myHandle);

%Loop over full ensemble
subplot(2,1,1)
%caxis(get(SWIfigHandle,'currentAxes'),[50,SW.powermax]);
caxis(get(SWIfigHandle,'currentAxes'),[SW.PowerMin,SW.PowerMax]);
% The size of IQData here is [nRows, nCols, nFrames, nPages]
for i = 2:SW.Ndetect-1 % for all combinations of 2 pages
    Power = EstimatePower(IQBuffer(:,:,1,i),IQBuffer(:,:,1,i+1),ROILA,Depth,Width,SW.ObserveLeftOfPush);
    set(myHandle,'CData',Power);
    MovieData(:,:,i) = Power;
    drawnow
end

assignin('base','MovieData',MovieData);

%When reprocessing data, save function workspace
if(evalin('base','ReprocessRcvData'))
    SWfilename = ['Luuk\SWdata_',datestr(now,'dd-mmmm-yyyy_HH-MM-SS'),'.mat'];
    save(SWfilename,'-v7.3'); %Save buffers (RF, IQ, image)
end
end

function Power = EstimatePower(IQ1,IQ2,ROILA,Depth,Width,ObserveLeftOfPush)
% kernel2D is used for 2D filter to smooth SWI
kernel2D = ...
    [   0.0073    0.0208    0.0294    0.0208    0.0073;
    0.0208    0.0589    0.0833    0.0589    0.0208;
    0.0294    0.0833    0.1179    0.0833    0.0294;
    0.0208    0.0589    0.0833    0.0589    0.0208;
    0.0073    0.0208    0.0294    0.0208    0.0073;];

ImMean = (imag(IQ1(ROILA)) + imag(IQ2(ROILA)))/2;
ReMean = (real(IQ1(ROILA)) + real(IQ2(ROILA)))/2;
Im = (imag(IQ1(ROILA))-ImMean) .* (real(IQ2(ROILA))-ReMean) - ...
    (imag(IQ2(ROILA))-ImMean) .* (real(IQ1(ROILA))-ReMean);
Re = (imag(IQ1(ROILA))-ImMean) .* (imag(IQ2(ROILA))-ImMean) + ...
    (real(IQ1(ROILA))-ReMean) .* (real(IQ2(ROILA))-ReMean);
Power = reshape((Im .* Im + Re .* Re).^0.125,Depth,Width);
buffer = filter2(Power,kernel2D,'full');
bufferS = size(buffer);
ind1 = (bufferS(1)-Depth)/2;
ind2 = (bufferS(2)-Width)/2;
Power = rot90(buffer(ind1+1:end-ind1,ind2+1:end-ind2),2);

%Modify ROI based on where we observe
if(ObserveLeftOfPush)
    Power(:,ceil(size(Power,2)/2:end)) = 0;
else
    Power(:,1:floor(size(Power,2)/2)) = 0;
end
end

function [m_mode,AxisSpatial,AxisTemporal] = EstimateAxialVelocities(IQ,HarmonicImaging,SeqControl,PData,lambda_mm,FocusZ_wavelength)
SW = evalin('base','SW');
fs = 1e6/SeqControl(11).argument; % Frame rate in Hz
if(HarmonicImaging)
    fs = fs/2;
end
dx = PData(2).PDelta(1)*lambda_mm; % Spatial resolution in mm (e.g., 0.5 mm)

%Clutter filter IQ data
IQ = squeeze(IQ);
if(isfield(SW,'ClutterFilter'))
    switch(SW.ClutterFilter)
        case "None"
            IQ = IQ;
        case "Poly"
            IQ = wfilt(IQ,'poly',SW.ClutterFilterN);
        case "SVD"
            IQ = wfilt(IQ,'svd',SW.ClutterFilterN);
    end
end

%Velocity estimation
IQ1 = IQ(:,:,1:end-1);
IQ2 = IQ(:,:,2:end);

phase_diff = angle(IQ2 .* conj(IQ1));
velocity = (fs / (4 * pi)) * phase_diff; % Approximate axial velocity

%Bandpass filter velocities. This removes constant and fast-changing
%velocities, effectively smoothing the velocity plot.
%[b, a] = butter(3, [5 50]/(fs/2), 'bandpass');
%velocity_filtered = filtfilt(b, a, velocity);
velocity_filtered = velocity;

%Modify ROI based on where we observe
if(SW.ObserveLeftOfPush)
    velocity_filtered(:,ceil(size(velocity_filtered,2)/2:end),:) = 0;
else
    velocity_filtered(:,1:floor(size(velocity_filtered,2)/2),:) = 0;
end

if(evalin('base','exist("M_lineX","var")'))
    disp("Using anatomical M-line from workspace")
    x = evalin('base','M_lineX');
    z = evalin('base','M_lineZ');
    [ind, ~] = GetMLine([z(1),x(1)],[z(2),x(2)],size(velocity,1:2));
    ind = ind(:);
    m_mode = velocity_filtered(ind + (0:size(velocity_filtered,3)-1).*(size(velocity_filtered,1)*size(velocity_filtered,2)));
else
    M_line = round((FocusZ_wavelength-PData(2).Origin(3))/PData(2).PDelta(3)); % Index of M-mode line (horizontal line through focal depth of push)
    m_mode = squeeze(velocity_filtered(M_line,:, :));
end

m_mode(:,1) = 0; %Remove first column
m_mode = medfilt2(m_mode,[20,3]); %Apply 2D median filter

AxisSpatial = (1:size(m_mode,1))*dx;
AxisTemporal = 1000.*(1:size(m_mode,2))/fs;

end

function [Line_time,Line_spatial,SWS] = EstimateSWSradon(m_mode,AxisSpatial,AxisTemporal)

%Set first columns to 0 (too high intensity)
m_mode(:,1:3) = 0;

%Prepare input data as struct
data = MakeDataStruct(AxisSpatial, AxisTemporal, m_mode);

%Apply Radon transform
theta = CalcTheta(data.dxdt);
radout = NormRadon(data.data, theta);

% Find Peak
peak = FindRadonPeaks(radout);

% % Calculate Trajectory
out = CalcTrajectory(data, peak);
%res = CalcResolution(data, radout, peak);

Line_time = out.times;
Line_spatial = data.xMm;
SWS = out.speed;

end

%% Other external functions used in callback functions

function moveFocus(varargin)

freeze = evalin('base','freeze');
vsExit = evalin('base','vsExit');

if freeze == 0 && vsExit == 0  % no response if freeze or exit
    Trans = evalin('base','Trans');
    TX = evalin('base','TX');
    lambda_mm = evalin('base','lambda_mm');

    %Get new push location in wavelengths and cm
    bmodeFigHandle = evalin('base','Resource.DisplayWindow(1).figureHandle');
    bmodeAxes = get(bmodeFigHandle,'currentAxes');
    currentPos = get(bmodeAxes,'CurrentPoint');

    SW = evalin('base','SW');

    SW.FocusX_wavelength = currentPos(1);
    SW.FocusX_cm = SW.FocusX_wavelength*1540/3.125e6*100;
    SW.FocusZ_wavelength = currentPos(3);
    SW.FocusZ_cm = SW.FocusZ_wavelength*1540/3.125e6*100;

    assignin('base','SW', SW);

    TX = CalculatePushTx(TX);
    OverlaySW_UI();

    %Update ROI
    PData = evalin('base','PData');
    PData(2).Region.Shape = struct(...
        'Name','Rectangle',...
        'Position',[SW.FocusX_wavelength,0,SW.FocusZ_wavelength-SW.ROI_height/2],...
        'width', SW.ROI_width,...
        'height', SW.ROI_height);

    PData(2).Region = computeRegions(PData(2));

    assignin('base','TX', TX);
    assignin('base','PData', PData);

    Control.Parameters = {'TX','PData','Recon'};
    Control.Command = 'update&Run';
    assignin('base','Control', Control);

    markHandle = evalin('base','markHandle');
    TransHandle = evalin('base','TransHandle');
    ROIHandle = evalin('base','ROIHandle');

    set(markHandle,'XData',SW.FocusX_wavelength);
    set(markHandle,'YData',SW.FocusZ_wavelength);
    set(TransHandle,'Position',[Trans.ElementPos(find(TX(end).Apod,1))/lambda_mm,4,TX(end).pushElements*Trans.spacing,4]);  %[x,y,width,height]
    set(ROIHandle,'Position',[SW.FocusX_wavelength-SW.ROI_width/2,SW.FocusZ_wavelength-SW.ROI_height/2,SW.ROI_width,SW.ROI_height]);  %[x,y,width,height]

    assignin('base','markHandle', markHandle);
    assignin('base','TransHandle', TransHandle);
    assignin('base','ROIHandle', ROIHandle);

end
end

function closeIQfig(varargin)

if evalin('base','vsExit')
    delete(gcf)
else
    SWIoff();
end

end

function SWIoff(varargin)

% if the movie is replay, no response
replay = evalin('base','replay');

if strcmp(replay,'on')
    msgbox('please stop reply first');

elseif strcmp(replay,'off')
    % change the start event
    nStart = evalin('base','nStartFlash');
    Control = evalin('base','Control');
    if isempty(Control(1).Command), n=1; else n=length(Control)+1; end
    Control(n).Command = 'set&Run';
    Control(n).Parameters = {'Parameters',1,'startEvent',nStart};
    evalin('base',['Resource.Parameters.startEvent =',num2str(nStart),';']);
    assignin('base','Control',Control);
end

end

function SWIon(varargin)

Resource = evalin('base','Resource');
SW = evalin('base','SW');
UI = evalin('base','UI');

% Handle for Shearwave figure
if ~evalin('base','exist(''SWIfigHandle'',''var'')')

    SWIfigHandle = figure('Name','ShearWaveVisulization',...
        'NumberTitle','off','Visible','on',...
        'Position',[495 177 552 403], ...
        'CloseRequestFcn',@closeIQfig);

    assignin('base','SWIfigHandle',SWIfigHandle);
else
    set(evalin('base','SWIfigHandle'),'Visible','on');
end

% newIQ is used for check whether new SWI axes is required
evalin('base','newIQ = 1;');

% change the start event
nStart = evalin('base','nStartAlternativeBmode');
Control = evalin('base','Control');
if isempty(Control(1).Command), n=1; else n=length(Control)+1; end
evalin('base',['Resource.Parameters.startEvent =',num2str(nStart),';']); %Update start event in Resource object
Control(n).Command = 'set&Run';
Control(n).Parameters = {'Parameters',1,'startEvent',nStart}; %Updates Resource.Parameters.startEvent
assignin('base','Control',Control);

%Make save button visible
set(UI(6).handle,'Visible','on');

end

function TX = CalculatePushTx(TX)

Trans = evalin('base','Trans');
PData = evalin('base','PData');
SW = evalin('base','SW');
Bmode_WB = evalin('base','Bmode_WB');
lambda_mm = evalin('base','lambda_mm');

%Observe
n = length(TX)-2;
if(SW.ObserveLeftOfPush)
    ObserveLocationX_wavelength = 0;
else
    ObserveLocationX_wavelength = 0;
end

ObserveAngle = atan(ObserveLocationX_wavelength/(SW.dApex+SW.FocusZ_wavelength));
%ObserveAngle = 0;
%warning('Using 0 angle in push observe')
TXorg = SW.dApex*tan(ObserveAngle);

[SW,TXwaveform1,TXwaveform2] = SelectTXwaveforms(SW);
assignin('base','SW',SW)
TX(n).waveform = TXwaveform1;
TX(n).Origin = [TXorg, 0.0, 0.0];

%Transmit using all elements
TX(n).Apod = Bmode_WB.TxWindow;
TX(n).Steer = [ObserveAngle,0.0];
TX(n).Delay = computeTXDelays(TX(n));
TX(n).TXPD = computeTXPD(TX(n),PData);

TX(n+1).waveform = TXwaveform2;
TX(n+1).Origin = TX(n).Origin;
TX(n+1).Apod = TX(n).Apod;
TX(n+1).Steer = TX(n).Steer;
TX(n+1).Delay = TX(n).Delay;
TX(n+1).TXPD = TX(n).TXPD;

%Push
n = length(TX);
TX(n).waveform = 3;
TX(n).focus = SW.FocusZ_wavelength;       % wavelength, can be changed in the GUI
TX(n).focusX = SW.FocusX_wavelength;      % can be changed in the GUI
TX(n).pushElements = SW.nb_push_elmts;       % can be changed in the GUI
TX(n).FocalPt = [SW.FocusX_wavelength,0,SW.FocusZ_wavelength];

%Select NTxBmode closest elements
d = sqrt((TX(n).FocalPt(1)-Trans.ElementPos(:,1)./lambda_mm).^2+(TX(n).FocalPt(3)).^2);
[~,ElementIndex] = sort(d);
TX(n).Apod = zeros(1,Trans.numelements);
TX(n).Apod(sort(ElementIndex(1:SW.nb_push_elmts))) = 1;

%Apply windowing
ce=ElementIndex(1);
lft = find(TX(n).Apod,1,"first");
rt = find(TX(n).Apod,1,"last");
switch(SW.PushApodization)
    case "kaiser"
        disp("Push transmit uses Kaiser window")
        TX(n).Apod(lft:rt) = kaiser(length(lft:rt),3)';
    case "none"
        disp("Push transmit uses rectangular window")
end

%Compute delays and TXPD
TX(n).Delay = computeTXDelays(TX(n));
TX(n).TXPD = computeTXPD(TX(n),PData);

end

function [PData,DisplayWindow] = CalculatePData(Bmode_WB,Bmode_FC,Bmode_DW,SW,mode)
% Set up PData structure. The number of Regions matches the number of rays, but the Regions
% are larger in angle.  The apex of the SectorFT region is positioned above the virtual apex of the
% scan, to obtain a wider sector at the transducer face.
PData(1).Coord = 'rectangular';
PData(1).PDelta = [0.5, 0, 0.5];
PData(1).Size(1) = 10 + ceil((Bmode_WB.endDepth-Bmode_WB.startDepth)/PData(1).PDelta(3));
PData(1).Size(2) = 10 + ceil(2*(Bmode_WB.endDepth + Bmode_WB.dApex)*sin(Bmode_WB.theta/2)/PData(1).PDelta(1));
PData(1).Size(3) = 1;
PData(1).Origin = [-(PData(1).Size(2)/2)*PData(1).PDelta(1),0,Bmode_WB.startDepth];
PData(1).Region = repmat(struct(...
    'Shape',struct('Name','SectorFT', ...
    'Position',[0,0,-Bmode_WB.dApex], ...
    'z',Bmode_WB.startDepth, ...
    'r',Bmode_WB.dApex+Bmode_WB.endDepth, ...
    'angle',Bmode_WB.RegionAngle, ... %Bmode_WB.rayDelta*5.5
    'steer',0, ...
    'andWithPrev',1)),1,Bmode_WB.na+1);

% First Region covers entire field of view
PData(1).Region(1).Shape.andWithPrev = 0;
PData(1).Region(1).Shape.angle = Bmode_WB.theta;

% Set ray line Regions Position and steering angles
Bmode_WB.Angles = (-Bmode_WB.theta/2):Bmode_WB.rayDelta:(Bmode_WB.theta/2);
dz = Bmode_WB.dApex/2;  % distance from virtual apex to SectorFT apex.
for n = 2:Bmode_WB.na+1
    PData(1).Region(n).Position = [-dz*tan(Bmode_WB.Angles(n-1)),0,(Bmode_WB.dApex+dz)];
    PData(1).Region(n).Shape.steer = Bmode_WB.Angles(n-1);
end
PData(1).Region = computeRegions(PData(1));
[Bmode_WB.x,Bmode_WB.z] = PDataToAxis(PData(1));
assignin('base','Bmode_WB',Bmode_WB)

% Specify PData(2) structure array for Shearwave visulization
PData(2) = PData(1);
PData(2).Region = PData(2).Region(1);

PData(2).Region.Shape = struct(...
    'Name','Rectangle',...
    'Position',[SW.FocusX_wavelength,0,SW.FocusZ_wavelength-SW.ROI_height/2],...
    'width', SW.ROI_width,...
    'height', SW.ROI_height);

PData(2).Region = computeRegions(PData(2));
[SW.x,SW.z] = PDataToAxis(PData(2));
assignin('base','SW',SW)

% Specify rectangular PData structure array with a SectorFT region for each ray.
PData(3).Coord = 'rectangular';
PData(3).PDelta = [0.5, 0, 0.5];  % x, y, z pdeltas
PData(3).Size(1) = 10 + ceil((Bmode_FC.endDepth-Bmode_FC.startDepth)/PData(3).PDelta(3));
PData(3).Size(2) = 10 + ceil(2*(Bmode_FC.endDepth + Bmode_FC.radius)*sin(-Bmode_FC.theta)/PData(3).PDelta(1));
PData(3).Size(3) = 1;      % single image page
PData(3).Origin = [-(PData(3).Size(2)/2)*PData(3).PDelta(1),0,Bmode_FC.startDepth];
% - specify 128 Region structures.
PData(3).Region = repmat(struct('Shape',struct( ...
    'Name','SectorFT',...
    'Position',[0,0,-Bmode_FC.radius],...
    'z',Bmode_FC.startDepth,...
    'r',Bmode_FC.radius+Bmode_FC.endDepth,...
    'angle',Bmode_FC.rayDelta,...
    'steer',0)),1,Bmode_FC.na);
% - set position of regions to correspond to beam spacing.
for i = 1:Bmode_FC.na
    PData(3).Region(i).Shape.steer(1) = Bmode_FC.theta + (i-1)*Bmode_FC.rayDelta;
end
PData(3).Region = computeRegions(PData(3));
[Bmode_FC.x,Bmode_FC.z] = PDataToAxis(PData(3));
assignin('base','Bmode_FC',Bmode_FC)

% Specify PData(4) structure array for DW
PData(4).Coord = 'rectangular';
PData(4).PDelta = [0.5,0,0.5];
PData(4).Size(1) = 10+ceil((Bmode_DW.endDepth-Bmode_DW.startDepth)/PData(4).PDelta(3)); % startDepth, endDepth and pdelta set PData(4).Size.
PData(4).Size(2) = 10 + ceil (2*(Bmode_DW.endDepth + Bmode_DW.radius) * sin(-Bmode_DW.theta)/PData(4).PDelta(1));
PData(4).Size(3) = 1;      % single image page
PData(4).Origin = [-PData(4).Size(2)/2*PData(4).PDelta(1),0,Bmode_DW.startDepth]; % x,y,z of upper lft crnr. (wls)

PData(4).Region = struct(...
    'Shape',struct('Name','SectorFT', ...
    'Position',[0,0,-Bmode_DW.radius], ...
    'z',Bmode_DW.startDepth, ...
    'r',Bmode_DW.radius+Bmode_DW.endDepth, ...
    'angle',-Bmode_DW.fullAngle, ...
    'steer',0));
PData(4).Region = computeRegions(PData(4));
[Bmode_DW.x,Bmode_DW.z] = PDataToAxis(PData(4));
assignin('base','Bmode_DW',Bmode_DW)

PData(5) = PData(1);
PData(6) = PData(1);

if(mode == "update")
    DisplayWindow = [];
    return;
end

DisplayWindow(1).Title = 'Widebeam_imaging';
DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/DisplayWindow(1).pdelta);
DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
    DwWidth, DwHeight];
DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
DisplayWindow(1).Type = 'Matlab'; %Verasonics or Matlab
DisplayWindow(1).numFrames = Bmode_WB.Nframes;
DisplayWindow(1).AxesUnits = 'mm'; %mm or wavelengths
DisplayWindow.Colormap = gray(256);
end

function ReceivePulse = ComputeReceivePulse(TW,S)

%For non-harmonic excitations, the receive pulse is simply the time-reversed
%transmit pulse
if(~S.HarmonicImaging)
    ReceivePulse = TransmitToReceivePulse(TW,S.demodFrequency);
else
    %For harmonic excitations, the receive pulse depends on the excitation
    %type
    switch(evalin('base','CodedExcitationType'))
        case "Chirp"
            %For a chirp, the matched filter is a time-reversed harmonic
            %chirp
            TW_harmonic = TW;
            TW_harmonic.envFrequency = linspace(2*TW.envFrequency(1),2*TW.envFrequency(end),TW.envNumCycles);
            [~, ~, ~, ~, TW_harmonic] = computeTWWaveform(TW_harmonic);
            ReceivePulse = TransmitToReceivePulse(TW_harmonic,S.demodFrequency);
        case "Barker"
            %For a Barker pulse, the matched filter is matched to the envelope shape
            TW_harmonic = TW;
            TW_harmonic.envFrequency = 2.*TW_harmonic.envFrequency;
            TW_harmonic.envPulseWidth = abs(TW_harmonic.envPulseWidth);
            [~, ~, ~, ~, TW_harmonic] = computeTWWaveform(TW_harmonic);
            ReceivePulse = TransmitToReceivePulse(TW_harmonic,S.demodFrequency);
        otherwise
            error("Undefined coded excitation type")
    end
end
end

function ReceivePulse = TransmitToReceivePulse(TW,demodFrequency)
%Transmitted pulse, fs = 250 MHz
TransmitPulse = TW.Wvfm2Wy;
ts_TransmitPulse = 1/250e6;
t_TransmitPulse = 0:ts_TransmitPulse:(length(TransmitPulse)-1)*ts_TransmitPulse;

%Received pulse, fs = 4*fdemod
ts_ReceivePulse = 1/(2e6*demodFrequency);
t_ReceivePulse = 0:ts_ReceivePulse:t_TransmitPulse(end);
ReceivePulse = interp1(t_TransmitPulse,TransmitPulse,t_ReceivePulse,"spline");
ReceivePulse = fliplr(ReceivePulse);
ReceivePulse = ReceivePulse';
end

function RF = decodeRF(RF_in,ReceivePulse)
RF = conv2(double(RF_in), ReceivePulse, 'same')./length(ReceivePulse);
RF = circshift(RF,-round(length(ReceivePulse)/2),1);
RF = int16(RF);
end

function RF = decodeRF_WB(RF_in)
ReceivePulse = evalin('base','Bmode_WB.ReceivePulse');
RF = decodeRF(RF_in,ReceivePulse);
end

function RF = decodeRF_SW(RF_in)
ReceivePulse = evalin('base','SW.ReceivePulse');
RF = decodeRF(RF_in,ReceivePulse);
end

function RF = decodeRF_FC(RF_in)
ReceivePulse = evalin('base','Bmode_FC.ReceivePulse');
RF = decodeRF(RF_in,ReceivePulse);
end

function RF = decodeRF_DW(RF_in)
ReceivePulse = evalin('base','Bmode_DW.ReceivePulse');
RF = decodeRF(RF_in,ReceivePulse);
end

function TX = ReplaceTXPD(TX,TW_index,PData)
TX_tmp = TX;
TX_tmp.waveform = TW_index;
TX_tmp.TXPD = [];
TX_tmp.TXPD = computeTXPD(TX_tmp,PData);
TX.TXPD = TX_tmp.TXPD;
end

function [x,z] = PDataToAxis(PData)
lambda_mm = evalin('base','lambda_mm');
x = PData.Origin(1):PData.PDelta(1):PData.Origin(1)+(PData.Size(2)-1)*PData.PDelta(1);
x = x.*lambda_mm;
z = PData.Origin(3):PData.PDelta(3):PData.Origin(3)+(PData.Size(1)-1)*PData.PDelta(3);
z = z.*lambda_mm;
end

function [S,TXwaveform1,TXwaveform2] = SelectTXwaveforms(S)

if(~isfield(S,'HarmonicOverlap'))
    S.HarmonicOverlap = false;
end

if(S.HarmonicImaging | S.HarmonicOverlap)
    S.demodFrequency = 2*evalin('base','HarmonicTransmitFrequency');
else
    S.demodFrequency = evalin('base','FundamentalTransmitFrequency');
end

bF1 = fir1(40,2e6*[S.demodFrequency-1.5, S.demodFrequency+1.5]./(4e6*S.demodFrequency),'bandpass');
S.InputBandpassFilter = bF1(1:21);

if(S.CodedExcitation)
    if(S.HarmonicImaging | S.HarmonicOverlap)
        TXwaveform1 = 6;
        TXwaveform2 = 7;
    else
        TXwaveform1 = 2;
        TXwaveform2 = 2;
    end
else
    if(S.HarmonicImaging | S.HarmonicOverlap)
        TXwaveform1 = 4;
        TXwaveform2 = 5;
    else
        TXwaveform1 = 1;
        TXwaveform2 = 1;
    end
end

if(S.CodedExcitation)
    TW = evalin('base','TW');
    S.ReceivePulse = ComputeReceivePulse(TW(TXwaveform1),S);
else
    S.ReceivePulse = [];
end

end

function [] = SetStartEventToPush()
nStart = evalin('base','nStartPushSeq');
Control = evalin('base','Control');
if isempty(Control(1).Command), n=1; else n=length(Control)+1; end
Control(n).Command = 'set&Run';
Control(n).Parameters = {'Parameters',1,'startEvent',nStart};
evalin('base',['Resource.Parameters.startEvent =',num2str(nStart),';']);
assignin('base','Control',Control);
end

function [UIPos,SG,UI] = CreateUI(Resource)
import vsv.seq.uicontrol.VsSliderControl;
import vsv.seq.uicontrol.VsButtonControl

% Define UIPos, which contains the default GUI positions - three columns of 10 controls. The x,y
%    locations increment up columns, with each column being a separate page. The origin
%    specified by UIPos is the lower left corner of a virtual box that encloses the control.
UIPos = zeros(10,2,3);
UIPos(:,1,1) = 0.0625;
UIPos(:,1,2) = 0.375;
UIPos(:,1,3) = 0.6875;
UIPos(:,2,1) = 0.0:0.1:0.9;
UIPos(:,2,2) = 0.0:0.1:0.9;
UIPos(:,2,3) = 0.0:0.1:0.9;

% Define slider group offsets and sizes. All units are normalized.
SG = struct('TO',[0.0,0.0975],...   % title offset
    'TS',[0.25,0.025],...   % title size
    'TF',0.8,...            % title font size
    'SO',[0.0,0.06],...     % slider offset
    'SS',[0.25,0.031],...   % slider size
    'EO',[0.075,0.031],...   % edit box offset
    'ES',[0.11,0.031]);     % edit box size

% - Push Elements Adjustment
TX = evalin('base', 'TX');
TX(end).oldElements = evalin('base','TX(end).pushElements');
assignin('base','TX',TX);
UI(4).Control = VsSliderControl('LocationCode','UserB4','Label','Push Elements',...
    'SliderMinMaxVal',[10,80,evalin('base','TX(end).pushElements')],...
    'SliderStep',[1/70,5/70],'ValueFormat','%3.0f',...
    'Callback', @NumPushElemCallback );

replay = 'off';
loopfps = 5;
assignin('base','replay',replay);
assignin('base','loopfps',loopfps);

UI(15).Control = VsSliderControl('LocationCode','UserB3','Label','Push cycles',...
    'SliderMinMaxVal',[50,3500,evalin('base','SW.pushCycle')],...
    'SliderStep',[10/1500,50/1500],'ValueFormat','%3.0f',...
    'Callback', @PushCyclesCallback );

%PData slider
UI(16).Control = VsSliderControl('LocationCode','UserC8','Label','Processing gain',...
    'SliderMinMaxVal',[1,21,evalin('base','pgn')],...
    'SliderStep',[1/20,4/20],'ValueFormat','%.2f',...
    'Callback', @PgainCallback);

end

function [] = ProcessTurnSWoff()
SWIoff();
end

function [] = Start3beatsAcq(varargin)
% change the start event
nStart = evalin('base','nStartAcqStrain');
Control = evalin('base','Control');
if isempty(Control(1).Command), n=1; else n=length(Control)+1; end
evalin('base',['Resource.Parameters.startEvent =',num2str(nStart),';']); %Update start event in Resource object
Control(n).Command = 'set&Run';
Control(n).Parameters = {'Parameters',1,'startEvent',nStart}; %Updates Resource.Parameters.startEvent
assignin('base','Control',Control);
end

function [] = SaveSW_RF(varargin)
if ~isempty(findobj('tag','UI')) % running VSX
    if evalin('base','freeze')==0   % no action if not in freeze
        msgbox('Please freeze VSX');
        return
    else
        Control.Command = 'copyBuffers';
        runAcq(Control); % NOTE:  If runAcq() has an error, it reports it then exits MATLAB.
    end
else % not running VSX
    if evalin('base','exist(''RcvData'',''var'');')
        RcvData = evalin('base','RcvData');
    else
        disp('RcvData does not exist!');
        return
    end
end

RFfilename = ['Luuk\RFdata_',datestr(now,'dd-mmmm-yyyy_HH-MM-SS'),'.mat'];
assignin('base','RFfilename',RFfilename)
%save(RFfilename,'-v7.3'); %Save buffers (RF, IQ, image)
% evalin('base','save(RFfilename,"-append")'); %Append settings from base workspace

evalin('base','save(RFfilename,"-v7.3")'); %Save settings from base workspace

%Resource exists in base and this function, but contains different fields.
%We want to merge these two structures.
Resource_func = Resource;
Resource = evalin('base','Resource');

mergestructs = @(x,y) cell2struct([struct2cell(x);struct2cell(y)],[fieldnames(x);fieldnames(y)]);

Resource.InterBuffer = mergestructs(Resource.InterBuffer,Resource_func.InterBuffer);
Resource.ImageBuffer = mergestructs(Resource.ImageBuffer,Resource_func.ImageBuffer);
Resource.DisplayWindow = mergestructs(Resource.DisplayWindow,Resource_func.DisplayWindow);
if(isfield(Resource_func,'RcvBuffer'))
    Resource.RcvBuffer = mergestructs(Resource.RcvBuffer,Resource_func.RcvBuffer);
end

clear Resource_func
save(RFfilename,'-append'); %Append buffers to mat file

fprintf('The RF data has been saved at %s \n',RFfilename);
end

function [] = SaveStrainData()
%Save IQ data
Control.Command = 'copyBuffers';
runAcq(Control);
RFfilename = ['Luuk\Bmode_strain_',datestr(now,'dd-mmmm-yyyy_HH-MM-SS'),'.mat'];

assignin('base','RFfilename',RFfilename)
%save(RFfilename,'-v7.3'); %Save buffers (RF, IQ, image)
% evalin('base','save(RFfilename,"-append")'); %Append settings from base workspace

evalin('base','save(RFfilename,"-v7.3")'); %Save settings from base workspace

%Resource exists in base and this function, but contains different fields.
%We want to merge these two structures.
Resource_func = Resource;
Resource = evalin('base','Resource');

mergestructs = @(x,y) cell2struct([struct2cell(x);struct2cell(y)],[fieldnames(x);fieldnames(y)]);

Resource.InterBuffer = mergestructs(Resource.InterBuffer,Resource_func.InterBuffer);
Resource.ImageBuffer = mergestructs(Resource.ImageBuffer,Resource_func.ImageBuffer);
Resource.DisplayWindow = mergestructs(Resource.DisplayWindow,Resource_func.DisplayWindow);
Resource.RcvBuffer = mergestructs(Resource.RcvBuffer,Resource_func.RcvBuffer);

clear Resource_func
save(RFfilename,'-append'); %Append buffers to mat file

fprintf('The RF data has been saved at %s \n',RFfilename);
% change the start event
nStart = evalin('base','nStartFlash');
Control = evalin('base','Control');
if isempty(Control(1).Command), n=1; else n=length(Control)+1; end
evalin('base',['Resource.Parameters.startEvent =',num2str(nStart),';']); %Update start event in Resource object
Control(n).Command = 'set&Run';
%Control(n).Command = 'set';
Control(n).Parameters = {'Parameters',1,'startEvent',nStart}; %Updates Resource.Parameters.startEvent
assignin('base','Control',Control);
end

function Event = ModEventObject(Event)
%Add recon and decode events
for k = 1:length(Event)
    if(isequal(Event(k).info,"recon FC"))
        Event(k).recon = 3;
    end

    if(isequal(Event(k).info,"recon DW"))
        Event(k).recon = 4;
    end

    if(isequal(Event(k).info,"recon and process for SWI"))
        Event(k).recon = 2;
        Event(k).process = 3;
    end

    if(isequal(Event(k).info,"recon") | isequal(Event(k).info,"recon B-mode SW")) %B-mode during SW
        if(Event(k).recon>0) %This is a regular B-mode event
            continue;
        else
            Event(k).info = 'recon B-mode SW';
            Event(k).recon = 5;
        end
    end

    if(isequal(Event(k).info,"recon RF strain")) %B-mode strain
        Event(k).recon = 6;
    end

    if(isequal(Event(k).info,"Decode pulse RF FC"))
        Event(k).process = 6;
    end

    if(isequal(Event(k).info,"Decode pulse RF DW"))
        Event(k).process = 7;
    end

    if(isequal(Event(k).info,"Decode pulse RF SW"))
        if(Event(k).process ==5) %Already decoded RF data
            Event(k).process = 0;
        else
            Event(k).process = 5;
        end
    end

    if(isequal(Event(k).info,"Decode pulse RF") | isequal(Event(k).info,"Decode pulse RF B-mode SW"))
        Event(k).info = 'Decode pulse RF B-mode SW';
        if(Event(k).process ==10) %Already decoded RF data
            Event(k).process = 0;
        else
            Event(k).process = 10;
        end
    end

    if(isequal(Event(k).info,"Decode pulse RF strain"))
        Event(k).process = 12;
    end

    if(isequal(Event(k).info,"Display image SW"))
        Event(k).process = 1;
    end

end

end

function [] = PlotMIvoltage(MI,TX_sel,TW)
%Given a transmit frequency, compute the minimum voltage over depth that
%corresponds to a MI of 1.9

%Transmit frequency [MHz] 
ind_waveform = TX_sel.waveform; 

if(TW(ind_waveform).type == "parametric")
    ft = TW(ind_waveform).Parameters(1);
else
    ft = mean(TW(ind_waveform).envFrequency); 
end

%Derate by 0.3 dB/cm/MHz
Depth = linspace(0,10,101); %Depth axis in cm
DerateFactor_dB = 0.3*Depth*ft; %dB
DerateFactor = db2mag(DerateFactor_dB);

%Calculate voltage that corresponds to 1.9 MI
Pr = MI*sqrt(ft); 
sens_total = GetHydrophoneSensitivity(ft);
V = Pr*sens_total;
V = V.*DerateFactor;

figure
plot(Depth,V,LineWidth=2,DisplayName="Curve")   
grid minor
xlabel('Depth (cm)')
ylabel('V-min limit (V)')
title('Allowed voltage to achieve MI=1.9')
hold on 
xline(0,DisplayName=sprintf("%.2f V at d=0 cm",GetVoltageAtDepth(Depth,V,0)))
xline(4,DisplayName=sprintf("%.2f V at d=4 cm",GetVoltageAtDepth(Depth,V,4)))
xline(6,DisplayName=sprintf("%.2f V at d=6 cm",GetVoltageAtDepth(Depth,V,6)))
xline(8,DisplayName=sprintf("%.2f V at d=8 cm",GetVoltageAtDepth(Depth,V,8)))
xline(10,DisplayName=sprintf("%.2f V at d=10 cm",GetVoltageAtDepth(Depth,V,10)))
legend('show',Location='northwest')
end

function V_depth = GetVoltageAtDepth(Depth,V,d)
    V_depth = V(find(Depth>=d,1));
end

function sens_total = GetHydrophoneSensitivity(ft)
%Returns the pressure in MPa that corresponds to the measured voltage.
%Based on calibration data June 2024.

persistent sens_hydrophone capacitance_hydrophone sens_preamp capacitance_preamp

%Read parameters as function of frequency
if(isempty(sens_hydrophone))
    load('HydrophoneTables.mat','AmplifierTable','HydrophoneTable')
    
    ind_hydrophone = find(HydrophoneTable.FREQ_MHz >= ft,1);
    sens_hydrophone = db2mag(HydrophoneTable{ind_hydrophone,"SENS_DB"})*1e6;
    capacitance_hydrophone = HydrophoneTable{ind_hydrophone,"CAP_PF"}*1e-12;
    
    ind_amplifier = find(AmplifierTable.FREQ_MHZ >= ft,1);
    sens_preamp = db2mag(AmplifierTable{ind_amplifier,"GAIN_DB"});
    capacitance_preamp = AmplifierTable{ind_amplifier,"CAP_PF"}*1e-12;
end

%Average/constant values
% sens_hydrophone = db2mag(-254)*1e6; %%From calibration: -254 dB re. 1V/microPa
% capacitance_hydrophone = 1.2314e-11; %From calibration
% sens_preamp = db2mag(20); %From calibration
% capacitance_preamp = 6.10e-12; %From calibration

sens_total = sens_preamp*sens_hydrophone*capacitance_hydrophone/(capacitance_preamp+capacitance_hydrophone);
sens_total = sens_total*1e6; %[V/MPa]
end


function sens_total = GetHydrophoneSensitivity_old()
sens_hydrophone = 2.034e-7; %This should be -254 dB re. 1V/microPa
capacitance_hydrophone = 1.2314e-11; %Where does this come from? Should be 30e-12 according to the datasheet
sens_preamp = db2mag(20); %From calibration
capacitance_preamp = 6.10e-12; %From calibration

sens_total = sens_preamp*sens_hydrophone*capacitance_hydrophone/(capacitance_preamp+capacitance_hydrophone);
sens_total = sens_total*1e6; %[V/MPa]

end