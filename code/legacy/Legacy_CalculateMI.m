close all;clear all;clc;

cd("D:\Luuk van Knippenberg\Onedrive folder\OneDrive - TU Eindhoven\PhD\Automatic MI\IntensityMeasurements")

%Constants
rho = 1000; %Density [kg/m3]
c = 1500; %Speed of sound [m/s]

%% B-mode 4 cm
load('BmodeTransmitPulse_4cm.mat')
ft = TW(1).Parameters(1); %[MHz]
Vmin_measured = 4.26; %[V]
Depth = 4; %cm

%Get pulse repetition frequency
PRI = SeqControl(1).argument; %[us]
PRF = 1e6/PRI; %[Hz]
%PRF = 10;

FN = "SafetyBmode4cm.mat";
CalculateSafety(Depth,ft,Vmin_measured,TW,rho,c,PRF,FN)
SafetyBmode4cm = load(FN);

%% B-mode 6 cm
load('BmodeTransmitPulse_6cm.mat')
ft = TW(1).Parameters(1); %[MHz]
Vmin_measured = 4.19; %[V]
Depth = 6; %cm

%Get pulse repetition frequency
PRI = SeqControl(1).argument; %[us]
PRF = 1e6/PRI; %[Hz]

FN = "SafetyBmode6cm.mat";
CalculateSafety(Depth,ft,Vmin_measured,TW,rho,c,PRF,FN)
SafetyBmode6cm = load(FN);

%% B-mode 8 cm
load('BmodeTransmitPulse_8cm.mat')
ft = TW(1).Parameters(1); %[MHz]
Vmin_measured = 3.53; %[V]
Depth = 8; %cm

%Get pulse repetition frequency
PRI = SeqControl(1).argument; %[us]
PRF = 1e6/PRI; %[Hz]

FN = "SafetyBmode8cm.mat";
CalculateSafety(Depth,ft,Vmin_measured,TW,rho,c,PRF,FN)
SafetyBmode8cm = load(FN);

%% Shear wave push
load('ShearWaveTransmitPulse.mat')
ft = TW(2).Parameters(1); %[MHz]
Vmin_measured = 3.19; %[V]
Depth = 4; %cm

%Get pulse repetition frequency
PRF = 5; 

FN = "SafetySW4cm.mat";
CalculateSafety(Depth,ft,Vmin_measured,TW(2),rho,c,PRF,FN)
SafetySW4cm = load(FN);

%% Shear wave push (reduce PRF to 7000 Hz to meet FDA requirement on I_spta)
load('ShearWaveTransmitPulse.mat')
ft = TW(1).Parameters(1); %[MHz]
Vmin_measured = 3.13; %[V]
Depth = 0; %cm

%Get pulse repetition frequency
% PRI = SeqControl(6).argument; %[us]
% PRF = 1e6/PRI; %[Hz]
PRF = 7000; 

FN = "SafetyDW.mat";
CalculateSafety(Depth,ft,Vmin_measured,TW(1),rho,c,PRF,FN)
SafetyDW = load(FN);


