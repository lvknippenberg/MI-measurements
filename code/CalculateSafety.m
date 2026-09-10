function [MI,I_sppa,I_spta] = CalculateSafety(Depth,ft,Vmin_measured,TW,rho,c,PRF,MI_method, FN, ScopeImpedance_Ohm, withPreamp, LoadCap_pF)

%Input arguments:
% - depth at which Vmin was measured (cm)
% - Transmit frequency (MHz)
% - Measured Vmin (V)
% - TW structure of transmitted pulse (used to simulate the transmitted
% pulse)
% - tissue density rho (kg/m3), typically 1000 kg/m3
% - speed of sound c (m/s), typically 1540
% - PRF of pulse
% - MI_method: "ParallelCircuit" (correct; applies the Onda open-circuit ->
%   loaded capacitive-divider correction, Eq. 2a of Onda_HydroCalMethod) or
%   "Simple". Defaults to "ParallelCircuit" if empty/omitted.
% - Filename used for saving
% - ScopeImpedance_Ohm: input impedance of the oscilloscope used to record
%   Vmin (50 or 1e6). See load-correction note below. Defaults to 1e6.
%   (only relevant WITH the preamp — its 50 ohm output stage is what the
%   resistive load correction accounts for.)
% - withPreamp: true (default) = AH-2010 preamp in the chain -> apply the
%   preamp gain, its 6.1 pF input-capacitance divider, and the scope-impedance
%   load correction. false = hydrophone connected DIRECTLY to the scope
%   (preamp REMOVED, e.g. to avoid its ~2 V-peak clipping at high transmit
%   voltage): NO gain and NO resistive load correction; instead the bare
%   capacitive hydrophone is loaded by the full cable+scope capacitance.
% - LoadCap_pF: total cable+scope capacitance [pF] seen by the hydrophone
%   when withPreamp=false. This DOMINATES the reading (a bare HGL-0400 is only
%   ~12 pF) and MUST be measured/known; there is no sensible default.

MakePlots = false;

%Default to the physically correct calibration combination
if(nargin < 8 || isempty(MI_method))
    MI_method = "ParallelCircuit";
end

%Filename for saving is optional
if(nargin < 9)
    FN = [];
end

%Preamp in the chain? (default yes)
if(nargin < 11 || isempty(withPreamp))
    withPreamp = true;
end
if(nargin < 12)
    LoadCap_pF = [];
end

if(Vmin_measured==0)
    MI = 0;
    I_sppa = 0;
    I_spta = 0;
    return
end

% --- Correct for oscilloscope input impedance ------------------------------
% The Onda AH-2010 preamp has a 50 ohm OUTPUT impedance, and its 20 dB gain
% (hence the Onda combined system sensitivity M_L) is referenced to a 50 ohm
% load. A 50 ohm source driving a high-impedance scope (1 MOhm) delivers ~2x
% the voltage it delivers into 50 ohm, so a high-Z reading must be scaled back
% to its 50-ohm-equivalent before applying the calibration:
%     V_50equiv = V_meas * [Rref/(Rref+Rout)] / [R_L/(R_L+Rout)]
% -> factor 0.5 for R_L = 1 MOhm, 1.0 for R_L = 50 ohm.
% (Only applies WITH the preamp: the resistive divider is between the preamp's
%  50 ohm output and the scope input. Without the preamp there is no 50 ohm
%  output stage, so no resistive correction — the capacitive load divider is
%  handled in VoltageToPressure instead.)
if(withPreamp)
    if(nargin < 10 || isempty(ScopeImpedance_Ohm))
        ScopeImpedance_Ohm = 1e6;   % default: high-Z scope input
        warning('CalculateSafety:scopeImp', ...
            ['ScopeImpedance_Ohm not supplied; assuming 1 MOhm (high-Z). ' ...
             'Pass 50 if the scope was terminated into 50 ohm.']);
    end
    Rout = 50;   % AH-2010 output impedance [ohm]
    Rref = 50;   % load the amplifier gain / calibration is referenced to [ohm]
    LoadCorr = (Rref/(Rref+Rout)) / (ScopeImpedance_Ohm/(ScopeImpedance_Ohm+Rout));
    Vmin_measured = Vmin_measured * LoadCorr;
end

%Derate by 0.3 dB/cm/MHz
DerateFactor_dB = 0.3*Depth*ft; %dB
DerateFactor = db2mag(DerateFactor_dB);
Vmin_measured = Vmin_measured/DerateFactor;

%Calculate MI (max 1.9)
Pr = VoltageToPressure(Vmin_measured,ft,MI_method,withPreamp,LoadCap_pF); %[MPa]
MI = Pr/sqrt(ft);

%Simulate transmit pulse
if(isfield(TW,"Wvfm2Wy"))
    TransmitPulse_sim = TW.Wvfm2Wy;
else
    [TransmitPulse_sim,~,~,~] = computeTWWaveform(TW(1));
end
dt = 1/250;
t = 0:dt:(length(TransmitPulse_sim)-1)*dt; %[us]

%Rescale pulse to match measured Vmin
Vmin_sim = min(TransmitPulse_sim);
TransmitPulse_voltage = TransmitPulse_sim*(-Vmin_measured/Vmin_sim); %[V]

%Plot transmit pulse
if(MakePlots)
    figure
    plot(t,TransmitPulse_sim)
    hold on
    plot(t,TransmitPulse_voltage)
    xlabel('Time (us)')
    ylabel('Voltage (V)')
    title('(Rescaled) transmit pulse')
    grid minor
    legend('Simulated','Rescaled')
end

if(nargout==1)
    return;
end

%Convert transmit pulse to pressure and intensity
TransmitPulse_pressure = VoltageToPressure(TransmitPulse_voltage,ft,MI_method,withPreamp,LoadCap_pF); %[MPa]
TransmitPulse_pressure = 1e6*TransmitPulse_pressure; %[Pa]
TransmitPulse_intensity = TransmitPulse_pressure.^2/(rho*c); %[W/m2]

%Calculate pulse duration based on relative intensity integral
RelativeIntensityIntegral = rescale(cumsum(TransmitPulse_intensity),0,1);
RelativeIntensityIntegral_10p = find(RelativeIntensityIntegral>0.1,1);
RelativeIntensityIntegral_90p = find(RelativeIntensityIntegral>0.9,1);
PulseDuration_us = 1.25*(t(RelativeIntensityIntegral_90p)-t(RelativeIntensityIntegral_10p));

%Plot relative intensity integral
if(MakePlots)
    figure
    plot(t,RelativeIntensityIntegral)
    hold on
    xline(t(RelativeIntensityIntegral_10p),'--')
    xline(t(RelativeIntensityIntegral_90p),'--')
    yline(0.1,'--')
    yline(0.9,'--')
    title(sprintf("Pulse duration %.2f us",PulseDuration_us))
    xlabel('Time (us)')
    ylabel('Relative intensity integral')
    grid minor

    figure
    plot(TransmitPulse_intensity)
    hold on
    xline(RelativeIntensityIntegral_10p)
    xline(RelativeIntensityIntegral_90p)
    grid minor
    title('PII calculation')
end

%Calculate pulse intensity interval (PII)
PII = trapz(TransmitPulse_intensity(RelativeIntensityIntegral_10p:RelativeIntensityIntegral_90p))*(dt/1e6); %J/m2

%Spatial peak pulse average intensity (I_sppa), max 190 W/cm2
I_sppa = PII/(PulseDuration_us/1e6); %[W/m2]
I_sppa = I_sppa/10000; %[W/cm2]

%Spatial peak time average intensity (I_spta), max 430 mW/cm2
I_spta = PII*PRF; %[W/m2]
I_spta = 1000*I_spta; %[mW/m2]
I_spta = I_spta/10000; %[mW/cm2]

if(~isempty(FN))
    save(FN)
end

end

function P = VoltageToPressure(V,ft,MI_method,withPreamp,LoadCap_pF)
%Returns the pressure in MPa that corresponds to the measured voltage.
%Based on calibration data June 2024.

if(nargin < 4 || isempty(withPreamp)), withPreamp = true; end
if(nargin < 5), LoadCap_pF = []; end

persistent sens_hydrophone capacitance_hydrophone sens_preamp capacitance_preamp

%Read parameters as function of frequency
if(isempty(sens_hydrophone))
    load('HydrophoneTables.mat','AmplifierTable','HydrophoneTable')

    ind_hydrophone = find(HydrophoneTable.FREQ_MHz >= ft,1);
    sens_hydrophone = db2mag(HydrophoneTable{ind_hydrophone,"SENS_DB"})*1e6; %V/Pa
    capacitance_hydrophone = HydrophoneTable{ind_hydrophone,"CAP_PF"}*1e-12; %Capicitance in F

    ind_amplifier = find(AmplifierTable.FREQ_MHZ >= ft,1);
    sens_preamp = db2mag(AmplifierTable{ind_amplifier,"GAIN_DB"}); %Gain factor
    capacitance_preamp = AmplifierTable{ind_amplifier,"CAP_PF"}*1e-12;
end

%Average/constant values
% sens_hydrophone = db2mag(-254)*1e6; %%From calibration: -254 dB re. 1V/microPa
% capacitance_hydrophone = 1.2314e-11; %From calibration
% sens_preamp = db2mag(20); %From calibration
% capacitance_preamp = 6.10e-12; %From calibration

if(~withPreamp)
    % Preamp REMOVED: hydrophone straight into the scope. No gain. The bare
    % hydrophone (C_h ~ 12 pF, open-circuit sensitivity M_c) is loaded by the
    % full cable+scope capacitance C_load, so the voltage divides by
    % C_h/(C_h + C_load). C_load dominates and must be supplied.
    if(isempty(LoadCap_pF))
        warning('CalculateSafety:noPreampLoadCap', ...
            ['withPreamp=false but LoadCap_pF (cable+scope) not supplied. The ' ...
             'load capacitance dominates a bare-hydrophone reading and MUST be ' ...
             'measured. Falling back to open-circuit (no divider) — this will ' ...
             'strongly UNDER-read the true pressure.']);
        C_load = 0;
    else
        C_load = LoadCap_pF*1e-12;
    end
    sens_total = sens_hydrophone*capacitance_hydrophone/(capacitance_hydrophone+C_load); %[V/Pa]
    sens_total = sens_total*1e6; %[V/MPa]
    P = V./sens_total; %[MPa]
    return
end

switch(MI_method)
    case "ParallelCircuit"
        sens_total = sens_preamp*sens_hydrophone*capacitance_hydrophone/(capacitance_preamp+capacitance_hydrophone);
        sens_total = sens_total*1e6; %[V/MPa]
        P = V./sens_total; %[MPa]

    case "Simple"
        V = V/10; %Gain from pre-amplifier
        sens_total = sens_hydrophone*1e6; %Only use sensitivity of hydrophone [V/MPa]
        P = V./sens_total; %[MPa]
end


end
