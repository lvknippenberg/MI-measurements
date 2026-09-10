function [] = CalculateSafety(Depth,ft,Vmin_measured,TW,rho,c,PRF,FN)
%Derate by 0.3 dB/cm/MHz
DerateFactor_dB = 0.3*Depth*ft %dB
DerateFactor = db2mag(DerateFactor_dB)
Vmin_measured = Vmin_measured/DerateFactor

%Calculate MI
Pr = VoltageToPressure(Vmin_measured) %[MHz]
MI = Pr/sqrt(ft)

%Simulate transmit pulse
[TransmitPulse_sim,~,~,~] = computeTWWaveform(TW(1));
dt = 1/250;
t = 0:dt:(length(TransmitPulse_sim)-1)*dt; %[us]

%Rescale pulse to match measured Vmin
Vmin_sim = min(TransmitPulse_sim);
TransmitPulse_voltage = TransmitPulse_sim*(-Vmin_measured/Vmin_sim); %[V]

%Plot transmit pulse
figure
plot(t,TransmitPulse_sim)
hold on
plot(t,TransmitPulse_voltage)
xlabel('Time (us)')
ylabel('Voltage (V)')
title('(Rescaled) transmit pulse')
grid minor
legend('Simulated','Rescaled')

%Convert transmit pulse to pressure and intensity
TransmitPulse_pressure = VoltageToPressure(TransmitPulse_voltage); %[MPa]
TransmitPulse_pressure = 1e6*TransmitPulse_pressure; %[Pa]
TransmitPulse_intensity = TransmitPulse_pressure.^2/(rho*c);

%Calculate pulse duration based on relative intensity integral
RelativeIntensityIntegral = rescale(cumsum(TransmitPulse_intensity),0,1);
RelativeIntensityIntegral_10p = find(RelativeIntensityIntegral>0.1,1);
RelativeIntensityIntegral_90p = find(RelativeIntensityIntegral>0.9,1);
PulseDuration_us = 1.25*(t(RelativeIntensityIntegral_90p)-t(RelativeIntensityIntegral_10p))

%Plot relative intensity integral
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

%Calculate pulse intensity interval (PII)
PII = trapz(TransmitPulse_intensity(RelativeIntensityIntegral_10p:RelativeIntensityIntegral_90p))*(dt/1e6) %J/m2

%Spatial peak pulse average intensity (I_sppa)
I_sppa = PII/(PulseDuration_us/1e6); %[W/m2]
I_sppa = I_sppa/10000 %[W/cm2]

%Spatial peak time average intensity (I_spta)
I_spta = PII*PRF; %[W/m2]
I_spta = 1000*I_spta; %[mW/m2]
I_spta = I_spta/10000 %[mW/cm2]

save(FN)

end

function P = VoltageToPressure(V)
%Returns the pressure in MPa that corresponds to the measured voltage.
%Based on the Excel sheet provided by Xufei.

sens_hydrophone = 2.034e-7; %This should be -254 dB re. 1V/microPa
capacitance_hydrophone = 1.2314e-11; %Where does this come from? Should be 30e-12 according to the datasheet
sens_preamp = db2mag(20); %From calibration
capacitance_preamp = 6.10e-12; %From calibration

sens_total = sens_preamp*sens_hydrophone*capacitance_hydrophone/(capacitance_preamp+capacitance_hydrophone);
sens_total = sens_total*1e6; %[V/MPa]

P = V./sens_total; %[MPa]
end
