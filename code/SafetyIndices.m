function S = SafetyIndices(y, dt, opts)
%SAFETYINDICES  MI, I_sppa.3 and I_spta.3 from one measured hydrophone waveform.
%
%   S = SafetyIndices(y, dt, Name=Value) converts a hydrophone voltage record
%   y [V] sampled at interval dt [s] into the derated acoustic-output indices.
%   This is the stand-alone core of the repository: no Verasonics, no
%   toolboxes, no assumed pulse shape -- everything comes from the capture.
%
%   S = SafetyIndices(hwsFile, [], Name=Value) reads the capture first.
%
%   Required options
%   ----------------
%     Sensitivity   chain sensitivity [V/MPa] at f_awf (see SystemSensitivity)
%     Depth_cm      distance from the transducer face to the measurement point
%     Fawf_MHz      acoustic working frequency (see AcousticWorkingFrequency)
%
%   Optional
%   --------
%     PRF_Hz        pulse repetition frequency used for I_spta.3. For a burst
%                   protocol use the EFFECTIVE rate (pushes per burst divided
%                   by burst + idle time), not the in-burst rate. Default NaN
%                   -> I_spta.3 is returned as NaN.
%     Rho           tissue density [kg/m^3], default 1000
%     C             speed of sound  [m/s],   default 1500
%     DerateCoef    derating coefficient [dB/cm/MHz], default 0.3 (FDA/IEC)
%     RemoveBaseline  subtract the median of y first, default true
%
%   Returned fields
%   ---------------
%     pr_MPa        peak-rarefactional pressure IN WATER
%     pr3_MPa       ... derated
%     MI            = pr3 / sqrt(f_awf)
%     PD_s          pulse duration, 1.25 x (t90 - t10) of the intensity integral
%     PII_J_per_m2  pulse intensity integral, in water
%     Isppa3_W_cm2  derated spatial-peak pulse-average intensity
%     Ispta3_mW_cm2 derated spatial-peak time-average intensity
%     derate        amplitude derating factor (>1); intensities use its square
%
%   Definitions follow FDA Track 3 / IEC 62359: derating 0.3 dB/cm/MHz,
%   MI = p_r.3/sqrt(f_awf), PD from the 10-90 % points of the cumulative
%   intensity integral, I_sppa.3 = PII.3/PD, I_spta.3 = PII.3 x PRF.
%   Limits for reference: MI 1.9, I_sppa.3 190 W/cm^2, I_spta.3 720 mW/cm^2.
%
%   See also: SystemSensitivity, AcousticWorkingFrequency, readHWS, RunDemo.

    arguments
        y
        dt = []
        opts.Sensitivity (1,1) double {mustBePositive}
        opts.Depth_cm    (1,1) double {mustBeNonnegative}
        opts.Fawf_MHz    (1,1) double {mustBePositive}
        opts.PRF_Hz      (1,1) double = NaN
        opts.Rho         (1,1) double = 1000
        opts.C           (1,1) double = 1500
        opts.DerateCoef  (1,1) double = 0.3
        opts.RemoveBaseline (1,1) logical = true
    end

    if (ischar(y) || isstring(y)) && isempty(dt)
        w  = readHWS(char(y));
        dt = w(1).dt;
        y  = w(1).y;
    end
    y = y(:);
    if opts.RemoveBaseline, y = y - median(y); end

    % --- pressure ---------------------------------------------------------
    p = y / opts.Sensitivity * 1e6;                       % Pa
    derate = 10.^(opts.DerateCoef*opts.Depth_cm*opts.Fawf_MHz/20);   % amplitude, >1

    S.pr_MPa  = -min(p)/1e6;
    S.pr3_MPa = S.pr_MPa / derate;
    S.MI      = S.pr3_MPa / sqrt(opts.Fawf_MHz);

    % --- pulse duration and PII (IEC 62359 10-90 % rule) ------------------
    I  = p.^2 / (opts.Rho*opts.C);                        % W/m^2, in water
    cR = rescale(cumsum(I), 0, 1);
    i10 = find(cR > 0.1, 1);
    i90 = find(cR > 0.9, 1);
    S.PD_s         = 1.25*(i90 - i10)*dt;
    S.PII_J_per_m2 = trapz(I(i10:i90))*dt;

    Ider = (1/derate)^2;                                  % intensity derating
    S.Isppa3_W_cm2  = S.PII_J_per_m2 * Ider / S.PD_s / 1e4;
    S.Ispta3_mW_cm2 = S.Isppa3_W_cm2 * S.PD_s * opts.PRF_Hz * 1e3;

    % --- provenance -------------------------------------------------------
    S.derate   = derate;
    S.i10      = i10;   S.i90 = i90;
    S.dt       = dt;
    S.Fawf_MHz = opts.Fawf_MHz;
    S.Depth_cm = opts.Depth_cm;
    S.PRF_Hz   = opts.PRF_Hz;
    S.Sensitivity_V_per_MPa = opts.Sensitivity;
end
