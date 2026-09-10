function [Cload_pF, info] = CalibrateLoadCapacitance(dataDir, patPreamp, patNoPreamp, V, f_MHz)
%CALIBRATELOADCAPACITANCE  Cable+scope capacitance seen by the bare hydrophone.
%
%   Cload_pF = CalibrateLoadCapacitance(dataDir, patPreamp, patNoPreamp, V, f_MHz)
%   uses paired captures of the SAME acoustic point taken with and without the
%   pre-amplifier over the transmit voltages V (chosen below the preamp's clip
%   ceiling, where both chains are linear). The ratio of the two readings fixes
%   the only unknown in the no-preamp divider:
%
%       V_pre / V_no = [G * C_H/(C_H+C_A)] / [C_H/(C_H+C_load)]
%   =>  C_load = <ratio> * C_H/(G * C_H/(C_H+C_A)) - C_H
%
%   patPreamp / patNoPreamp are sprintf patterns taking the transmit voltage,
%   e.g. 'S5-1_%dV_Preamp_peak.hws' and 'S5-1_%dV_NoPreamp_opt.hws'.
%
%   This value is a property of the CABLE + ADAPTERS + SCOPE INPUT only. It has
%   to be redone whenever any of those change, and it is what makes a
%   no-preamp measurement quantitative at all.

    C  = miCalibration();
    Moc = 10.^(interp1(C.HydrophoneTable.FREQ_MHz, C.HydrophoneTable.SENS_DB, f_MHz)/20)*1e6;
    Ch  =       interp1(C.HydrophoneTable.FREQ_MHz, C.HydrophoneTable.CAP_PF,  f_MHz)*1e-12;
    sensPre = SystemSensitivity(f_MHz, Chain="preamp", ScopeImpedance_Ohm=50)/1e6;  % V/Pa

    r = nan(size(V));
    for i = 1:numel(V)
        fp = fullfile(dataDir, sprintf(patPreamp,   V(i)));
        fn = fullfile(dataDir, sprintf(patNoPreamp, V(i)));
        if isfile(fp) && isfile(fn), r(i) = peakNeg(fp)/peakNeg(fn); end
    end
    Cload_pF = (mean(r,'omitnan')*(Moc*Ch)/sensPre - Ch)*1e12;
    info = struct('V',V,'ratio',r,'meanRatio',mean(r,'omitnan'), ...
                  'spread_pct',100*std(r,'omitnan')/mean(r,'omitnan'),'f_MHz',f_MHz);
end

function v = peakNeg(f)
    w = readHWS(f);  y = w(1).y - median(w(1).y);  v = -min(y);
end
