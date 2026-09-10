function C = miCalibration()
%MICALIBRATION  The Onda calibration tables, cached.
%
%   C = miCalibration() returns a struct with
%       C.HydrophoneTable   FREQ_MHz, SENS_DB (open circuit), CAP_PF
%       C.AmplifierTable    FREQ_MHZ, GAIN_DB, PHASE_DEG, CAP_PF
%       C.file              where it came from
%
%   Source of truth are the two Onda certificates in calibration/; the .mat is
%   the digested form of exactly those files (see calibration/README.md).

    persistent C_
    if isempty(C_)
        f  = miRoot('calibration','HydrophoneTables.mat');
        S  = load(f,'HydrophoneTable','AmplifierTable');
        C_ = struct('HydrophoneTable',S.HydrophoneTable, ...
                    'AmplifierTable', S.AmplifierTable, 'file', f);
    end
    C = C_;
end
