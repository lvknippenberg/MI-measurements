function [sens, info] = SystemSensitivity(f_MHz, opts)
%SYSTEMSENSITIVITY  Loaded sensitivity [V/MPa] of the hydrophone measuring chain.
%
%   sens = SystemSensitivity(f_MHz) returns the sensitivity of the default
%   chain -- Onda HGL-0400 into the AH-2010-025 pre-amplifier, read on a
%   50 ohm-terminated scope -- at the acoustic working frequency f_MHz, so that
%
%       p [MPa] = V_measured [V] / sens
%
%   Name-value options
%   ------------------
%     Chain        "preamp"   (default) AH-2010 in the chain
%                  "nopreamp" hydrophone cable straight into the scope
%     ScopeImpedance_Ohm   50 (default) or 1e6. Only meaningful for "preamp":
%                  the AH-2010 has a 50 ohm OUTPUT impedance and its gain -- so
%                  the Onda combined sensitivity -- is referenced to a 50 ohm
%                  load. A high-Z scope therefore reads ~2x too much, and the
%                  returned sensitivity is doubled to absorb that.
%     LoadCap_pF   total cable+scope capacitance seen by the bare hydrophone.
%                  Required for "nopreamp"; it dominates the reading (the bare
%                  HGL-0400 is only ~12 pF) and must be measured, not guessed.
%                  See CalibrateLoadCapacitance.
%
%   Physics
%   -------
%   The Onda certificate is an OPEN-CIRCUIT sensitivity M_c(f) (the file says
%   ELECTRICAL_LOAD OpenCircuit), so Onda's loaded-sensitivity formula
%   (Onda_HydroCalMethod Eq. 2a) requires the capacitive divider term:
%
%       with preamp:   sens = G(f) * M_c(f) * C_H/(C_H + C_A)
%       no preamp:     sens =        M_c(f) * C_H/(C_H + C_load)
%
%   Dropping the divider ("Simple" method) under-reads the pressure by 1.5x.
%
%   See also: SafetyIndices, AcousticWorkingFrequency, miCalibration.

    arguments
        f_MHz (1,1) double {mustBePositive}
        opts.Chain (1,1) string {mustBeMember(opts.Chain,["preamp","nopreamp"])} = "preamp"
        opts.ScopeImpedance_Ohm (1,1) double = 50
        opts.LoadCap_pF double = []
    end

    C  = miCalibration();
    HT = C.HydrophoneTable;  AT = C.AmplifierTable;

    Moc = 10.^(interp1(HT.FREQ_MHz, HT.SENS_DB, f_MHz)/20) * 1e6;   % V/Pa, open circuit
    Ch  =       interp1(HT.FREQ_MHz, HT.CAP_PF,  f_MHz) * 1e-12;    % F

    switch opts.Chain
        case "preamp"
            gain = 10.^(interp1(AT.FREQ_MHZ, AT.GAIN_DB, f_MHz)/20);
            Ca   =       interp1(AT.FREQ_MHZ, AT.CAP_PF,  f_MHz) * 1e-12;
            sens = gain * Moc * Ch/(Ch + Ca) * 1e6;                 % V/MPa @ 50 ohm

            % A 50 ohm source into R_L delivers R_L/(R_L+50) of its open-circuit
            % swing; the calibration is referenced to R_L = 50 ohm.
            Rout = 50;  RL = opts.ScopeImpedance_Ohm;
            sens = sens * (RL/(RL+Rout)) / (Rout/(Rout+Rout));
            info = struct('Moc_V_per_Pa',Moc,'Ch_pF',Ch*1e12,'gain',gain, ...
                          'Ca_pF',Ca*1e12,'LoadCorr',(RL/(RL+Rout))/0.5);

        case "nopreamp"
            if isempty(opts.LoadCap_pF)
                error('SystemSensitivity:noLoadCap', ...
                    ['Chain="nopreamp" needs LoadCap_pF (cable+scope capacitance). ' ...
                     'It dominates a bare-hydrophone reading and must be measured -- ' ...
                     'see CalibrateLoadCapacitance.']);
            end
            Cl   = opts.LoadCap_pF * 1e-12;
            sens = Moc * Ch/(Ch + Cl) * 1e6;                        % V/MPa
            info = struct('Moc_V_per_Pa',Moc,'Ch_pF',Ch*1e12,'gain',1, ...
                          'Ca_pF',opts.LoadCap_pF,'LoadCorr',1);
    end
    info.f_MHz = f_MHz;  info.chain = opts.Chain;  info.sens_V_per_MPa = sens;
end
