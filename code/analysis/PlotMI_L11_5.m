function PlotMI_L11_5(dataDir, options)
%PLOTMI_L11_5 MI vs TX voltage for the L11-5 pulse-inversion sequence.
%   PlotMI_L11_5(dataDir) reads the hydrophone captures in dataDir
%   (normal/inverted pulse x 50/1e6 ohm scope impedance x TX-voltage sweep),
%   converts peak-negative voltage to MI with the same calibration path as
%   CalculateSafety.m (Onda parallel-circuit sensitivity + scope-impedance
%   load correction + 0.3 dB/cm/MHz derating at the measured depth), and
%   compares against the Verasonics-reported MI.
%
%   The hydrophone/preamp path SATURATES at high drive (a hard ceiling in the
%   peak voltage), so the true MI at high TX voltage is obtained by LINEAR
%   EXTRAPOLATION of the unsaturated (low/mid voltage) region, not from the
%   saturated readings.
%
%   Default dataDir is the "L11-5 realigned" folder.
%
%   Optional name-value arguments:
%       SpatialAvgCorr - finite-aperture (spatial-averaging) correction
%                        factor (>=1) multiplied into the MI to recover the
%                        true spatial-peak pressure the 400 um tip averages
%                        out. Default 1 (off). NOTE: this INCREASES MI.
%       Fnum           - if given (and SpatialAvgCorr is not), the factor is
%                        computed from the beam geometry via spatialAvgFactor
%                        (jinc focal-plane model) at the measured f_awf.
%       ApertureUm     - hydrophone aperture for the model [um]. Default 400.
%       c              - speed of sound [m/s]. Default 1500.

    arguments
        dataDir {mustBeTextScalar} = ...
            miData('L11-5','2026-08-13_realigned')
        options.SpatialAvgCorr (1,1) double = NaN
        options.Fnum           (1,1) double = NaN
        options.ApertureUm     (1,1) double = 400
        options.c              (1,1) double = 1500
        options.fMI            (1,1) double = NaN   % [MHz] freq used in MI = pr/sqrt(f)
                                                    % and derating; default = measured f_awf
        options.Method (1,1) string {mustBeMember(options.Method, ...
                          ["singlefreq","broadband"])} = "singlefreq"
        options.flo   (1,1) double = 1     % broadband passband low edge [MHz]
        options.fhi   (1,1) double = 18    % broadband passband high edge [MHz]
        options.taper (1,1) double = 1     % broadband band-edge taper [MHz]
    end
    dataDir = char(dataDir);
    addpath(fileparts(mfilename('fullpath')));  % co-located analysis code
    repo = miRoot('calibration');

    % --- sweep + Verasonics reference ----------------------------
    V    = [2 4 6 8 10 15 20 25 30 35 40 45 50];
    Vera = [0.01 0.04 0.06 0.09 0.11 0.18 0.25 0.31 0.38 0.44 0.51 0.57 0.64];

    imp      = {'50ohm','50ohm','1e6ohm','1e6ohm'};
    pol      = {'neg','pos','neg','pos'};
    loadcorr = [1 1 0.5 0.5];      % 50 ohm ref -> x1 ; 1 MOhm high-Z -> x0.5
    lbl      = {'50 \Omega neg','50 \Omega pos','1 M\Omega neg (\div2)','1 M\Omega pos (\div2)'};

    % --- acoustic working frequency (IEC: -6 dB centre of the fundamental,
    %     from a LOW-drive capture so harmonics do not bias it) ------------
    w  = readHWS(fullfile(dataDir,'L11-5_10V_50ohm_neg.hws'));
    y  = w(1).y - median(w(1).y); dt = w(1).dt;
    Yf = abs(fft(y.*hann(numel(y)))); fax = (0:numel(Yf)-1)'/(numel(Yf)*dt);
    hb = 1:floor(numel(Yf)/2); fh = fax(hb); Yh = Yf(hb);
    bf = fh>2e6 & fh<8e6;                         % fundamental lobe only
    pk = max(Yh.*bf); idx = find(Yh>=pk/2 & bf);  % -6 dB edges
    ft = (fh(idx(1)) + fh(idx(end)))/2 / 1e6;     % MHz

    % --- measurement depth from the note -> derating -------------
    [~, cfg0] = readHWS(fullfile(dataDir,'L11-5_20V_50ohm_neg.hws'));
    tok = regexp(cfg0.notes, 'Pos\s*=\s*\(\s*[-\d.]+\s*,\s*[-\d.]+\s*,\s*([-\d.]+)', 'tokens', 'once');
    z_mm = str2double(tok{1});  Depth_cm = z_mm/10;

    % frequency used in the MI formula & derating (Verasonics uses the
    % programmed transmit frequency; default = measured acoustic working freq)
    fMI = options.fMI; if isnan(fMI), fMI = ft; end
    DerateFactor = db2mag(0.3*Depth_cm*fMI);

    % --- parallel-circuit sensitivity at ft ----------------------
    sens_total = sensitivityVperMPa(repo, ft);   % V/MPa
    tables     = load(fullfile(repo,'HydrophoneTables.mat'), ...
                      'AmplifierTable','HydrophoneTable');   % for broadband
    bandopts   = struct('flo',options.flo,'fhi',options.fhi,'taper',options.taper);

    % --- read all captures -> raw peak-negative pressure ---------
    %   prRaw [MPa] is the in-water peak-rarefactional pressure implied by
    %   each file's raw waveform (before load correction). Two methods:
    %     singlefreq - Vpk / sens_total   (original, default; reproduces prior work)
    %     broadband  - frequency-domain deconvolution (magnitude + amp phase)
    Vpk   = nan(numel(V), 4);
    prRaw = nan(numel(V), 4);
    for s = 1:4
        for i = 1:numel(V)
            f = fullfile(dataDir, sprintf('L11-5_%dV_%s_%s.hws', V(i), imp{s}, pol{s}));
            ww = readHWS(f);
            yb = ww(1).y - median(ww(1).y);
            Vpk(i,s) = -min(yb);                          % |peak negative| [V]
            if options.Method == "broadband"
                prRaw(i,s) = voltageToPressureBroadband(ww(1).t, yb, tables, bandopts);
            else
                prRaw(i,s) = Vpk(i,s) / sens_total;       % [MPa]
            end
        end
    end

    % --- finite-aperture (spatial-averaging) correction ----------
    %   The 400 um tip averages the focal field and under-reads the true
    %   spatial peak; CF (>=1) recovers it. This RAISES MI.
    if ~isnan(options.SpatialAvgCorr)
        CF = options.SpatialAvgCorr;  cfSrc = 'user-supplied';
    elseif ~isnan(options.Fnum)
        CF = spatialAvgFactor(ft, options.Fnum, options.ApertureUm, options.c);
        cfSrc = sprintf('model: F#=%.1f, a=%.0f um', options.Fnum, options.ApertureUm);
    else
        CF = 1;  cfSrc = 'none';
    end

    % --- derated MI ----------------------------------------------
    %   MI = (in-water pr)/derate/sqrt(fMI), with load correction and the
    %   optional spatial-averaging factor CF. Identical to the original
    %   expression when Method="singlefreq" (prRaw = Vpk/sens_total).
    MI    = (prRaw .* loadcorr) / DerateFactor / sqrt(fMI) * CF;
    MIunc =  prRaw(:,3:4)       / DerateFactor / sqrt(fMI) * CF;  % 1MOhm, no /2

    % --- saturation-aware linear extrapolation -------------------
    %   Fit the unsaturated region (combined polarities) and extend.
    [MIx50,  Vsat50 ] = extrapLinear(V, [MI(:,1); MI(:,2)], V);
    [MIx1e6, Vsat1e6] = extrapLinear(V, [MI(:,3); MI(:,4)], V);

    % --- console summary -----------------------------------------
    fprintf('f_awf(measured)=%.2f MHz | f_MI=%.2f MHz | depth=%.1f cm | derate=%.3f | sens_total=%.3f V/MPa\n', ...
        ft, fMI, Depth_cm, DerateFactor, sens_total);
    fprintf('Method: %s | Spatial-averaging correction: CF=%.3f (%s)\n', ...
        options.Method, CF, cfSrc);
    fprintf('Saturation onset: 50 ohm >= %g V, 1 MOhm >= %g V\n', Vsat50, Vsat1e6);
    rr = MIx50(:) ./ Vera(:);
    fprintf('Extrapolated(50 ohm)/Verasonics ratio (V>=15): mean %.2f\n', mean(rr(V>=15)));
    fprintf('\n TXV | MI50(raw) MI50(extrap) | 1e6/2(raw) | 1e6 UNcorr | Verasonics\n');
    for i = 1:numel(V)
        fprintf('%4d | %8.3f %11.3f  | %8.3f   | %8.3f   | %5.3f\n', ...
            V(i), mean(MI(i,1:2)), MIx50(i), mean(MI(i,3:4)), mean(MIunc(i,:)), Vera(i));
    end

    % --- plot ----------------------------------------------------
    fig = figure('Color','w','Position',[100 100 820 600]); hold on
    col = [0.85 0.33 0.10; 0.93 0.69 0.13; 0.00 0.45 0.74; 0.30 0.75 0.93];
    ucol = [0.60 0.0 0.0];

    % saturation shading
    xl = [0 max(V)+2];
    ysat = max([MIunc(:); MIx50(:)])*1.12;
    hSat = patch([min(Vsat50,Vsat1e6) xl(2) xl(2) min(Vsat50,Vsat1e6)], ...
                 [0 0 ysat ysat], [0.6 0.6 0.6], 'FaceAlpha',0.07, 'EdgeColor','none');

    hVera = plot(V, Vera, 'k-', 'LineWidth', 2.2, 'Marker','.', 'MarkerSize',14);

    % uncorrected 1 MOhm (old method) reference
    hu = plot(V, mean(MIunc,2), '^', 'Color', ucol, 'MarkerFaceColor', ucol, ...
              'MarkerSize',6, 'LineWidth',0.8);

    % extrapolation lines (the saturation-corrected MI)
    plot(V, MIx50,  '--', 'Color', col(1,:), 'LineWidth', 1.6, 'HandleVisibility','off');
    plot(V, MIx1e6, '--', 'Color', col(3,:), 'LineWidth', 1.6, 'HandleVisibility','off');

    % measured markers (raw, saturating)
    h = gobjects(1,4); mrk = {'o','s','o','s'};
    for s = 1:4
        if s <= 2
            h(s) = plot(V, MI(:,s), mrk{s}, 'Color', col(s,:), ...
                'MarkerFaceColor', col(s,:), 'MarkerSize',7);
        else
            h(s) = plot(V, MI(:,s), mrk{s}, 'Color', col(s,:), ...
                'MarkerSize',7, 'LineWidth',1.3);
        end
    end

    grid on; box on; xlim(xl); ylim([0 ysat]);
    xlabel('TX voltage (V)'); ylabel('Mechanical Index (MI, derated)');
    title(sprintf(['L11-5 pulse inversion (realigned, z=%.0f mm): MI vs Verasonics' ...
        ', f_{awf}=%.1f MHz, %s'], z_mm, ft, options.Method));
    legend([hVera h hu hSat], ['Verasonics reference', lbl, ...
        '1 M\Omega (UNcorrected, old method)', 'saturated region'], ...
        'Location','northwest');
    cfNote = ''; if CF ~= 1, cfNote = sprintf('\nspatial-avg correction CF=%.3f applied', CF); end
    text(0.98,0.03, sprintf(['thick dashed = linear extrapolation of unsaturated V\n' ...
        'extrapolated(50 \\Omega) \\approx %.2f\\times Verasonics%s'], mean(rr(V>=15)), cfNote), ...
        'Units','normalized','HorizontalAlignment','right','VerticalAlignment','bottom', ...
        'FontSize',9,'Color',[0.3 0.3 0.3]);

    suffix = ''; if options.Method == "broadband", suffix = '_broadband'; end
    outdir = miRoot('results'); if ~isfolder(outdir), mkdir(outdir); end
    out = fullfile(outdir, sprintf('MI_vs_TXvoltage_L11_5_realigned%s.png', suffix));
    exportgraphics(fig, out, 'Resolution', 150);
    fprintf('\nSaved figure: %s\n', out);
end

% ===================================================================
function [MIx, Vsat] = extrapLinear(V, MIstack, Vout)
%EXTRAPLINEAR Robust linear fit over the unsaturated region, evaluated at Vout.
%   MIstack may stack several polarities of the same channel (V repeated).
    Vs = [V(:); V(:)]; Vs = Vs(1:numel(MIstack)); MIstack = MIstack(:);

    % anchor on the clean low/mid region, past turn-on and pre-saturation
    anc = Vs>=4 & Vs<=20;
    p   = polyfit(Vs(anc), MIstack(anc), 1);
    pred = polyval(p, Vs);

    % saturation onset = first V>20 that drops >5% below the linear trend
    % for two consecutive voltages (robust to single noisy points)
    Vu = unique(V(V>20));
    Vsat = max(V)+1;
    for k = 1:numel(Vu)
        m  = Vs==Vu(k);
        if all(MIstack(m) < 0.95*pred(m))
            if k==numel(Vu) || all(MIstack(Vs==Vu(k+1)) < 0.95*polyval(p,Vu(k+1)))
                Vsat = Vu(k); break;
            end
        end
    end

    % refit over the unsaturated set, dropping any residual low outliers
    mask = Vs>=4 & Vs<Vsat;
    p2   = polyfit(Vs(mask), MIstack(mask), 1);
    keep = mask & (MIstack >= 0.9*polyval(p2,Vs));
    p3   = polyfit(Vs(keep), MIstack(keep), 1);
    MIx  = polyval(p3, Vout(:));
end

% ===================================================================
function sens_total = sensitivityVperMPa(repo, ft)
%SENSITIVITYVPERMPA Onda parallel-circuit loaded sensitivity [V/MPa] at ft [MHz].
    S  = load(fullfile(repo,'HydrophoneTables.mat'),'AmplifierTable','HydrophoneTable');
    HT = S.HydrophoneTable; AT = S.AmplifierTable;
    ih = find(HT.FREQ_MHz  >= ft, 1);  sens_h = db2mag(HT{ih,"SENS_DB"})*1e6; cap_h = HT{ih,"CAP_PF"}*1e-12;
    ia = find(AT.FREQ_MHZ  >= ft, 1);  gain   = db2mag(AT{ia,"GAIN_DB"});     cap_a = AT{ia,"CAP_PF"}*1e-12;
    sens_total = gain * sens_h * cap_h/(cap_a+cap_h) * 1e6;   % V/MPa
end
