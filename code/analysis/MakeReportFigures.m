function MakeReportFigures(dataRoot, outDir)
%MAKEREPORTFIGURES  Generate the figures of the MI-measurement how-to report.
%
%   MakeReportFigures() uses the default data/output locations.
%   MakeReportFigures(dataRoot, outDir) overrides them.
%
%   Produces, in <outDir>:
%     calibration.png      hydrophone / preamp calibration vs frequency and the
%                          combined system sensitivity with and without preamp
%     waveform_metrics.png one raw .hws capture and every quantity read from it
%                          (p_r, f_awf from the -6 dB band, pulse duration, PII)
%     preamp_clipping.png  preamp clipping vs transmit voltage, and the
%                          no-preamp chain that replaces it above the ceiling
%     axial_deration.png   axial sweep: why the DERATED maximum sits proximal
%                          of the in-water maximum
%     reflection_check.png does the long burst show tank reflections? (envelope
%                          step test + its sensitivity, calibrated by injection)
%
%   Everything it needs is in the repository: readHWS on the path, the calibration
%   in calibration/, and the raw .hws captures in data/ (resolved via miData, so
%   nothing here has to be edited after a clone).
%
%   Requires MATLAB (developed on R2025b). No toolboxes beyond base + Signal.

    if nargin < 1 || isempty(dataRoot)
        dataRoot = miData('S5-1','2026-08-17_sessionB_preamp_vs_nopreamp');
    end
    if nargin < 2 || isempty(outDir)
        outDir = miRoot('results');
    end
    here = fileparts(mfilename('fullpath'));
    addpath(here);
    if ~isfolder(outDir), mkdir(outDir); end

    C = loadCalibration(miRoot('calibration','HydrophoneTables.mat'));

    % --- the S5-1 push chain: acoustic working frequency and sensitivities ---
    fawf = measureFawf(fullfile(dataRoot, 'S5-1_8V_Preamp_peak.hws'), [1.3 3.5]);
    [Moc, Ch, gain, Ca] = calAt(C, fawf);
    sensPre = gain * Moc * Ch / (Ch + Ca) * 1e6;         % V/MPa, with preamp
    Cload   = estimateLoadCapacitance(dataRoot, C, fawf); % pF, cable+scope
    sensNo  = Moc * Ch / (Ch + Cload*1e-12) * 1e6;        % V/MPa, no preamp

    fprintf('f_awf = %.2f MHz | sens(preamp) = %.3f V/MPa | C_load = %.1f pF | sens(no preamp) = %.2f mV/MPa\n', ...
        fawf, sensPre, Cload, sensNo*1e3);

    figCalibration(C, Cload, [fawf 5.0], fullfile(outDir,'calibration.png'));
    figWaveformMetrics(dataRoot, sensPre, fawf, fullfile(outDir,'waveform_metrics.png'));
    figPreampClipping(dataRoot, sensPre, sensNo, fullfile(outDir,'preamp_clipping.png'));
    figAxialDeration(dataRoot, sensPre, fawf, fullfile(outDir,'axial_deration.png'));
    figReflectionCheck(dataRoot, fawf, fullfile(outDir,'reflection_check.png'));
end

% =====================================================================
%  Figure 5 -- is the long burst contaminated by tank reflections?
% =====================================================================
function figReflectionCheck(d, fawf, outFile)
%FIGREFLECTIONCHECK  The envelope test of the reflection section, run over the
%   axial sweep, plus a calibration of its sensitivity by injecting a synthetic
%   reflection of known relative amplitude.
%
%   Principle: a reflection arriving Dt after the direct field adds coherently
%   from that moment on, so it STEPS the burst envelope by up to +-r, where r
%   is its relative amplitude. Its phase depends on position, so across the
%   axial stations the steps would vary in sign; a flat envelope at every
%   station bounds r. Independently, reverberation would keep arriving for Dt
%   AFTER the transmit stops, giving a decaying tail.

    ks = 1:10; n = numel(ks);
    z = nan(n,1); rmsv = nan(n,1); tail = nan(n,1); droop = nan(n,1);
    E = cell(n,1); T = cell(n,1);
    for i = 1:n
        f = fullfile(d,sprintf('S5-1_4V_Preamp_%d.hws',ks(i)));
        [w,cfg] = readHWS(f);
        t = w(1).t*1e6; y = w(1).y - median(w(1).y); dt = w(1).dt;
        z(i) = depthFromNotes(cfg);
        [T{i},E{i},rmsv(i),tail(i)] = burstEnvelope(y,t,dt,fawf);
        % head-to-tail decay: smooth, monotonic and position-independent, so it
        % is the transmit supply sagging over the burst -- NOT a reflection,
        % which would step at a fixed delay with a position-dependent sign.
        droop(i) = mean(E{i}(1:3)) - mean(E{i}(end-2:end));
    end

    % --- sensitivity calibration: inject a known reflection at 200 us -------
    [w,~] = readHWS(fullfile(d,'S5-1_4V_Preamp_peak.hws'));
    y0 = w(1).y - median(w(1).y); dt = w(1).dt; t0 = w(1).t*1e6;
    lag = round(200e-6/dt);
    rInj = [0 0.01 0.02 0.03 0.05 0.10 0.20];
    rmsInj = nan(size(rInj));
    for j = 1:numel(rInj)
        yr = y0 + rInj(j)*[zeros(lag,1); y0(1:end-lag)];
        [~,~,rmsInj(j),~] = burstEnvelope(yr,t0,dt,fawf);
    end
    % invert the calibration at the worst observed station
    rBound = interp1(rmsInj(2:end), rInj(2:end), max(rmsv), 'linear', 'extrap');

    fig = figure('Color','w','Position',[60 60 1020 430]);
    tl = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    % ---- (a) the measured envelopes -------------------------------------
    nexttile; hold on; grid on; box on
    cmap = parula(n);
    for i = 1:n
        plot(T{i}-T{i}(1), E{i}, 'Color', [cmap(i,:) 0.85], 'LineWidth', 0.9);
    end
    yline(1,'k--','LineWidth',1);
    ylim([0.95 1.05]); xlim([0 850]);
    xlabel('time since start of the burst (\mus)');
    ylabel('burst envelope, normalised to its median');
    title(sprintf(['(a) %d axial stations, 4 V: no step at any station\n' ...
        'the common %.1f %% droop is smooth and position-independent (supply sag)'], ...
        n, 100*mean(droop)));
    cb = colorbar; cb.Label.String = 'distance from transducer (mm)';
    clim([min(z) max(z)]);

    % ---- (b) how big a reflection would have been visible ---------------
    nexttile; hold on; grid on; box on
    plot(100*rInj, 100*rmsInj, 'o-','Color',[.85 .33 .10], ...
        'MarkerFaceColor',[.85 .33 .10],'LineWidth',1.5, ...
        'DisplayName','synthetic reflection injected');
    yl = [0 1.05*100*max(rmsInj)];
    patch([0 22 22 0],[100*min(rmsv) 100*min(rmsv) 100*max(rmsv) 100*max(rmsv)], ...
        [0 .45 .74],'FaceAlpha',.16,'EdgeColor','none', ...
        'DisplayName','observed, all stations');
    plot([0 100*rBound],[100*max(rmsv) 100*max(rmsv)],'k:','HandleVisibility','off');
    plot([100*rBound 100*rBound],[0 100*max(rmsv)],'k:','HandleVisibility','off');
    plot(100*rBound, 100*max(rmsv),'kp','MarkerSize',14,'MarkerFaceColor','y', ...
        'DisplayName',sprintf('upper bound: r < %.0f %%',100*rBound));
    xlim([0 22]); ylim(yl);
    xlabel('reflection amplitude r, relative to direct field (%)');
    ylabel('envelope rms (%)');
    title(sprintf(['(b) sensitivity of the test\n' ...
        'post-burst tail decays in %.1f-%.1f \\mus (%.0f-%.0f cycles)'], ...
        min(tail), max(tail), min(tail)*fawf, max(tail)*fawf));
    legend('Location','northwest','FontSize',8);

    title(tl, ['Reflection check on the 1900-cycle push: a reflection would step the burst ' ...
        'envelope' newline 'when it arrives, and would ring on after the transmit stops. ' ...
        'Neither is seen.'], 'FontWeight','bold','FontSize',10);
    exportgraphics(fig, outFile, 'Resolution', 200); close(fig);
    fprintf(['saved: %s\n  envelope rms %.1f-%.1f %%, tail %.1f-%.1f us ' ...
             '=> reflection amplitude r < %.0f %% of the direct field\n'], ...
        outFile, 100*min(rmsv), 100*max(rmsv), min(tail), max(tail), 100*rBound);
end

function [te,e,rmsv,tail] = burstEnvelope(y,t,dt,fawf)
%BURSTENVELOPE  Block-RMS envelope over the steady part of the burst, its rms
%   deviation, and the time for the post-burst tail to fall below 5 %.
%
%   Block RMS rather than a peak envelope: with an 8-bit digitizer a peak
%   envelope carries +-1 code of quantisation dither (~2 % here), which buries
%   the effect being looked for. Averaging |y|^2 over a 20 us block (2000
%   samples) suppresses that by ~sqrt(N) while leaving any step from a
%   reflection --- a persistent shift, not a fluctuation --- fully intact.
    ncyc = max(2,round(2/(fawf*1e6)/dt));
    pk   = movmax(abs(y), ncyc);            % only to locate the burst edges
    thr  = 0.2*max(pk);
    on   = find(pk>thr,1); off = find(pk>thr,1,'last');
    i1 = on + round(0.02*(off-on)); i2 = off - round(0.02*(off-on));

    ns = round(20e-6/dt); seg = y(i1:i2); ts = t(i1:i2);
    nb = floor(numel(seg)/ns);
    e  = sqrt(mean(reshape(seg(1:nb*ns), ns, nb).^2, 1))';
    te = mean(reshape(ts(1:nb*ns), ns, nb), 1)';
    e  = e/median(e);
    rmsv = std(e);

    lvl = median(pk(i1:i2)); post = pk(off:end); tp = t(off:end);
    b = find(post < 0.05*lvl, 1);
    tail = NaN; if ~isempty(b), tail = tp(b)-tp(1); end
end

% =====================================================================
%  Figure 1 -- calibration as a function of frequency
% =====================================================================
function figCalibration(C, Cload, fMark, outFile)
    f = C.HT.FREQ_MHz;                       % 1 - 20 MHz, the hydrophone grid
    Moc  = db2mag(C.HT.SENS_DB) * 1e6;       % V/Pa (open circuit)
    Ch   = C.HT.CAP_PF;                      % pF
    gain = db2mag(interp1(C.AT.FREQ_MHZ, C.AT.GAIN_DB,   f));
    gdB  =         interp1(C.AT.FREQ_MHZ, C.AT.GAIN_DB,   f);
    phs  =         interp1(C.AT.FREQ_MHZ, C.AT.PHASE_DEG, f);
    Ca   =         interp1(C.AT.FREQ_MHZ, C.AT.CAP_PF,    f);

    sensPre = gain .* Moc .* Ch ./ (Ch + Ca) * 1e6;      % V/MPa
    sensNo  =         Moc .* Ch ./ (Ch + Cload) * 1e6;   % V/MPa

    fig = figure('Color','w','Position',[60 60 1000 720]);
    tl  = tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

    nexttile; hold on; grid on; box on
    plot(f, C.HT.SENS_DB, 'LineWidth',1.5,'Color',[0 .45 .74]);
    markFreqs(fMark);
    xlabel('frequency (MHz)'); ylabel('M_c (dB re 1 V/\muPa)');
    title('(a) hydrophone open-circuit sensitivity');
    yyaxis right; plot(f, Moc*1e9,'--','LineWidth',1.2); ylabel('M_c (nV/Pa)');
    ax = gca; ax.YAxis(2).Color = [.4 .4 .4];

    nexttile; hold on; grid on; box on
    plot(f, gdB, 'LineWidth',1.5,'Color',[.85 .33 .10]); markFreqs(fMark);
    xlabel('frequency (MHz)'); ylabel('gain (dB)');
    title('(b) AH-2010-025 preamp gain and phase');
    yyaxis right; plot(f, phs, '--','LineWidth',1.2); ylabel('phase (deg)');
    ax = gca; ax.YAxis(2).Color = [.4 .4 .4];

    nexttile; hold on; grid on; box on
    plot(f, Ch, 'LineWidth',1.5,'Color',[0 .45 .74],'DisplayName','C_H  hydrophone');
    plot(f, Ca, 'LineWidth',1.5,'Color',[.85 .33 .10],'DisplayName','C_A  preamp input');
    yline(Cload,'k--','LineWidth',1.2,'DisplayName',sprintf('C_{load} cable+scope = %.0f pF',Cload));
    markFreqs(fMark);
    set(gca,'YScale','log'); ylim([4 300]);
    xlabel('frequency (MHz)'); ylabel('capacitance (pF)');
    title('(c) capacitances of the divider'); legend('Location','east','FontSize',8);

    nexttile; hold on; grid on; box on
    h1 = plot(f, sensPre, 'LineWidth',1.8,'Color',[.85 .33 .10]);
    markFreqs(fMark);
    xlabel('frequency (MHz)'); ylabel('with preamp (V/MPa)'); ylim([0 2]);
    yyaxis right
    h2 = plot(f, sensNo*1e3,'LineWidth',1.8,'Color',[0 .45 .74]);
    ylabel('no preamp (mV/MPa)'); ylim([0 25]);
    ax = gca; ax.YAxis(1).Color = [.85 .33 .10]; ax.YAxis(2).Color = [0 .45 .74];
    title('(d) combined system sensitivity M_L(f)');
    legend([h1 h2], {'with preamp (50 \Omega ref.)','no preamp (into C_{load})'}, ...
        'Location','southeast','FontSize',8);

    title(tl, sprintf(['Onda HGL-0400 SN 1037 + AH-2010-025 SN 1109, calibration 10-Jun-2024\n' ...
        'dashed vertical lines: f_{awf} of the S5-1 push (%.2f MHz) and the L11-5 pulse (%.2f MHz)'], ...
        fMark(1), fMark(2)), 'FontWeight','bold');
    exportgraphics(fig, outFile, 'Resolution', 200); close(fig);
    fprintf('saved: %s\n', outFile);
end

function markFreqs(fMark)
    for k = 1:numel(fMark)
        xline(fMark(k), 'k:', 'LineWidth', 1, 'HandleVisibility','off');
    end
end

% =====================================================================
%  Figure 2 -- every quantity we read out of one saved waveform
% =====================================================================
function figWaveformMetrics(d, sens, fawf, outFile)
    rho = 1000; c = 1500;
    [w, cfg] = readHWS(fullfile(d,'S5-1_8V_Preamp_peak.hws'));
    t = w(1).t*1e6; y = w(1).y - median(w(1).y); dt = w(1).dt;
    p = y / sens;                                       % MPa
    depth_mm = depthFromNotes(cfg);                     % mm from transducer face

    % pulse duration and PII, exactly as in CalculateSafety.m
    I  = (p*1e6).^2/(rho*c);                            % W/m^2
    cR = rescale(cumsum(I),0,1);
    i10 = find(cR>0.1,1); i90 = find(cR>0.9,1);
    PD  = 1.25*(t(i90)-t(i10));                         % us
    PII = trapz(I(i10:i90))*dt;                         % J/m^2

    % spectrum and the -6 dB acoustic working frequency
    Y = abs(fft(y.*hann(numel(y)))); fax = (0:numel(Y)-1)'/(numel(Y)*dt)/1e6;
    hb = 1:floor(numel(Y)/2); fh = fax(hb); Yh = Yf2db(Y(hb));
    band = fh>1.3 & fh<3.5; pk = max(Yh(band)); idx = find(Yh>=pk-6 & band);
    f1 = fh(idx(1)); f2 = fh(idx(end));

    [~, imin] = min(p); tzoom = t(imin) + [-2 2];

    fig = figure('Color','w','Position',[60 60 1000 720]);
    tl = tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

    nexttile; hold on; grid on; box on
    plot(t, p, 'Color',[0 .45 .74]);
    xlabel('time (\mus)'); ylabel('pressure (MPa)');
    title(sprintf('(a) full record: %.0f \\mus at %.0f MS/s', t(end)-t(1), 1/dt/1e6));
    xlim([t(1) t(end)]);

    nexttile; hold on; grid on; box on
    plot(t, p, 'Color',[0 .45 .74],'LineWidth',1.2);
    yline(min(p),'r--','LineWidth',1.2);
    plot(t(imin), min(p), 'rv','MarkerFaceColor','r','MarkerSize',8);
    text(t(imin), min(p), sprintf('  p_r = %.2f MPa', -min(p)), ...
        'Color','r','VerticalAlignment','bottom','FontWeight','bold');
    xlim(tzoom); xlabel('time (\mus)'); ylabel('pressure (MPa)');
    title('(b) peak rarefactional pressure p_r (in water)');

    axSpec = nexttile; hold on; grid on; box on
    plot(fh, Yh-pk, 'Color',[.85 .33 .10],'LineWidth',1.2);
    yline(-6,'k--','-6 dB','LineWidth',1);
    xlim([0 8]); ylim([-60 3]); xlabel('frequency (MHz)'); ylabel('magnitude (dB re peak)');
    title('(c) acoustic working frequency: -6 dB centre of the fundamental');
    text(4.5, -13, '2^{nd}','HorizontalAlignment','center','FontSize',8,'Color',[.3 .3 .3]);
    text(6.7, -27, '3^{rd}','HorizontalAlignment','center','FontSize',8,'Color',[.3 .3 .3]);
    nexttile; hold on; grid on; box on
    plot(t, cR, 'Color',[.47 .67 .19],'LineWidth',1.4);
    yline(0.1,'k--'); yline(0.9,'k--');
    xline(t(i10),'k:','LineWidth',1); xline(t(i90),'k:','LineWidth',1);
    xlabel('time (\mus)'); ylabel('normalised cumulative \intI dt');
    title(sprintf('(d) pulse duration PD = 1.25(t_{90}-t_{10}) = %.0f \\mus', PD));
    xlim([t(1) t(end)]);

    der  = db2mag(0.3*depth_mm/10*fawf);
    MI   = (-min(p))/der/sqrt(fawf);
    Ider = (1/der)^2;
    Isppa = PII*Ider/(PD*1e-6)/1e4;
    title(tl, sprintf(['S5-1 push, 61 elements, 8 V, %.0f mm from the transducer\n' ...
        'p_r = %.2f MPa   f_{awf} = %.2f MHz   derating %.2f dB (x%.2f)   ' ...
        'MI = %.2f   I_{sppa.3} = %.0f W/cm^2'], ...
        depth_mm, -min(p), fawf, 0.3*depth_mm/10*fawf, der, MI, Isppa), 'FontWeight','bold');

    % Inset on (c): the fundamental lobe of an 1900-cycle burst is far narrower
    % than the 0-8 MHz scale can show. Added last -- a bare axes() call before
    % the layout title would invalidate the tiledlayout handle.
    drawnow; pos = axSpec.Position;
    axIn = axes(fig,'Position',[pos(1)+0.14*pos(3) pos(2)+0.44*pos(4) 0.40*pos(3) 0.42*pos(4)]);
    hold(axIn,'on'); grid(axIn,'on'); box(axIn,'on');
    plot(axIn, fh, Yh-pk, 'Color',[.85 .33 .10],'LineWidth',1.2);
    yline(axIn,-6,'k--','LineWidth',1);
    xline(axIn,f1,'k:','LineWidth',1); xline(axIn,f2,'k:','LineWidth',1);
    xline(axIn,(f1+f2)/2,'b-','LineWidth',1.5);
    xlim(axIn,[(f1+f2)/2-0.12 (f1+f2)/2+0.12]); ylim(axIn,[-20 2]);
    axIn.FontSize = 7;
    title(axIn, sprintf('zoom: f_{awf} = %.2f MHz',(f1+f2)/2),'Color','b','FontSize',8);

    exportgraphics(fig, outFile, 'Resolution', 200); close(fig);
    fprintf('saved: %s\n', outFile);
end

function db = Yf2db(Y)
    db = 20*log10(max(Y, eps));
end

% =====================================================================
%  Figure 3 -- preamp clipping and the no-preamp alternative
% =====================================================================
function figPreampClipping(d, sensPre, sensNo, outFile)
    Vp = [2 3 4 5 6 7 8 9 10 12 15 20];
    Vn = [2 3 4 5 6 7 8 9 10 12 15 20 25 30 35 40 45 50];
    vnegP = arrayfun(@(v) vneg(fullfile(d,sprintf('S5-1_%dV_Preamp_peak.hws',v))),  Vp);
    vnegN = arrayfun(@(v) vneg(fullfile(d,sprintf('S5-1_%dV_NoPreamp_opt.hws',v))), Vn);

    % linear fit through the region where the preamp is demonstrably linear
    lin = Vp <= 8;
    kP  = Vp(lin)' \ vnegP(lin)';
    linN = Vn <= 20;
    kN  = Vn(linN)' \ vnegN(linN)';

    fig = figure('Color','w','Position',[60 60 1020 420]);
    tl = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    % (a) waveform zooms: the flat top appears
    nexttile; hold on; grid on; box on
    show = [4 8 15 20]; cmap = lines(numel(show));
    for k = 1:numel(show)
        w = readHWS(fullfile(d,sprintf('S5-1_%dV_Preamp_peak.hws',show(k))));
        t = w(1).t*1e6; y = w(1).y - median(w(1).y);
        [~,i0] = min(y);
        plot(t - t(i0), y, 'Color', cmap(k,:), 'LineWidth',1.2, ...
            'DisplayName', sprintf('%d V', show(k)));
    end
    xlim([-1.5 1.5]); xlabel('time relative to the negative peak (\mus)');
    ylabel('preamp output (V)');
    title('(a) preamp output clips at \approx -2.2 V');
    legend('Location','southeast','FontSize',8);

    % (b) peak-negative voltage vs transmit voltage, both chains
    nexttile; hold on; grid on; box on
    plot(Vp, vnegP, 'o-','Color',[.85 .33 .10],'MarkerFaceColor',[.85 .33 .10], ...
        'LineWidth',1.4,'DisplayName','with preamp (left axis)');
    plot([0 22], kP*[0 22], 'k--','LineWidth',1,'DisplayName','linear fit, \leq 8 V');
    ylabel('preamp output |V_{neg}| (V)'); ylim([0 3]);
    xlabel('transmit voltage (V)');
    yyaxis right
    plot(Vn, vnegN, 's-','Color',[0 .45 .74],'MarkerFaceColor',[0 .45 .74], ...
        'LineWidth',1.4,'DisplayName','no preamp (right axis)');
    plot([0 52], kN*[0 52], ':','Color',[0 .45 .74],'LineWidth',1.2,'HandleVisibility','off');
    ylabel('bare hydrophone |V_{neg}| (V)');
    ax = gca; ax.YAxis(2).Color = [0 .45 .74];
    xlim([0 52]);
    title('(b) the two chains, same acoustic point');
    legend('Location','northwest','FontSize',8);

    title(tl, sprintf(['S5-1 push, 61 elements: the preamp is the ceiling, not the scope\n' ...
        'sensitivity with preamp %.2f V/MPa, without preamp %.1f mV/MPa (%.0fx smaller)'], ...
        sensPre, sensNo*1e3, sensPre/sensNo), 'FontWeight','bold');
    exportgraphics(fig, outFile, 'Resolution', 200); close(fig);
    fprintf('saved: %s\n', outFile);
end

% =====================================================================
%  Figure 4 -- why the derated maximum is proximal of the in-water maximum
% =====================================================================
function figAxialDeration(d, sens, fawf, outFile)
%FIGAXIALDERATION  Two axial sweeps: a shallow push (40 mm focus, small shift)
%   and the deep focused imaging beam (79 mm focus, large shift). The point is
%   the same in both: the maximum you must report is the maximum of the DERATED
%   profile, and it is never distal of the in-water maximum.
    base = miData('S5-1','2026-08-13_sessionA_preamp_50ohm');

    fig = figure('Color','w','Position',[60 60 1020 430]);
    tl = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    axialPanel(d, 'S5-1_4V_Preamp_%d.hws', 1:10, sens, fawf, ...
        'push, 61 el, 4 V, focus 40 mm');
    axialPanel(fullfile(base,'focused_imaging_axial'), 'S5-1_10V_50ohm_pos_foc_%d.hws', 1:15, ...
        sens, NaN, 'focused imaging, 10 V, focus 79 mm');

    title(tl, ['Axial sweeps: at every depth the lateral/elevation maximum was located first' newline ...
        'Report the maximum of the DERATED profile (blue): it sits at the' newline ...
        'proximal edge of the in-water plateau'], 'FontWeight','bold','FontSize',10);
    exportgraphics(fig, outFile, 'Resolution', 200); close(fig);
    fprintf('saved: %s\n', outFile);
end

function axialPanel(d, pat, ks, sens, fawf, ttl)
    if isnan(fawf)
        fawf = measureFawf(fullfile(d,sprintf(pat,ks(1))), [1.3 3.5]);
    end
    n = numel(ks); dep = nan(n,1); pr = nan(n,1);
    for i = 1:n
        [w,cfg] = readHWS(fullfile(d,sprintf(pat,ks(i))));
        y = w(1).y - median(w(1).y);
        dep(i) = depthFromNotes(cfg);
        pr(i)  = -min(y)/sens;                          % MPa, in water
    end
    [dep,o] = sort(dep); pr = pr(o);
    der = db2mag(0.3*dep/10*fawf);
    MIw = pr/sqrt(fawf);  MI3 = pr./der/sqrt(fawf);
    [~,i3] = max(MI3);

    % The in-water profile is a broad PLATEAU, not a sharp peak (and the 8-bit
    % digitizer quantises V_neg, so several depths read identically). Quoting a
    % single "in-water maximum" would be meaningless; the plateau is the honest
    % object to compare the derated maximum against.
    onPlateau = MIw >= 0.99*max(MIw);
    zPlat = dep(onPlateau);
    zMid  = (zPlat(1)+zPlat(end))/2;

    nexttile; hold on; grid on; box on
    yl = [0 1.14];
    patch([zPlat(1) zPlat(end) zPlat(end) zPlat(1)], [yl(1) yl(1) yl(2) yl(2)], ...
        [.5 .5 .5],'FaceAlpha',.12,'EdgeColor','none','DisplayName','in-water plateau (top 1 %)');
    plot(dep, MIw/max(MIw), 'o-','Color',[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5], ...
        'LineWidth',1.4,'DisplayName','in water (not derated)');
    plot(dep, MI3/max(MI3), 's-','Color',[0 .45 .74],'MarkerFaceColor',[0 .45 .74], ...
        'LineWidth',1.6,'DisplayName','derated, 0.3 dB/cm/MHz');
    plot(dep(i3), 1,'p','MarkerSize',16,'MarkerFaceColor',[0 .45 .74], ...
        'MarkerEdgeColor','k','DisplayName','derated maximum');
    ylabel('MI, normalised to own maximum'); ylim(yl);
    yyaxis right
    plot(dep, 1./der, 'k--','LineWidth',1.2,'DisplayName','derating factor');
    ylabel('derating factor 10^{-0.3 z f/20}');
    ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = [.3 .3 .3];
    yyaxis left
    xlabel('distance from the transducer (mm)');
    title(sprintf(['%s\nf_{awf} = %.2f MHz\nderated max %.0f mm | plateau %.0f-%.0f mm ' ...
        '(centre %.0f mm)'], ttl, fawf, dep(i3), zPlat(1), zPlat(end), zMid), 'FontSize', 9);
    legend('Location','south','FontSize',7.5);
end

% =====================================================================
%  helpers
% =====================================================================
function C = loadCalibration(matFile)
    S = load(matFile, 'HydrophoneTable', 'AmplifierTable');
    C.HT = S.HydrophoneTable; C.AT = S.AmplifierTable;
end

function [Moc, Ch, gain, Ca] = calAt(C, ft)
%CALAT  Calibration values at frequency ft [MHz]. Moc [V/Pa], Ch/Ca [F], gain [-].
    Moc  = db2mag(interp1(C.HT.FREQ_MHz, C.HT.SENS_DB, ft))*1e6;
    Ch   =        interp1(C.HT.FREQ_MHz, C.HT.CAP_PF,  ft)*1e-12;
    gain = db2mag(interp1(C.AT.FREQ_MHZ, C.AT.GAIN_DB, ft));
    Ca   =        interp1(C.AT.FREQ_MHZ, C.AT.CAP_PF,  ft)*1e-12;
end

function Cload = estimateLoadCapacitance(d, C, fawf)
%ESTIMATELOADCAPACITANCE  Cable+scope capacitance [pF] seen by the bare
%   hydrophone, from the preamp / no-preamp voltage ratio at the same point
%   over the range where the preamp is unsaturated.
    [~, Ch, gain, Ca] = calAt(C, fawf);
    V = [3 4 5 6 7 8 9 10 12];
    r = arrayfun(@(v) vneg(fullfile(d,sprintf('S5-1_%dV_Preamp_peak.hws',v))) / ...
                      vneg(fullfile(d,sprintf('S5-1_%dV_NoPreamp_opt.hws',v))), V);
    Cload = (mean(r)*(Ch+Ca)/gain - Ch)*1e12;
end

function v = vneg(file)
    w = readHWS(file); v = -min(w(1).y - median(w(1).y));
end

function fawf = measureFawf(file, bandMHz)
%MEASUREFAWF  IEC acoustic working frequency: centre of the -6 dB band of the
%   fundamental, from a LOW-drive capture (harmonics would bias it upward).
    w = readHWS(file); y = w(1).y - median(w(1).y); dt = w(1).dt;
    Y = abs(fft(y.*hann(numel(y)))); f = (0:numel(Y)-1)'/(numel(Y)*dt);
    hb = 1:floor(numel(Y)/2);
    bf = f(hb) > bandMHz(1)*1e6 & f(hb) < bandMHz(2)*1e6;
    pk = max(Y(hb).*bf); idx = find(Y(hb) >= pk/2 & bf);
    fawf = (f(idx(1)) + f(idx(end)))/2/1e6;
end

function dep = depthFromNotes(cfg)
%DEPTHFROMNOTES  Distance [mm] from the transducer centre (`Home`) to `Pos`,
%   both written into the .hws note at acquisition time.
    h = regexp(cfg.notes,'Home\s*=\s*\(\s*([-\d.]+)\s*,\s*([-\d.]+)\s*,\s*([-\d.]+)','tokens','once');
    p = regexp(cfg.notes,'Pos\s*=\s*\(\s*([-\d.]+)\s*,\s*([-\d.]+)\s*,\s*([-\d.]+)','tokens','once');
    dep = norm(str2double(p) - str2double(h));
end
