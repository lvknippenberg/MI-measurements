function S = DemoSingleCapture(hwsFile, opts)
%DEMOSINGLECAPTURE  One hydrophone capture in; the safety numbers and one figure out.
%
%   DemoSingleCapture                 run on a bundled example capture
%   DemoSingleCapture(file)           run on any .hws capture
%   S = DemoSingleCapture(...)        also return every quantity as a struct
%
%   This is the whole measurement chain at its smallest: read a capture, work
%   out how it was acquired from the note stored inside it, convert volts to
%   pressure, and report MI, I_sppa.3 and I_spta.3 with a figure that shows
%   where each number came from. For the fuller demonstration -- calibration,
%   a voltage sweep, validation against Verasonics, a regression check -- see
%   RunDemo.
%
%   Name-value options
%   ------------------
%     PRF_Hz      pulse repetition frequency, needed for I_spta.3 only. This is
%                 a property of the SEQUENCE, not of the capture, so it cannot
%                 be read from the file. For a burst protocol give the
%                 EFFECTIVE rate (pushes per burst / (burst + idle)), not the
%                 rate inside the burst. Default: NaN -> I_spta.3 is NaN.
%     Band_MHz    where to look for the fundamental when measuring f_awf.
%                 Default [1 8]; widen or narrow it for a different probe.
%     LoadCap_pF  cable+scope capacitance, used only when the capture was taken
%                 WITHOUT the pre-amplifier. Default 123.3 pF, the value
%                 measured for this setup by CalibrateLoadCapacitance -- redo
%                 that measurement if the cabling has changed.
%     Home_mm     transducer-centre coordinate, if the note records only Pos.
%     Depth_cm, ScopeImpedance_Ohm, WithPreamp
%                 supply what the capture note does not record. Without them
%                 an incomplete note is prompted for interactively, and under
%                 `matlab -batch` it is an error rather than a guess.
%     OutDir      where to write the figure. Default results/demo.
%
%   The figure opens on screen and is saved as a .png at the same time; its
%   handle comes back in S.figure.
%
%   Example
%   -------
%     % the S5-1 focused imaging beam, 10 V, scanned at 71 lines x 270 us
%     S = DemoSingleCapture(miData('S5-1','2026-08-13_sessionA_preamp_50ohm', ...
%             'focused_imaging_sweep','S5-1_10V_50ohm_pos_foc.hws'), ...
%             PRF_Hz = 1/(71*270e-6));
%
%   See also: RunDemo, SafetyIndices, SystemSensitivity, readHWS.

    arguments
        hwsFile          {mustBeTextScalar} = ""
        opts.PRF_Hz      (1,1) double = NaN
        opts.Band_MHz    (1,2) double = [1 8]
        opts.LoadCap_pF  (1,1) double = 123.3
        opts.Home_mm           double = []
        opts.Depth_cm    (1,1) double = NaN
        opts.ScopeImpedance_Ohm double = []
        opts.WithPreamp        = []
        opts.OutDir      {mustBeTextScalar} = ""
    end

    % ---- the default example: a short focused imaging pulse, easy to read ----
    usingDefault = (strlength(hwsFile) == 0);
    if usingDefault
        hwsFile = miData('S5-1','2026-08-13_sessionA_preamp_50ohm', ...
                         'focused_imaging_sweep','S5-1_10V_50ohm_pos_foc.hws');
        if isnan(opts.PRF_Hz)
            opts.PRF_Hz = 1/(71*270e-6);      % 71 lines per frame, 270 us apart
        end
    end
    hwsFile = char(hwsFile);
    outDir  = char(opts.OutDir);
    if isempty(outDir), outDir = miRoot('results','demo'); end
    if ~isfolder(outDir), mkdir(outDir); end

    % ================================================================== read
    [w, cfg] = readHWS(hwsFile);
    t  = w(1).t;
    y  = w(1).y - median(w(1).y);             % remove the DC offset
    dt = w(1).dt;

    % ---- how was it acquired? ---------------------------------------------
    % The note the operator wrote into the file says. Anything it does not say
    % is asked for rather than assumed -- see CaptureConditions for which
    % fields are required and what happens when one is absent.
    ccArgs = {};
    if ~isnan(opts.Depth_cm),           ccArgs = [ccArgs {'Depth_cm', opts.Depth_cm}]; end
    if ~isempty(opts.Home_mm),          ccArgs = [ccArgs {'Home_mm', opts.Home_mm}]; end
    if ~isempty(opts.ScopeImpedance_Ohm)
        ccArgs = [ccArgs {'ScopeImpedance_Ohm', opts.ScopeImpedance_Ohm}];
    end
    if ~isempty(opts.WithPreamp),       ccArgs = [ccArgs {'WithPreamp', opts.WithPreamp}]; end

    C = CaptureConditions(cfg, ccArgs{:});
    chain    = C.chain;
    Rscope   = C.scopeImpedance_Ohm;
    depth_cm = C.depth_cm;

    % ================================================================ convert
    fawf = AcousticWorkingFrequency(y, dt, opts.Band_MHz);

    if chain == "preamp"
        sens = SystemSensitivity(fawf, Chain="preamp", ScopeImpedance_Ohm=Rscope);
        chainTxt = sprintf('with pre-amplifier, scope on %s', ohmStr(Rscope));
    else
        sens = SystemSensitivity(fawf, Chain="nopreamp", LoadCap_pF=opts.LoadCap_pF);
        chainTxt = sprintf('no pre-amplifier, C_load = %.1f pF', opts.LoadCap_pF);
    end

    S = SafetyIndices(y, dt, Sensitivity=sens, Depth_cm=depth_cm, ...
                      Fawf_MHz=fawf, PRF_Hz=opts.PRF_Hz, RemoveBaseline=false);
    S.file  = hwsFile;
    S.chain = chain;

    % ================================================================= report
    fprintf('\n%s\n', repmat('=',1,68));
    fprintf('  %s\n', shorten(hwsFile));
    fprintf('%s\n', repmat('=',1,68));
    fprintf('  ACQUISITION   (read from the capture)\n');
    % Not every capture stores the scope configuration block; the waveform's own
    % time axis always carries the sampling, so derive it from there rather than
    % printing blanks.
    dig = strtrim(cfg.digitizer);
    if isempty(dig), dig = '(not recorded)'; end
    fprintf('    digitizer        %s, %.0f MS/s, %d samples (%.1f us)\n', ...
            dig, 1/dt/1e6, numel(y), numel(y)*dt*1e6);
    fprintf('    chain            %s%s\n', chainTxt, srcTag(C.source.chain));
    fprintf('    depth            %.2f cm from the transducer face%s\n', ...
            depth_cm, srcTag(C.source.depth_cm));
    if C.probe ~= ""
        fprintf('    probe            %s', C.probe);
        if ~isnan(C.pushElements)
            fprintf(', %d elements, %d cycles', C.pushElements, C.pushCycles);
        end
        fprintf('\n');
    end
    fprintf('\n  CONVERSION\n');
    fprintf('    f_awf            %.3f MHz   (-6 dB centre of the fundamental)\n', fawf);
    fprintf('    sensitivity      %.4f V/MPa\n', sens);
    fprintf('    derating         %.3f x     (0.3 dB/cm/MHz over %.2f cm)\n', S.derate, depth_cm);
    fprintf('\n  RESULT\n');
    fprintf('    p_r  in water    %8.3f MPa\n', S.pr_MPa);
    fprintf('    p_r.3 derated    %8.3f MPa\n', S.pr3_MPa);
    fprintf('    pulse duration   %8.3f us\n', S.PD_s*1e6);
    fprintf('    PII.3            %8.4f J/m2\n', S.PII_J_per_m2/S.derate^2);
    fprintf('    %s\n', limitLine('MI', S.MI, 1.9, '%8.3f', ''));
    fprintf('    %s\n', limitLine('I_sppa.3', S.Isppa3_W_cm2, 190, '%8.1f', 'W/cm2'));
    if isnan(S.Ispta3_mW_cm2)
        fprintf('    I_spta.3           -- pass PRF_Hz to get it (it is a property\n');
        fprintf('                          of the sequence, not of this capture)\n');
    else
        fprintf('    %s   at PRF %.3f Hz\n', ...
                limitLine('I_spta.3', S.Ispta3_mW_cm2, 720, '%8.1f', 'mW/cm2'), opts.PRF_Hz);
    end
    fprintf('%s\n', repmat('-',1,68));
    fprintf('  Limits shown are the FDA peripheral-vessel row; pick the row for the\n');
    fprintf('  intended application (see docs/report/main.pdf, section 1.1).\n');

    % ================================================================= figure
    [~, stem] = fileparts(hwsFile);
    out = fullfile(outDir, ['capture_' matlab.lang.makeValidName(stem) '.png']);
    S.figure = makeFigure(t, y, sens, S, fawf, opts.Band_MHz, chainTxt, stem, out);
    fprintf('  figure: %s\n\n', shorten(out));

    if nargout == 0, clear S; end
end

% =====================================================================
function fig = makeFigure(t, y, sens, S, fawf, band, chainTxt, stem, out)
%MAKEFIGURE  The waveform with every quantity that was read off it.
    tus = t*1e6;
    p   = y/sens;                                   % MPa
    [~, imin] = min(p);
    [~, imax] = max(p);

    blue = [0 .45 .74];  red = [.85 .10 .10];  green = [.30 .60 .20];
    % Visible: this is a demonstration -- the figure is meant to be looked at.
    % It is still written to disk as well. Under `matlab -batch` there is no
    % display, and MATLAB handles that by itself.
    fig = figure('Color','w','Position',[60 60 1150 760],'Name', ...
                 ['MI measurement -- ' stem], 'NumberTitle','off');
    tl  = tiledlayout(fig,2,2,'Padding','compact','TileSpacing','compact');

    % ---------------- (a) the whole record ----------------------------
    ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    yl = [min(p) max(p)]*1.35;
    % the 10-90 % window that defines the pulse duration
    patch(ax, tus([S.i10 S.i90 S.i90 S.i10]), [yl(1) yl(1) yl(2) yl(2)], ...
          green,'FaceAlpha',.10,'EdgeColor','none','DisplayName','10-90 % of \int I dt');
    plot(ax, tus, p, 'Color', blue, 'LineWidth', 1.0, 'DisplayName','pressure');
    plot(ax, tus(imin), p(imin), 'v', 'Color', red, 'MarkerFaceColor', red, ...
         'MarkerSize', 8, 'DisplayName', sprintf('p_r = %.2f MPa', -p(imin)));
    plot(ax, tus(imax), p(imax), '^', 'Color', [.9 .6 0], 'MarkerFaceColor', [.9 .6 0], ...
         'MarkerSize', 8, 'DisplayName', sprintf('p_+ = %.2f MPa', p(imax)));
    yline(ax, p(imin), '--', 'Color', red, 'HandleVisibility','off');
    ylim(ax, yl); xlim(ax, [tus(1) tus(end)]);
    xlabel(ax, 'time (\mus)'); ylabel(ax, 'pressure in water (MPa)');
    title(ax, '(a) the capture, and what is read off it');
    legend(ax, 'Location','southeast','FontSize',8);

    % ---------------- (b) zoom on the rarefactional peak --------------
    ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    halfWin = 1.5/fawf;                              % 1.5 cycles either side, in us
    plot(ax, tus, p, 'Color', blue, 'LineWidth', 1.6);
    plot(ax, tus(imin), p(imin), 'v', 'Color', red, 'MarkerFaceColor', red, 'MarkerSize', 9);
    yline(ax, p(imin), '--', 'Color', red, 'LineWidth', 1.1);
    yline(ax, 0, ':', 'Color', [.5 .5 .5]);
    text(ax, tus(imin), p(imin), sprintf(' p_r = %.3f MPa ', -p(imin)), ...
         'Color', red, 'FontWeight','bold', 'VerticalAlignment','bottom', ...
         'HorizontalAlignment','left', 'BackgroundColor','w', 'Margin', 1);
    xlim(ax, tus(imin) + [-halfWin halfWin]);
    xlabel(ax, 'time (\mus)'); ylabel(ax, 'pressure (MPa)');
    title(ax, sprintf('(b) the rarefactional peak -- MI uses this, not V_{pp}/2'));

    % ---------------- (c) spectrum and f_awf --------------------------
    ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    n   = numel(y);
    win = 0.5 - 0.5*cos(2*pi*(0:n-1)'/(n-1));
    Y   = abs(fft(y.*win));
    f   = (0:n-1)'/(n*dtOf(t))/1e6;
    h   = 1:floor(n/2);
    f = f(h); Y = Y(h);
    in  = f > band(1) & f < band(2);
    Yd  = 20*log10(max(Y,eps));  Yd = Yd - max(Yd(in));
    edges = f(Yd >= -6 & in);
    plot(ax, f, Yd, 'Color', [.85 .33 .10], 'LineWidth', 1.2);
    yline(ax, -6, 'k--', 'LineWidth', 1);
    xline(ax, edges(1),   'k:', 'LineWidth', 1);
    xline(ax, edges(end), 'k:', 'LineWidth', 1);
    xline(ax, fawf, '-', 'Color', blue, 'LineWidth', 1.8);
    text(ax, fawf, -3, sprintf('  f_{awf} = %.2f MHz', fawf), 'Color', blue, 'FontWeight','bold');
    xlim(ax, [0 min(band(2)*2, f(end))]); ylim(ax, [-50 4]);
    xlabel(ax, 'frequency (MHz)'); ylabel(ax, 'magnitude (dB re peak)');
    title(ax, sprintf('(c) f_{awf}: centre of the -6 dB band (%.2f - %.2f MHz)', ...
                      edges(1), edges(end)));

    % ---------------- (d) the intensity integral ----------------------
    ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    I  = (p*1e6).^2/(1000*1500);
    cR = rescale(cumsum(I),0,1);
    plot(ax, tus, cR, 'Color', green, 'LineWidth', 1.6);
    yline(ax, 0.1, 'k--'); yline(ax, 0.9, 'k--');
    xline(ax, tus(S.i10), 'k:', 'LineWidth', 1);
    xline(ax, tus(S.i90), 'k:', 'LineWidth', 1);
    plot(ax, tus([S.i10 S.i90]), cR([S.i10 S.i90]), 'ko', 'MarkerFaceColor','k','MarkerSize',5);
    % bottom-right, so the annotation never sits on top of the curve
    text(ax, 0.97, 0.06, sprintf('PD = 1.25(t_{90}-t_{10}) = %s', usStr(S.PD_s*1e6)), ...
         'Units','normalized','HorizontalAlignment','right', ...
         'VerticalAlignment','bottom','FontSize',9,'BackgroundColor','w','Margin',2);
    ylim(ax, [0 1.05]); xlim(ax, [tus(1) tus(end)]);
    xlabel(ax, 'time (\mus)'); ylabel(ax, 'normalised cumulative \int I dt');
    title(ax, '(d) pulse duration, and PII from the same window');

    title(tl, sprintf(['%s  --  %s,  %.2f cm deep\n' ...
        'MI = %.2f   |   I_{sppa.3} = %.1f W/cm^2   |   I_{spta.3} = %s   ' ...
        '(derated %.2f\\times)'], ...
        strrep(stem,'_','\_'), chainTxt, S.Depth_cm, S.MI, S.Isppa3_W_cm2, ...
        nanStr(S.Ispta3_mW_cm2), S.derate), 'FontWeight','bold','FontSize',11);

    exportgraphics(fig, out, 'Resolution', 150);
    drawnow;                       % make sure it is actually on screen
end

% =====================================================================
function s = limitLine(name, value, lim, fmt, unit)
    s = sprintf(['%-16s ' fmt ' %-7s'], name, value, unit);
    if value > lim
        s = sprintf('%s  EXCEEDS the %g limit', s, lim);
    else
        s = sprintf('%s  (limit %g, %.0f %% of it)', s, lim, 100*value/lim);
    end
end

function s = srcTag(src)
%SRCTAG  Say so when a value did not come from the capture note.
    switch string(src)
        case "note",     s = '';
        case "argument", s = '   [supplied as an argument]';
        case "prompt",   s = '   [ASKED FOR -- not in the note]';
        case "default",  s = '   [assumed -- not in the note]';
        otherwise,       s = '';
    end
end

function s = ohmStr(R)
    if R >= 1e5, s = '1 MOhm'; else, s = sprintf('%g Ohm', R); end
end

function s = nanStr(v)
    if isnan(v), s = 'not evaluated'; else, s = sprintf('%.1f mW/cm^2', v); end
end

function s = usStr(v)
%USSTR  Pulse durations here span 0.3 us to 850 us; pick a sane precision.
    if v < 10, s = sprintf('%.3f \\mus', v);
    else,      s = sprintf('%.1f \\mus', v);
    end
end

function d = dtOf(t)
    d = t(2) - t(1);
end

function s = shorten(p)
    s = strrep(p, [miRoot() filesep], '');
end
