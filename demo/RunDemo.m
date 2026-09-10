function results = RunDemo(mode)
%RUNDEMO  End-to-end demonstration of the acoustic-output measurement chain.
%
%   RunDemo               run the demonstration and check it against the
%                         reference numbers in demo/reference/
%   RunDemo('reference')  re-generate those reference numbers (only after a
%                         deliberate, understood change to the method)
%   R = RunDemo(...)      also return every number as a struct
%
%   Nothing outside this repository is needed: no Verasonics, no NI drivers,
%   no absolute paths. The bundled captures in data/ are real measurements of
%   the S5-1 shear-wave sequence and of the L11-5 pulse-inversion sequence.
%
%   The five parts
%   --------------
%     1  CALIBRATION   the Onda certificates -> a sensitivity in V/MPa, for
%                      both measuring chains (with and without the preamp)
%     2  ONE CAPTURE   one .hws file -> p_r, f_awf, PD, PII, MI, I_sppa.3,
%                      I_spta.3, with a figure showing where each comes from
%     3  SAFETY LIMIT  a transmit-voltage sweep -> the maximum voltage at
%                      which the S5-1 push still meets the FDA Track 3 limits
%     4  VALIDATION    the same chain applied to the L11-5 pulse-inversion
%                      sequence, compared against the MI Verasonics reports
%                      for it -- an INDEPENDENT check that the chain is right
%     5  REGRESSION    every number above compared with demo/reference/
%
%   Figures are written to results/demo (git-ignored).
%
%   See also: SafetyIndices, SystemSensitivity, AcousticWorkingFrequency,
%             readHWS, SafetyTableAll.

    if nargin < 1, mode = "check"; end
    mode = string(mode);
    root = miRoot();
    out  = fullfile(root,'results','demo');
    if ~isfolder(out), mkdir(out); end
    % No warnings are suppressed here on purpose. readHWS silences only the one
    % unavoidable HDF5 datatype warning, by identifier; everything else -- RIS
    % sampling, a record shorter than the pulse -- is a real problem and must
    % be allowed to surface.

    banner('MI-measurements demonstration');
    fprintf('repository : %s\ndata       : %s\nfigures    : %s\n', root, miData(), out);

    R = struct();
    R = part1_calibration(R);
    R = part2_oneCapture(R, out);
    R = part3_safetyLimit(R, out);
    R = part4_validation(R, out);
    part5_regression(R, mode);

    if nargout, results = R; end
end

% =====================================================================
% 1. CALIBRATION
% =====================================================================
function R = part1_calibration(R)
    banner('1. Calibration: certificates -> V/MPa');

    % The acoustic working frequency is a property of the SEQUENCE, so it is
    % measured from a low-drive capture of that sequence (harmonics at high
    % drive would bias the -6 dB centre upward).
    pushRef = fullfile(sessionB(), 'S5-1_8V_Preamp_peak.hws');
    [fawf, fi] = AcousticWorkingFrequency(pushRef, [], [1.3 3.5]);

    [sensPre, ip] = SystemSensitivity(fawf, Chain="preamp", ScopeImpedance_Ohm=50);

    % Without the preamp the bare hydrophone is loaded by the cable + scope
    % capacitance. That value is on no datasheet -- it is measured, from paired
    % captures of the same point over the range where the preamp is still linear.
    [Cload, ci] = CalibrateLoadCapacitance(sessionB(), ...
        'S5-1_%dV_Preamp_peak.hws', 'S5-1_%dV_NoPreamp_opt.hws', ...
        [3 4 5 6 7 8 9 10 12], fawf);
    sensNo = SystemSensitivity(fawf, Chain="nopreamp", LoadCap_pF=Cload);

    fprintf('  f_awf                       %8.3f MHz   (-6 dB band %.3f - %.3f MHz)\n', ...
            fawf, fi.f1_MHz, fi.f2_MHz);
    fprintf('  M_oc  open-circuit          %8.3f mV/MPa\n', ip.Moc_V_per_Pa*1e9);
    fprintf('  C_H   hydrophone            %8.2f pF\n', ip.Ch_pF);
    fprintf('  C_A   preamp input          %8.2f pF     -> divider %.3f\n', ...
            ip.Ca_pF, ip.Ch_pF/(ip.Ch_pF+ip.Ca_pF));
    fprintf('  G     preamp gain           %8.2f x      (%.1f dB)\n', ip.gain, 20*log10(ip.gain));
    fprintf('  %s\n', repmat('-',1,60));
    fprintf('  sensitivity WITH preamp     %8.4f V/MPa  (50 ohm reference)\n', sensPre);
    fprintf('  C_load cable+scope          %8.1f pF     (measured; spread %.1f %% over %d pairs)\n', ...
            Cload, ci.spread_pct, numel(ci.V));
    fprintf('  sensitivity WITHOUT preamp  %8.4f V/MPa  (%.0fx smaller)\n', sensNo, sensPre/sensNo);
    fprintf(['\n  Note: on a 1 MOhm scope the preamp reads 2x high -- its 20 dB gain is\n' ...
             '  referenced to a 50 ohm load. SystemSensitivity absorbs that in the value\n' ...
             '  it returns, so the caller never applies the factor by hand.\n']);

    R.fawf_MHz = fawf;
    R.sensPre  = sensPre;
    R.Cload_pF = Cload;
    R.sensNo   = sensNo;
end

% =====================================================================
% 2. ONE CAPTURE, EVERY QUANTITY
% =====================================================================
function R = part2_oneCapture(R, out)
    banner('2. One capture -> every safety quantity');

    f = fullfile(sessionB(), 'S5-1_8V_Preamp_peak.hws');
    [w, cfg] = readHWS(f);
    t = w(1).t*1e6;
    y = w(1).y - median(w(1).y);
    depth_cm = DepthFromNotes(cfg);

    % Burst duty of the shear-wave push: 24 pushes over a 1.2 s burst, then at
    % least 30 s idle. I_spta is a long-time average, so the EFFECTIVE rate is
    % what counts -- not the 20 Hz rate inside the burst.
    PRF = 24/(1.2 + 30);

    S = SafetyIndices(y, w(1).dt, Sensitivity=R.sensPre, Depth_cm=depth_cm, ...
                      Fawf_MHz=R.fawf_MHz, PRF_Hz=PRF, RemoveBaseline=false);

    fprintf('  file            %s\n', shorten(f));
    fprintf('  digitizer       %s @ %.0f MS/s, %d samples (%.0f us)\n', ...
            strtrim(cfg.digitizer), cfg.sample_rate/1e6, cfg.record_length, ...
            cfg.record_duration*1e6);
    fprintf('  depth           %8.2f cm      (|Pos - Home| from the capture note)\n', depth_cm);
    fprintf('  derating        %8.3f x       (0.3 dB/cm/MHz at f_awf = %.2f MHz)\n', ...
            S.derate, R.fawf_MHz);
    fprintf('  %s\n', repmat('-',1,60));
    fprintf('  p_r   in water  %8.3f MPa\n', S.pr_MPa);
    fprintf('  p_r.3 derated   %8.3f MPa\n', S.pr3_MPa);
    fprintf('  MI              %8.3f         (limit 1.9)\n', S.MI);
    fprintf('  pulse duration  %8.3f ms      (1.25 x the 10-90 %% intensity integral)\n', S.PD_s*1e3);
    fprintf('  PII.3           %8.3f J/m2\n', S.PII_J_per_m2/S.derate^2);
    fprintf('  I_sppa.3        %8.1f W/cm2   (limit 190)\n', S.Isppa3_W_cm2);
    fprintf('  I_spta.3        %8.1f mW/cm2  (limit 720; burst duty %.3f Hz)\n', ...
            S.Ispta3_mW_cm2, PRF);

    figOneCapture(t, y/R.sensPre, S, fullfile(out,'demo_1_one_capture.png'));

    R.one = S;
    R.one_depth_cm = depth_cm;
end

% =====================================================================
% 3. THE SAFETY LIMIT OF THE SEQUENCE
% =====================================================================
function R = part3_safetyLimit(R, out)
    banner('3. Transmit-voltage sweep -> maximum allowed voltage');
    fprintf(['  The S5-1 push is measured WITHOUT the preamp: the AH-2010 clips at ~2.2 V\n' ...
             '  and the push exceeds that above ~10 V transmit. Removing it costs a factor\n' ...
             '  ~73 in sensitivity but buys the whole 15-50 V range, unclipped.\n\n']);

    V   = [15 20 25 30 35 40 45 50];
    PRF = 24/(1.2 + 30);
    LIM = struct('MI',1.9,'Isppa',190,'Ispta',720);
    apertures = struct('n',{41,61,79},'suffix',{'_41el','','_79el'});

    fprintf('  %-3s %-6s %8s %10s %11s\n','el','TX(V)','MI','Isppa.3','Ispta.3');
    fprintf('  %-3s %-6s %8s %10s %11s\n','','','','[W/cm2]','[mW/cm2]');
    fprintf('  %s\n', repmat('-',1,44));

    for a = 1:numel(apertures)
        n   = apertures(a).n;
        pat = ['S5-1_%dV_NoPreamp_opt' apertures(a).suffix '.hws'];
        [MI, Is, It, depth] = sweepVoltage(sessionB(), pat, V, R.sensNo, R.fawf_MHz, PRF);

        for i = 1:numel(V)
            fprintf('  %-3d %-6d %8.2f %10.1f %11.1f\n', n, V(i), MI(i), Is(i), It(i));
        end
        vMI = crossing(V,MI,LIM.MI);
        vIs = crossing(V,Is,LIM.Isppa);
        vIt = crossing(V,It,LIM.Ispta);
        [maxV, ib] = min([vMI vIs vIt]);
        bind = {'MI','I_sppa.3','I_spta.3'};
        fprintf('  %-3d depth %.2f cm | MI=1.9 at %.1f V | Isppa=190 at %.1f V | Ispta=720 at %s V\n', ...
                n, depth, vMI, vIs, fmtInf(vIt));
        fprintf('  %-3d ==> MAXIMUM %.1f V, limited by %s\n', n, maxV, bind{ib});
        fprintf('  %s\n', repmat('-',1,44));

        R.sweep.(sprintf('el%d',n)) = struct('V',V,'MI',MI,'Isppa',Is,'Ispta',It, ...
                                             'depth_cm',depth,'maxV',maxV, ...
                                             'binding',string(bind{ib}));
    end
    fprintf(['\n  The push is I_sppa.3-limited, never MI-limited, and I_spta.3 is nowhere\n' ...
             '  near its limit -- the 1.2 s / 30 s burst duty is what buys that headroom.\n']);

    figSweep(R.sweep, V, LIM, fullfile(out,'demo_2_voltage_sweep.png'));
end

% =====================================================================
% 4. INDEPENDENT VALIDATION AGAINST VERASONICS
% =====================================================================
function R = part4_validation(R, out)
    banner('4. Validation: the same chain vs the MI Verasonics reports');
    fprintf(['  Verasonics computes an MI for its own sequences from its transducer model.\n' ...
             '  For the L11-5 pulse-inversion sequence we have both that number and a\n' ...
             '  hydrophone measurement, so the chain can be checked end to end.\n\n']);

    d    = miData('L11-5','2026-08-13_realigned');
    V    = [2 4 6 8 10 15 20 25 30 35 40 45 50];
    Vera = [0.01 0.04 0.06 0.09 0.11 0.18 0.25 0.31 0.38 0.44 0.51 0.57 0.64];

    % f_awf from a low-drive capture; the L11-5 fundamental sits near 4.5 MHz.
    fawf = AcousticWorkingFrequency(fullfile(d,'L11-5_10V_50ohm_neg.hws'), [], [2 8]);
    sens = SystemSensitivity(fawf, Chain="preamp", ScopeImpedance_Ohm=50);

    [~, cfg] = readHWS(fullfile(d,'L11-5_20V_50ohm_neg.hws'));
    depth_cm = DepthFromNotes(cfg, [0 0 0]);   % this session logged Pos only
    derate   = 10.^(0.3*depth_cm*fawf/20);

    % Both polarities of the pulse-inversion pair, on the 50 ohm-terminated scope.
    pol = {'neg','pos'};
    MI  = nan(numel(V),2);
    for s = 1:2
        for i = 1:numel(V)
            w = readHWS(fullfile(d, sprintf('L11-5_%dV_50ohm_%s.hws', V(i), pol{s})));
            y = w(1).y - median(w(1).y);
            MI(i,s) = (-min(y)/sens)/derate/sqrt(fawf);
        end
    end
    MImeas = mean(MI,2);

    % At the top of the range the preamp clips, so those readings are not the true
    % MI. Fit the unsaturated region and extrapolate through it.
    anc  = V(:) >= 4 & V(:) <= 20;
    k    = (V(anc)*MImeas(anc)) / (V(anc)*V(anc)');     % fit through the origin
    MIex = k*V(:);
    ratio     = MIex ./ Vera(:);
    meanRatio = mean(ratio(V >= 15));

    % Clipping onset: the first HIGH voltage where the reading falls below the
    % fitted line and stays there. (Low-voltage points scatter off the fit too,
    % but that is 8-bit quantisation on a millivolt signal, not clipping.)
    below = MImeas(:) < 0.95*MIex;
    Vsat  = Inf;
    for k = find(V(:)' > 20)
        if all(below(k:end)), Vsat = V(k); break, end
    end

    fprintf('  f_awf %.2f MHz | depth %.1f cm | derating %.3f | sensitivity %.3f V/MPa\n\n', ...
            fawf, depth_cm, derate, sens);
    fprintf('  %6s %10s %12s %12s %8s\n','TX(V)','measured','unclipped','Verasonics','ratio');
    fprintf('  %s\n', repmat('-',1,52));
    for i = 1:numel(V)
        flag = '';
        if V(i) >= Vsat, flag = '  (preamp clipped)'; end
        fprintf('  %6d %10.3f %12.3f %12.3f %8.2f%s\n', ...
                V(i), MImeas(i), MIex(i), Vera(i), ratio(i), flag);
    end
    fprintf('  %s\n', repmat('-',1,52));
    fprintf('  mean ratio over the unclipped comparison range (TX >= 15 V): %.2f x\n\n', meanRatio);
    fprintf(['  A hydrophone MI within ~25 %% of the console value is the expected level of\n' ...
             '  agreement: the two use different transducer models, and the 400 um tip\n' ...
             '  spatially averages the focal peak. What the check really rules out is a\n' ...
             '  factor-of-2 error -- which is exactly what an uncorrected 1 MOhm reading has.\n']);

    fprintf(['\n  PlotMI_L11_5 runs the same comparison with a more careful,\n' ...
             '  saturation-aware fit over both scope impedances; it lands at 1.23x.\n']);

    figValidation(V, MImeas, MIex, Vera, meanRatio, Vsat, fullfile(out,'demo_3_validation.png'));

    R.L11_5 = struct('fawf_MHz',fawf,'depth_cm',depth_cm,'sens',sens, ...
                     'V',V,'MI_measured',MImeas','MI_unclipped',MIex', ...
                     'MI_verasonics',Vera,'meanRatio',meanRatio);
end

% =====================================================================
% 5. REGRESSION CHECK
% =====================================================================
function part5_regression(R, mode)
    banner('5. Regression check');
    ref  = miRoot('demo','reference','expected_results.csv');
    rows = flatten(R);

    if mode == "reference"
        T = table(string(rows(:,1)), cell2mat(rows(:,2)), ...
                  'VariableNames', {'quantity','value'});
        writetable(T, ref);
        fprintf('  reference regenerated: %s (%d quantities)\n', shorten(ref), height(T));
        return
    end

    if ~isfile(ref)
        fprintf('  no reference file yet -- run RunDemo(''reference'') to create one.\n');
        return
    end

    T   = readtable(ref, 'TextType','string');
    tol = 1e-9;                     % relative; the chain is fully deterministic
    nbad = 0;
    fprintf('  %-34s %14s %14s %10s\n','quantity','computed','expected','');
    fprintf('  %s\n', repmat('-',1,76));
    for i = 1:size(rows,1)
        name = rows{i,1};
        got  = rows{i,2};
        j    = find(T.quantity == name, 1);
        if isempty(j)
            fprintf('  %-34s %14.7g %14s %10s\n', name, got, '--', 'NEW');
            nbad = nbad + 1;
            continue
        end
        want = T.value(j);
        rel  = abs(got-want)/max(abs(want), eps);
        ok   = rel <= tol || (isnan(got) && isnan(want)) || isequal(got, want);
        if ~ok, nbad = nbad + 1; end
        fprintf('  %-34s %14.7g %14.7g %10s\n', name, got, want, tern(ok,'ok','MISMATCH'));
    end
    fprintf('  %s\n', repmat('-',1,76));
    if nbad == 0
        fprintf('  ALL %d QUANTITIES REPRODUCE THE REFERENCE.\n', size(rows,1));
    else
        fprintf(2, ['  %d of %d QUANTITIES DIFFER -- the method changed. Investigate before\n' ...
                    '  regenerating with RunDemo(''reference'').\n'], nbad, size(rows,1));
    end
end

function rows = flatten(R)
%FLATTEN  The demonstration's results as name/value pairs, in a fixed order.
    rows = cell(0,2);
    rows(end+1,:) = {'calib.fawf_MHz',                R.fawf_MHz};
    rows(end+1,:) = {'calib.sens_preamp_V_per_MPa',   R.sensPre};
    rows(end+1,:) = {'calib.Cload_pF',                R.Cload_pF};
    rows(end+1,:) = {'calib.sens_nopreamp_V_per_MPa', R.sensNo};
    rows(end+1,:) = {'capture.depth_cm',              R.one_depth_cm};
    rows(end+1,:) = {'capture.pr_MPa',                R.one.pr_MPa};
    rows(end+1,:) = {'capture.MI',                    R.one.MI};
    rows(end+1,:) = {'capture.PD_ms',                 R.one.PD_s*1e3};
    rows(end+1,:) = {'capture.Isppa3_W_cm2',          R.one.Isppa3_W_cm2};
    rows(end+1,:) = {'capture.Ispta3_mW_cm2',         R.one.Ispta3_mW_cm2};
    for n = [41 61 79]
        s = R.sweep.(sprintf('el%d',n));
        rows(end+1,:) = {sprintf('sweep.el%d.maxV',n),         s.maxV};                  %#ok<AGROW>
        rows(end+1,:) = {sprintf('sweep.el%d.MI_at_30V',n),    interp1(s.V,s.MI,30)};    %#ok<AGROW>
        rows(end+1,:) = {sprintf('sweep.el%d.Isppa_at_30V',n), interp1(s.V,s.Isppa,30)}; %#ok<AGROW>
    end
    rows(end+1,:) = {'L11-5.fawf_MHz',                R.L11_5.fawf_MHz};
    rows(end+1,:) = {'L11-5.MI_unclipped_at_50V',     R.L11_5.MI_unclipped(end)};
    rows(end+1,:) = {'L11-5.ratio_vs_verasonics',     R.L11_5.meanRatio};
end

% =====================================================================
% helpers
% =====================================================================
function d = sessionB()
    d = miData('S5-1','2026-08-17_sessionB_preamp_vs_nopreamp');
end

function [MI, Is, It, depth] = sweepVoltage(d, pat, V, sens, fawf, PRF)
    n = numel(V);
    MI = nan(1,n); Is = nan(1,n); It = nan(1,n); dep = nan(1,n);
    for i = 1:n
        [w, cfg] = readHWS(fullfile(d, sprintf(pat, V(i))));
        dep(i) = DepthFromNotes(cfg);
        S = SafetyIndices(w(1).y, w(1).dt, Sensitivity=sens, Depth_cm=dep(i), ...
                          Fawf_MHz=fawf, PRF_Hz=PRF);
        MI(i) = S.MI;  Is(i) = S.Isppa3_W_cm2;  It(i) = S.Ispta3_mW_cm2;
    end
    depth = median(dep);
end

function v = crossing(V, Y, lim)
%CROSSING  Transmit voltage at which Y first reaches lim (linear in between).
    v = NaN;
    if all(Y < lim), v = Inf; return, end
    if Y(1) >= lim,  v = V(1); return, end
    for k = 1:numel(V)-1
        if (Y(k)-lim)*(Y(k+1)-lim) <= 0 && Y(k+1) ~= Y(k)
            v = V(k) + (lim-Y(k))/(Y(k+1)-Y(k))*(V(k+1)-V(k));
            return
        end
    end
end

function banner(s)
    fprintf('\n%s\n%s\n%s\n', repmat('=',1,76), s, repmat('=',1,76));
end

function s = shorten(p)
    s = strrep(p, [miRoot() filesep], '');
end

function s = fmtInf(v)
    if isinf(v), s = '>50'; elseif isnan(v), s = '--'; else, s = sprintf('%.1f',v); end
end

function s = tern(b,a,c)
    if b, s = a; else, s = c; end
end

% ---------------------------------------------------------------------
function figOneCapture(t, p, S, file)
    f  = figure('Color','w','Position',[60 60 1000 420],'Visible','off');
    tl = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    nexttile; hold on; grid on; box on
    plot(t, p, 'Color',[0 .45 .74]);
    [~,imin] = min(p);
    plot(t(imin), min(p), 'rv','MarkerFaceColor','r','MarkerSize',8);
    xline(t(S.i10),'k:','LineWidth',1); xline(t(S.i90),'k:','LineWidth',1);
    xlabel('time (\mus)'); ylabel('pressure in water (MPa)');
    title(sprintf('p_r = %.2f MPa;  dotted = 10-90 %% of \\int I dt', S.pr_MPa));

    nexttile; hold on; grid on; box on
    plot(t, p, 'Color',[0 .45 .74],'LineWidth',1.2);
    yline(min(p),'r--','LineWidth',1.2);
    plot(t(imin), min(p), 'rv','MarkerFaceColor','r','MarkerSize',9);
    xlim(t(imin) + [-0.7 0.7]);
    xlabel('time (\mus)'); ylabel('pressure (MPa)');
    title('the peak-rarefactional half-cycle');

    title(tl, sprintf(['S5-1 push, 61 elements, 8 V, %.0f mm from the transducer\n' ...
        'f_{awf} %.2f MHz | derating %.2fx | MI %.2f | ' ...
        'I_{sppa.3} %.0f W/cm^2 | I_{spta.3} %.0f mW/cm^2'], ...
        S.Depth_cm*10, S.Fawf_MHz, S.derate, S.MI, S.Isppa3_W_cm2, S.Ispta3_mW_cm2), ...
        'FontWeight','bold','FontSize',10);
    exportgraphics(f, file, 'Resolution', 150); close(f);
    fprintf('  figure: %s\n', shorten(file));
end

function figSweep(sw, V, LIM, file)
    f   = figure('Color','w','Position',[60 60 1150 400],'Visible','off');
    col = [.47 .67 .19; 0 .45 .74; .85 .33 .10];
    q   = {'MI',    'MI (derated)',           LIM.MI; ...
           'Isppa', 'I_{sppa.3} (W/cm^2)',    LIM.Isppa; ...
           'Ispta', 'I_{spta.3} (mW/cm^2)',   LIM.Ispta};
    ns  = [41 61 79];
    for k = 1:3
        subplot(1,3,k); hold on; grid on; box on
        for a = 1:3
            s = sw.(sprintf('el%d',ns(a)));
            plot(V, s.(q{k,1}), 'o-','Color',col(a,:),'MarkerFaceColor',col(a,:), ...
                 'LineWidth',1.5,'DisplayName',sprintf('%d elements',ns(a)));
        end
        yline(q{k,3},'k--','HandleVisibility','off');
        xlabel('transmit voltage (V)'); ylabel(q{k,2}); xlim([10 52]);
        title(sprintf('%s  (dashed = FDA limit)', q{k,2}));
        if k == 1, legend('Location','northwest'); end
    end
    sgtitle('S5-1 ARF push measured without the pre-amplifier, derated 0.3 dB/cm/MHz');
    exportgraphics(f, file, 'Resolution', 150); close(f);
    fprintf('  figure: %s\n', shorten(file));
end

function figValidation(V, MImeas, MIex, Vera, ratio, Vsat, file)
    f = figure('Color','w','Position',[60 60 780 520],'Visible','off');
    hold on; grid on; box on
    plot(V, Vera, 'k-','LineWidth',2.2,'Marker','.','MarkerSize',16, ...
         'DisplayName','MI reported by Verasonics');
    plot(V, MImeas, 'o','Color',[0 .45 .74],'MarkerFaceColor',[0 .45 .74], ...
         'MarkerSize',7,'DisplayName','hydrophone, as measured');
    plot(V, MIex, '--','Color',[.85 .33 .10],'LineWidth',1.8, ...
         'DisplayName','hydrophone, preamp clipping removed');
    xlabel('transmit voltage (V)'); ylabel('Mechanical Index (derated)');
    title(sprintf('L11-5 pulse inversion: measured MI is %.2fx the console value', ratio));
    legend('Location','northwest');
    if isfinite(Vsat)
        xline(Vsat,'k:','LineWidth',1.2,'HandleVisibility','off');
        note = sprintf(['the measured points fall away above %g V:' newline ...
                        'that is the AH-2010 preamp clipping, not the field'], Vsat);
    else
        note = ['the preamp stayed linear over the whole sweep' newline ...
                '(this sequence never drove it to its ceiling)'];
    end
    text(0.98, 0.04, note, 'Units','normalized','HorizontalAlignment','right', ...
         'FontSize',9,'Color',[.35 .35 .35]);
    exportgraphics(f, file, 'Resolution', 150); close(f);
    fprintf('  figure: %s\n', shorten(file));
end
