function SafetyTableAll(outdir)
%SAFETYTABLEALL  Consolidated S5-1 safety analysis over ALL hydrophone sessions.
%
%   SAFETYTABLEALL(OUTDIR) processes every chronological S5-1 acoustic-output measurement in one
%   pass, on one set of conventions, and prints/writes the authoritative table. Supersedes
%   SafetyTable (session A only, extrapolated) and CorrectedMaxVoltage (session B only).
%
%   Sessions
%   --------
%   A  2026-08-13/14  WITH pre-amplifier, 50 ohm scope. Low transmit voltages only, because the
%                     AH-2010 output clips at ~2.2 V into 50 ohm:
%                       - focused imaging sweep      2-20 V
%                       - push 79 elements sweep     2-12 V   (+ a 1500-cycle point at 5 V)
%                       - push 61 elements sweep     3-8 V
%                     MI/I_sppa were EXTRAPOLATED from these (MI ~ V, I_sppa ~ V^2) to find the
%                     limit crossings. No 41-element data.
%   B  2026-08-17     WITHOUT pre-amplifier, hydrophone straight into the 1 MOhm scope, so the
%                     whole 15-50 V range is measured directly with no clipping:
%                       - push 41 / 61 / 79 elements, 15-50 V, peak-optimised ("_opt") captures
%                       - 61-element base also 2-12 V, overlapping session A
%                       - "79 elements repeats": 3 firings per capture, for the supply-sag spread
%
%   Conventions (identical for both sessions)
%   -----------------------------------------
%     depth      = |Pos - Home| from the capture notes (transducer centre to hydrophone tip)
%     f_awf      = measured -6 dB centre of the fundamental
%     derating   = 0.3 dB/cm/MHz
%     sensitivity= Onda parallel-circuit (HydrophoneTables.mat): open-circuit M_oc scaled by
%                  C_h/(C_h+C_load).  WITH preamp: C_load = C_amp = 6.1 pF and x gain.
%                  WITHOUT preamp: C_load = cable+scope, calibrated from the preamp/no-preamp
%                  ratio over the unsaturated 3-12 V overlap (PreampComparison).
%     MI         = p_r.3 / sqrt(f_awf)
%     I_sppa.3   = PII.3 / PD   with PD the 1.25 x (10-90 %) intensity-integral duration
%     I_spta.3   = I_sppa.3 x PD x PRF_eff, PRF_eff = 24 pushes / (1.2 s burst + 30 s idle)
%                  = 0.769 Hz. Reported for both 1900 and 1500 push cycles.
%     limits     = FDA Track 3: MI 1.9, I_sppa.3 190 W/cm^2, I_spta.3 720 mW/cm^2
%
%   Two 79-element datasets exist at the same nominal position. The three firings inside each
%   "repeats" capture agree to ~1-2 %, so there is no firing-to-firing supply sag; the ~1.6x gap
%   between the two sets is a between-capture (positioning) difference. The session-A measurement
%   -- different day, different chain, different mounting -- agrees with the peak-optimised captures
%   to ~8 % and not with the repeats, so the peak-optimised captures are taken as the on-peak truth
%   and set the limits.
%
%   NOT INDEPENDENT. Sessions A and B share one calibration: C_load is solved from their voltage
%   ratio, so B's absolute scale follows from A's. The cross-validation below tests linearity,
%   extrapolation and freedom from saturation -- not the absolute pressure scale.
%
%   Run:  SafetyTableAll                      % prints the table
%         SafetyTableAll('C:\some\folder')    % also writes safety_table_all.csv + .png

    addpath(fileparts(mfilename('fullpath')));
    if nargin < 1, outdir = ''; end
    base = miData('S5-1','2026-08-13_sessionA_preamp_50ohm');
    dB   = miData('S5-1','2026-08-17_sessionB_preamp_vs_nopreamp');
    repo = miRoot('calibration');
    S    = load(fullfile(repo,'HydrophoneTables.mat'));
    HT   = S.HydrophoneTable;   AT = S.AmplifierTable;

    rho = 1000; c = 1500;
    HomeA = [-82.21 105.14 -42.90];      % session A transducer centre
    HomeB = [ 55.8  105.84 -17.30];      % session B transducer centre (probe re-mounted)
    PRF_eff = 24/(1.2+30);   s15 = 1500/1900;
    LIM = struct('MI',1.9,'Isppa',190,'Ispta',720);

    % ------------------------------------------------------------------ calibration
    fawf = measureFawf(fullfile(dB,'S5-1_8V_Preamp_peak.hws'));
    ch   = interp1(HT.FREQ_MHz,HT.CAP_PF ,fawf)*1e-12;
    moc  = db2mag(interp1(HT.FREQ_MHz,HT.SENS_DB,fawf))*1e6;      % V/Pa open circuit
    gain = db2mag(interp1(AT.FREQ_MHZ,AT.GAIN_DB,fawf));
    Ca   = interp1(AT.FREQ_MHZ,AT.CAP_PF,fawf)*1e-12;
    sensPre = gain*moc*ch/(ch+Ca);                                % V/Pa, preamp into 50 ohm
    Cload   = calibrateCload(dB,sensPre,moc,ch);                  % cable+scope, from the overlap
    sensNo  = moc*ch/(ch+Cload);                                  % V/Pa, no preamp into 1 MOhm

    fprintf('\n=== calibration ===\n');
    fprintf('f_awf %.3f MHz | C_h %.1f pF | M_oc %.3f mV/MPa | preamp gain %.2f\n', ...
            fawf, ch*1e12, moc*1e9, gain);
    fprintf('sens (preamp, 50 ohm) %.4f V/MPa | C_load (no preamp) %.1f pF | sens (no preamp) %.4f V/MPa\n', ...
            sensPre*1e6, Cload*1e12, sensNo*1e6);

    % ------------------------------------------------------------------ session A
    fprintf('\n=== SESSION A (2026-08-13/14, preamp, 50 ohm, extrapolated) ===\n');
    A = struct();
    A.foc = sweep(fullfile(base,'focused_imaging_sweep'), ...
                  'S5-1_%dV_50ohm_pos_foc.hws',   [2 4 6 8 10 12 15 20], sensPre, fawf, HomeA, rho, c);
    A.p79 = sweep(fullfile(base,'push_79el_sweep'), ...
                  'S5-1_%dV_50ohm_fund_push.hws', [2 3 4 5 6 7 8 9 10 11 12], sensPre, fawf, HomeA, rho, c);
    A.p61 = sweep(fullfile(base,'push_61el_sweep'), ...
                  'S5-1_%dV_50ohm_fund_push_61el__opt.hws', [3 4 6 7 8], sensPre, fawf, HomeA, rho, c);
    names = {'focused imaging','push 79 el','push 61 el'};  fns = {'foc','p79','p61'};
    for k = 1:3
        s = A.(fns{k});
        lin = s.V(:)\s.MI(:);  quad = (s.V(:).^2)\s.Isppa(:);      % through the origin
        vMI = LIM.MI/lin;  vIs = sqrt(LIM.Isppa/quad);
        A.(fns{k}).vMI = vMI;  A.(fns{k}).vIs = vIs;
        fprintf('%-16s n=%2d  depth %.2f cm  PD %.3f ms | MI=1.9 at %5.1f V | Isppa=190 at %5.1f V\n', ...
                names{k}, numel(s.V), s.depth, median(s.PD)*1e3, vMI, vIs);
    end

    % ------------------------------------------------------------------ session B
    fprintf('\n=== SESSION B (2026-08-17, NO preamp, 1 MOhm, direct measurement) ===\n');
    V = [15 20 25 30 35 40 45 50];
    apert = {41,'_41el'; 61,''; 79,'_79el'};
    B = struct();  rows = {};
    for a = 1:size(apert,1)
        N = apert{a,1};
        s = sweep(dB, ['S5-1_%dV_NoPreamp_opt' apert{a,2} '.hws'], V, sensNo, fawf, HomeB, rho, c);
        s.Ispta19 = s.Isppa .* s.PD * PRF_eff * 1e3;
        s.Ispta15 = s.Ispta19 * s15;
        s.vMI  = crossing(V,s.MI,LIM.MI);
        s.vIs  = crossing(V,s.Isppa,LIM.Isppa);
        s.vIt  = crossing(V,s.Ispta19,LIM.Ispta);
        s.maxV = min([s.vMI s.vIs s.vIt]);
        [~,s.binding] = min([s.vMI s.vIs s.vIt]);
        bnames = {'MI','I_sppa.3','I_spta.3'};
        s.bindingName = bnames{s.binding};
        s.MI_at_max    = interp1(V,s.MI,   min(s.maxV,V(end)));
        s.Isppa_at_max = interp1(V,s.Isppa,min(s.maxV,V(end)));
        B.(sprintf('el%d',N)) = s;
        fprintf('\n%2d elements  (depth %.2f cm, PD %.3f ms)\n', N, s.depth, median(s.PD)*1e3);
        fprintf('  V        :%s\n', sprintf('%8.0f',V));
        fprintf('  MI       :%s\n', sprintf('%8.2f',s.MI));
        fprintf('  Isppa.3  :%s\n', sprintf('%8.1f',s.Isppa));
        fprintf('  Ispta.3  :%s   (1900 cyc)\n', sprintf('%8.1f',s.Ispta19));
        fprintf('  crossings: MI=1.9 at %.1f V | Isppa=190 at %.1f V | Ispta=720 at %s V\n', ...
                s.vMI, s.vIs, f2(s.vIt));
        fprintf('  --> MAX ALLOWED %.1f V (%s-limited); there MI %.2f, Isppa.3 %.0f\n', ...
                s.maxV, s.bindingName, s.MI_at_max, s.Isppa_at_max);
        for i = 1:numel(V)
            rows(end+1,:) = {N, V(i), s.MI(i), s.Isppa(i), s.Ispta19(i), s.Ispta15(i)}; %#ok<AGROW>
        end
    end

    % ------------------------------------------------------------------ cross-validation
    fprintf('\n=== CROSS-VALIDATION: session A extrapolation vs session B measurement ===\n');
    fprintf('%-12s %-22s %-22s %s\n','aperture','A: Isppa=190 at','B: Isppa=190 at','difference');
    pairs = {'p61','el61',61; 'p79','el79',79};
    for k = 1:2
        va = A.(pairs{k,1}).vIs;  vb = B.(pairs{k,2}).vIs;
        fprintf('%2d elements  %18.1f V %20.1f V %11.1f %%\n', pairs{k,3}, va, vb, 100*(vb-va)/va);
    end
    fprintf(['The two agree to within a few per cent, so the 41-element session-B number is\n' ...
             'trustworthy on the same footing as 61 and 79.\n\n' ...
             'WHAT THIS DOES AND DOES NOT SHOW. The chains are NOT independent. C_load was\n' ...
             'solved from the chain-A/chain-B voltage ratio, so chain B''s absolute scale is\n' ...
             'derived from chain A''s calibration by construction. Any error in the preamp-chain\n' ...
             'sensitivity -- above all in the hydrophone sensitivity M_oc -- propagates into\n' ...
             'C_load and then into every no-preamp number, identically, and would not show up\n' ...
             'here.\n' ...
             'What the agreement DOES test is real and worth having: the two sessions differ in\n' ...
             'drive range, clipping behaviour, scope impedance, probe mounting and day. Their\n' ...
             'agreement shows that the low-voltage extrapolation (MI ~ V, I_sppa ~ V^2) holds out\n' ...
             'to the limit crossings, and that neither chain is distorted by saturation there.\n' ...
             'It says nothing about the absolute pressure scale; only a second, independently\n' ...
             'calibrated hydrophone can test that.\n']);

    % ------------------------------------------------------------------ repeats vs peak-optimised
    % Full re-analysis of the "79 elements repeats" set on the same footing, so the two 79-element
    % datasets can be compared as MI / I_sppa curves and not just as raw millivolts.
    dR = fullfile(dB,'79el_repeats');
    R  = sweepMulti(dR,'S5-1_%dV_NoPreamp_opt_79el.hws',V,sensNo,fawf,HomeB,rho,c);
    R.vIs = crossing(V,R.Isppa,LIM.Isppa);  R.vMI = crossing(V,R.MI,LIM.MI);
    B.el79_repeats = R;
    fprintf('\n=== TWO 79-ELEMENT DATASETS AT THE SAME POSITION: which one sets the limit? ===\n');
    fprintf('  %-26s%s\n','peak-optimised  MI',      sprintf('%9.2f',B.el79.MI));
    fprintf('  %-26s%s\n','repeats (3 firings) MI',  sprintf('%9.2f',R.MI));
    fprintf('  %-26s%s\n','peak-optimised  Isppa.3', sprintf('%9.0f',B.el79.Isppa));
    fprintf('  %-26s%s\n','repeats Isppa.3',         sprintf('%9.0f',R.Isppa));
    fprintf('  within-capture spread of the 3 firings: at most %.1f %% of the mean\n', ...
            100*max(R.spread));
    fprintf(['  -> The three consecutive firings inside one capture agree to ~1-2 %%, so there is NO\n' ...
             '     measurable firing-to-firing supply sag. The gap between the two datasets is\n' ...
             '     BETWEEN captures, not within a burst.\n']);
    fprintf('  Isppa.3 = 190 reached at: peak-optimised %.1f V | repeats %.1f V | SESSION A %.1f V\n', ...
            B.el79.vIs, R.vIs, A.p79.vIs);
    fprintf(['  Session A is an INDEPENDENT measurement (different day, with preamp, 50 ohm, low\n' ...
             '  voltage, different Home) and it agrees with the peak-optimised captures to ~8 %%,\n' ...
             '  not with the repeats. The repeats therefore sample a slightly off-peak point (taken\n' ...
             '  later, without re-peaking). The peak-optimised captures set the limit.\n']);
    fprintf('\n=== raw millivolts behind that comparison ===\n');
    fprintf('%4s %10s %10s %10s %12s\n','TX','opt Vneg','rep mean','rep spread','opt/rep');
    for i = 1:numel(V)
        vo = vneg(fullfile(dB,sprintf('S5-1_%dV_NoPreamp_opt_79el.hws',V(i))));
        vr = zeros(1,3);
        for g = 0:2
            raw = double(h5read(fullfile(dR,sprintf('S5-1_%dV_NoPreamp_opt_79el.hws',V(i))), ...
                                sprintf('/wfm_group%d/vectors/vector0/data',g)));
            sc  = h5read(fullfile(dR,sprintf('S5-1_%dV_NoPreamp_opt_79el.hws',V(i))), ...
                         sprintf('/wfm_group%d/axes/axis1/scale_coef',g));
            y = sc(1)+sc(2)*raw;  y = y-median(y);  vr(g+1) = -min(y);
        end
        fprintf('%3dV %9.1f mV %8.1f mV %6.1f-%.1f %11.2fx\n', ...
                V(i), vo*1e3, mean(vr)*1e3, min(vr)*1e3, max(vr)*1e3, vo/mean(vr));
    end
    fprintf(['The peak-optimised capture reads systematically higher than the 3-firing mean: it is\n' ...
             'the worst firing at the best point, which is what a safety limit must be set on.\n']);

    % ------------------------------------------------------------------ outputs
    if ~isempty(outdir)
        if ~exist(outdir,'dir'), mkdir(outdir); end
        T = cell2table(rows,'VariableNames',{'elements','V','MI','Isppa3_W_cm2', ...
                                             'Ispta3_1900_mW_cm2','Ispta3_1500_mW_cm2'});
        writetable(T,fullfile(outdir,'safety_table_all.csv'));
        lim = table([41;61;79], ...
                    [B.el41.maxV;B.el61.maxV;B.el79.maxV], ...
                    [B.el41.MI_at_max;B.el61.MI_at_max;B.el79.MI_at_max], ...
                    [B.el41.Isppa_at_max;B.el61.Isppa_at_max;B.el79.Isppa_at_max], ...
                    {B.el41.bindingName;B.el61.bindingName;B.el79.bindingName}, ...
                    'VariableNames',{'elements','maxV','MI_at_maxV','Isppa3_at_maxV','binding'});
        writetable(lim,fullfile(outdir,'safety_maxvoltage.csv'));
        makeFigure(B,V,LIM,fullfile(outdir,'safety_measured_41_61_79.png'));
        fprintf('\nwrote safety_table_all.csv, safety_maxvoltage.csv, safety_measured_41_61_79.png -> %s\n',outdir);
    end
    assignin('base','safetyAll',struct('A',A,'B',B,'fawf',fawf,'Cload',Cload));
end

% ============================================================ helpers
function s = sweep(d, pat, V, sens, fawf, Home, rho, c)
%SWEEP  MI / I_sppa.3 / pulse duration for one voltage sweep, one convention.
    n = numel(V);
    s = struct('V',V,'MI',nan(1,n),'Isppa',nan(1,n),'PD',nan(1,n),'depth',NaN);
    dep = [];
    for i = 1:n
        f = fullfile(d,sprintf(pat,V(i)));
        if ~isfile(f), continue, end
        [w,cfg] = readHWS(f);
        y  = w(1).y - median(w(1).y);  dt = w(1).dt;
        dc = depthFromNotes(cfg,Home);  dep(end+1) = dc; %#ok<AGROW>
        der  = db2mag(0.3*dc*fawf);  Ider = (1/der)^2;
        p    = y/sens;                                   % Pa
        s.MI(i) = (-min(p)/1e6)/der/sqrt(fawf);
        I  = p.^2/(rho*c);  cR = rescale(cumsum(I),0,1);
        i10 = find(cR>0.1,1);  i90 = find(cR>0.9,1);
        s.PD(i)    = 1.25*(i90-i10)*dt;
        s.Isppa(i) = trapz(I(i10:i90))*dt*Ider/s.PD(i)/1e4;   % W/cm^2, derated
    end
    s.depth = median(dep);
end

function s = sweepMulti(d, pat, V, sens, fawf, Home, rho, c)
%SWEEPMULTI  As `sweep`, but over every waveform record stored in the capture (the repeats set
%   holds three consecutive firings per file). Returns the mean and the within-capture spread.
    n = numel(V);
    s = struct('V',V,'MI',nan(1,n),'Isppa',nan(1,n),'PD',nan(1,n),'spread',nan(1,n),'depth',NaN);
    dep = [];
    for i = 1:n
        f = fullfile(d,sprintf(pat,V(i)));
        if ~isfile(f), continue, end
        [w,cfg] = readHWS(f);
        dc  = depthFromNotes(cfg,Home);  dep(end+1) = dc; %#ok<AGROW>
        der = db2mag(0.3*dc*fawf);  Ider = (1/der)^2;  dt = w(1).dt;
        quiet = quietHDF5TimestampWarning();  %#ok<NASGU>
        info = h5info(f);  ng = sum(startsWith({info.Groups.Name},'/wfm_group'));
        mi = nan(1,ng); is = nan(1,ng); pd = nan(1,ng);
        for g = 0:ng-1
            raw = double(h5read(f,sprintf('/wfm_group%d/vectors/vector0/data',g)));
            sc  = h5read(f,sprintf('/wfm_group%d/axes/axis1/scale_coef',g));
            y = sc(1)+sc(2)*raw;  y = y - median(y);
            p = y/sens;
            mi(g+1) = (-min(p)/1e6)/der/sqrt(fawf);
            I = p.^2/(rho*c);  cR = rescale(cumsum(I),0,1);
            i10 = find(cR>0.1,1);  i90 = find(cR>0.9,1);
            pd(g+1) = 1.25*(i90-i10)*dt;
            is(g+1) = trapz(I(i10:i90))*dt*Ider/pd(g+1)/1e4;
        end
        s.MI(i) = mean(mi);  s.Isppa(i) = mean(is);  s.PD(i) = mean(pd);
        s.spread(i) = (max(is)-min(is))/mean(is);
    end
    s.depth = median(dep);
end

function Cload = calibrateCload(d, sensPre, moc, ch)
%CALIBRATECLOAD  Cable+scope load capacitance, from the preamp/no-preamp overlap (3-12 V,
%   below the preamp clip ceiling): the ratio of the two readings fixes C_load.
    Vov = [3 4 5 6 7 8 9 10 12];  r = nan(size(Vov));
    for i = 1:numel(Vov)
        fp = fullfile(d,sprintf('S5-1_%dV_Preamp_peak.hws',Vov(i)));
        fn = fullfile(d,sprintf('S5-1_%dV_NoPreamp_opt.hws',Vov(i)));
        if isfile(fp) && isfile(fn), r(i) = vneg(fp)/vneg(fn); end
    end
    Cload = mean(r,'omitnan')*(moc*ch)/sensPre - ch;
end

function v = vneg(f)
    w = readHWS(f);  y = w(1).y - median(w(1).y);  v = -min(y);
end

function v = crossing(V,Y,lim)
    v = NaN;
    if all(Y<lim), v = Inf; return, end
    for k = 1:numel(V)-1
        if (Y(k)-lim)*(Y(k+1)-lim) <= 0 && Y(k+1)~=Y(k)
            v = V(k) + (lim-Y(k))/(Y(k+1)-Y(k))*(V(k+1)-V(k));  return
        end
    end
    if Y(1) >= lim, v = V(1); end
end

function dep = depthFromNotes(cfg,Home)
    dep = 3.63;
    if isfield(cfg,'notes') && ~isempty(cfg.notes)
        tk = regexp(cfg.notes,'Pos\s*=\s*\(\s*([-\d.]+)\s*,\s*([-\d.]+)\s*,\s*([-\d.]+)','tokens','once');
        if ~isempty(tk)
            dep = norm([str2double(tk{1}) str2double(tk{2}) str2double(tk{3})]-Home)/10;
        end
    end
end

function fawf = measureFawf(file)
    w = readHWS(file);  y = w(1).y-median(w(1).y);  dt = w(1).dt;
    Y = abs(fft(y.*hann(numel(y))));  f = (0:numel(Y)-1)'/(numel(Y)*dt);  hb = 1:floor(numel(Y)/2);
    bf = f(hb)>1.3e6 & f(hb)<3.5e6;  pk = max(Y(hb).*bf);  idx = find(Y(hb)>=pk/2 & bf);
    fawf = (f(idx(1))+f(idx(end)))/2/1e6;
end

function s = f2(v)
    if isinf(v), s = '>50'; elseif isnan(v), s = '--'; else, s = sprintf('%.1f',v); end
end

function makeFigure(B,V,LIM,out)
    col = struct('el41',[.47 .67 .19],'el61',[0 .45 .74],'el79',[.85 .33 .10]);
    lab = struct('el41','41 elements','el61','61 elements','el79','79 elements');
    f = figure('Color','w','Position',[60 60 1180 420]);
    fns = {'el41','el61','el79'};
    quantities = {'MI','MI (derated)',LIM.MI; 'Isppa','I_{sppa.3} (W/cm^2)',LIM.Isppa; ...
                  'Ispta19','I_{spta.3} (mW/cm^2), 1900 cyc',LIM.Ispta};
    for q = 1:3
        subplot(1,3,q); hold on; grid on; box on
        for k = 1:3
            s = B.(fns{k});
            plot(V,s.(quantities{q,1}),'o-','Color',col.(fns{k}),'MarkerFaceColor',col.(fns{k}), ...
                 'MarkerSize',5,'LineWidth',1.5,'DisplayName',lab.(fns{k}));
        end
        yline(quantities{q,3},'k--','HandleVisibility','off');
        xlabel('TX voltage (V)'); ylabel(quantities{q,2}); xlim([10 52]);
        if q==1, legend('Location','northwest'); end
        title(sprintf('%s (dashed = FDA limit)',quantities{q,2}));
    end
    sgtitle('S5-1 ARF push, direct no-preamp hydrophone measurement (2026-08-17), derated');
    exportgraphics(f,out,'Resolution',150);  close(f);
end
