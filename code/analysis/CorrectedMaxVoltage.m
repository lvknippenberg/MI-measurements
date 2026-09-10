function CorrectedMaxVoltage()
%CORRECTEDMAXVOLTAGE  Max allowed TX voltage per push config from DIRECT no-preamp data.
%
%   Independent cross-check of SafetyTable (which extrapolated a linear MI~V / quadratic
%   Isppa~V^2 fit of the low-voltage 2-8 V *with-preamp* sweep). Here we instead read the
%   direct no-preamp push captures (15-50 V, preamp removed so no clipping) for the
%   peak-optimised ("_opt") 41 / 61 / 79 element files and evaluate MI.3 / Isppa.3 / Ispta.3
%   at each measured voltage, then find where each FDA limit is actually crossed.
%
%   RESULT: the direct measurement REPRODUCES SafetyTable almost exactly (41el 32 V, 61el 23 V,
%   79el 20 V, all Isppa.3-limited). The low-voltage extrapolation is therefore validated: the
%   limits are crossed at 20-32 V, BELOW the region where nonlinear saturation / supply sag bite
%   (that only shows above ~40 V, e.g. 79el Vneg peaks ~40 V then dips at 50 V), so extrapolation
%   and direct measurement agree at the crossing voltages.
%
%   NB there are two no-preamp 79el datasets at the SAME position/depth: these peak-optimised
%   "_opt_79el.hws" single captures (WORST CASE, Vneg ~102 mV @40 V -> used here for the limit),
%   and the later "79 elements repeats\" set (TYPICAL, Vneg ~66 mV @40 V, monotonic, matches the
%   operator's multi-second observation -- see RepeatVariance.m). Safety uses the worst case, so
%   this script and the table use the "_opt" captures. Do NOT use the lower repeat values for a
%   limit -- they read the supply-sagged typical, not the peak.
%
%   Conditions (Notes.txt): push at depth 36.3 mm, 1 MOhm scope, no preamp, 1900 cycles,
%   f_awf measured per aperture; burst duty 20 Hz for 1.2 s then >=30 s off -> 0.769 Hz.

    addpath(fileparts(mfilename('fullpath')));
    d = miData('S5-1','2026-08-17_sessionB_preamp_vs_nopreamp');
    repo = miRoot('calibration');
    S = load(fullfile(repo,'HydrophoneTables.mat')); HT = S.HydrophoneTable;
    rho = 1000; c = 1500; Cload = 123.3e-12; Home = [55.8 105.84 -17.3];
    PRF_burst = 24/(1.2+30);  s15 = 1500/1900;  LIM = [1.9 190 720];
    V = [15 20 25 30 35 40 45 50];
    apert = {41,'_41el'; 61,'';  79,'_79el'};

    fprintf('\n%-8s %-6s %-6s %-6s %-6s | %-8s\n','config','MI=1.9','Is=190','Ita15','Ita19','max V');
    fprintf('%s\n',repmat('-',1,52));
    results = struct();
    for a = 1:size(apert,1)
        N = apert{a,1}; suf = apert{a,2};
        MI = nan(size(V)); Is = nan(size(V)); PD = nan(size(V));
        fawf = measureFawf(fullfile(d,sprintf('S5-1_%dV_NoPreamp_opt%s.hws',V(1),suf)));
        ch = interp1(HT.FREQ_MHz,HT.CAP_PF,fawf)*1e-12;
        moc = db2mag(interp1(HT.FREQ_MHz,HT.SENS_DB,fawf))*1e6;      % V/Pa (open circuit)
        sens = moc*ch/(ch+Cload);                                    % V/Pa (loaded, no gain)
        for i = 1:numel(V)
            f = fullfile(d,sprintf('S5-1_%dV_NoPreamp_opt%s.hws',V(i),suf));
            [w,cfg] = readHWS(f); y = w(1).y-median(w(1).y); dt = w(1).dt;
            dep = depthFromNotes(cfg,Home);                          % cm
            der = db2mag(0.3*dep*fawf); Ider = (1/der)^2;
            p = y/sens;                                              % Pa
            pr = -min(p);                                            % Pa (peak rarefac.)
            MI(i) = (pr/1e6)/der/sqrt(fawf);
            I = p.^2/(rho*c); cR = rescale(cumsum(I),0,1);
            i10 = find(cR>0.1,1); i90 = find(cR>0.9,1);
            PD(i) = 1.25*(i90-i10)*dt; PII = trapz(I(i10:i90))*dt;
            Is(i) = PII*Ider/PD(i)/1e4;                              % W/cm^2 derated
        end
        PDm = median(PD);
        Ita19 = Is.*PD*PRF_burst*1000;                              % mW/cm^2, 1900 cyc
        Ita15 = Is.*PD*s15*PRF_burst*1000;                          % mW/cm^2, 1500 cyc
        vMI  = crossing(V,MI,LIM(1));
        vIs  = crossing(V,Is,LIM(2));
        vI15 = crossing(V,Ita15,LIM(3));
        vI19 = crossing(V,Ita19,LIM(3));
        maxV15 = min([vMI vIs vI15]); maxV19 = min([vMI vIs vI19]);
        fprintf('%2del     %-6s %-6s %-6s %-6s | 1500:%.0f 1900:%.0f\n', N, ...
            f2(vMI), f2(vIs), f2(vI15), f2(vI19), maxV15, maxV19);
        results.(sprintf('el%d',N)) = struct('V',V,'MI',MI,'Isppa',Is,'Ita15',Ita15,'Ita19',Ita19,...
            'vMI',vMI,'vIs',vIs,'maxV15',maxV15,'maxV19',maxV19,'fawf',fawf,'depth_cm',median(depthAll(d,suf,V,Home)));
    end
    fprintf(['\nDirect no-preamp measurement (nonlinear saturation captured); supersedes the\n' ...
             'low-voltage-extrapolated SafetyTable. 1500 vs 1900 cyc differ only via Ispta.\n']);
    assignin('base','maxVresults',results);
end

% --- helpers ----------------------------------------------------------------
function v = crossing(V,Y,lim)
%CROSSING  First voltage where monotone-ish Y crosses lim (linear interp); NaN if never.
    v = NaN;
    if all(Y<lim), v = Inf; return; end          % never reached in range -> ">max"
    for k = 1:numel(V)-1
        if (Y(k)-lim)*(Y(k+1)-lim) <= 0 && Y(k+1)~=Y(k)
            v = V(k) + (lim-Y(k))/(Y(k+1)-Y(k))*(V(k+1)-V(k)); return
        end
    end
    if Y(1) >= lim, v = V(1); end                 % already over at lowest V
end

function dep = depthFromNotes(cfg,Home)
    dep = 3.63;  % cm default
    if isfield(cfg,'notes') && ~isempty(cfg.notes)
        tk = regexp(cfg.notes,'Pos\s*=\s*\(\s*([-\d.]+)\s*,\s*([-\d.]+)\s*,\s*([-\d.]+)','tokens','once');
        if ~isempty(tk)
            dep = norm([str2double(tk{1}) str2double(tk{2}) str2double(tk{3})]-Home)/10; % cm
        end
    end
end

function ds = depthAll(d,suf,V,Home)
    ds = nan(size(V));
    for i=1:numel(V)
        try; [~,cfg]=readHWS(fullfile(d,sprintf('S5-1_%dV_NoPreamp_opt%s.hws',V(i),suf)));
             ds(i)=depthFromNotes(cfg,Home); catch; end
    end
end

function fawf = measureFawf(file)
    w = readHWS(file); y = w(1).y-median(w(1).y); dt = w(1).dt;
    Y = abs(fft(y.*hann(numel(y)))); f = (0:numel(Y)-1)'/(numel(Y)*dt); hb = 1:floor(numel(Y)/2);
    bf = f(hb)>1.3e6 & f(hb)<3.5e6; pk = max(Y(hb).*bf); idx = find(Y(hb)>=pk/2 & bf);
    fawf = (f(idx(1))+f(idx(end)))/2/1e6;
end

function s = f2(v)
    if isinf(v), s = '>50'; elseif isnan(v), s = '--'; else, s = sprintf('%.1f',v); end
end
