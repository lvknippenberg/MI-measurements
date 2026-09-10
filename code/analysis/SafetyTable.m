function SafetyTable()
%SAFETYTABLE  MI / Isppa.3 / Ispta.3 for the S5-1 configurations, from .hws.
%
%   Reproduces the safety table for:
%     1. Focused imaging, 30 V                          (scanned, PRF 52.2 Hz)
%     2. Focused push 79 el, 1900 cyc, max allowed V    (burst duty)
%     3. Focused push 61 el, 1900 cyc, max allowed V    (burst duty)
%     4. Focused push 61 el, 1500 cyc, max allowed V    (burst duty)
%     5. Focused push 41 el, 1500 cyc, 30 V   (ESTIMATED, not measured)
%
%   Method (identical to the rest of the analysis):
%     - depth = distance from Home (transducer centre) to the peak
%     - derating 0.3 dB/cm/MHz at the measured f_awf (-6 dB fundamental)
%     - Onda parallel-circuit sensitivity (HydrophoneTables.mat); 50 ohm data
%     - MI = pr.3 / sqrt(f_awf);  MI is linear in V
%     - Isppa.3 = PII.3 / PD  [W/cm2];  ~ V^2
%     - Ispta.3 = PII.3 * PRF [mW/cm2]
%         imaging PRF = 1/(71 lines * 270 us PRI)         = 52.2 Hz (scanned)
%         push PRF (burst) = 24 pushes /(1.2 s + 30 s off) = 0.769 Hz
%     - "max allowed V" = min voltage at which MI=1.9, Isppa=190 or Ispta=720
%     - pulse length scales Ispta (~cycles); MI and Isppa are unaffected
%     - 41 el is scaled from the 61/79 el fit (focal gain ~ N^p); MEASURE it
%       before use. Values within ~15%.
%
%   Caveat: the burst Ispta is the regulatory long-time average; the transient
%   temperature rise during each 1.2 s burst should be checked separately (TI).
%
%   Just run:  SafetyTable

    addpath(fileparts(mfilename('fullpath')));  % co-located analysis code
    base = miData('S5-1','2026-08-13_sessionA_preamp_50ohm');
    Home = [-82.21 105.14 -42.90];
    rho = 1000; c = 1500;
    LIM = [1.9 190 720];                       % MI, Isppa [W/cm2], Ispta [mW/cm2]
    PRF_img  = 1/(71*270e-6);                  % scanned imaging: 71 lines, 270 us PRI
    PRF_burst = 24/(1.2+30);                   % push: 20 Hz for 1.2 s, then >=30 s off

    % --- read each voltage sweep -> linear coefficients ---
    imgDir = fullfile(base,'focused_imaging_sweep');
    p79Dir = fullfile(base,'push_79el_sweep');
    p61Dir = fullfile(base,'push_61el_sweep');
    [kMIi,aIsi,PDi] = readsweep(imgDir,'S5-1_%dV_50ohm_pos_foc.hws',      [2 4 6 8 10 12 15 20], Home,rho,c);
    [kMI9,aIs9,PD9] = readsweep(p79Dir,'S5-1_%dV_50ohm_fund_push.hws',    [2 3 4 5 6 7 8],       Home,rho,c);
    [kMI6,aIs6,PD6] = readsweep(p61Dir,'S5-1_%dV_50ohm_fund_push_61el__opt.hws',[3 4 6 7 8],     Home,rho,c);

    % --- 41 el: scale from the 61/79 el focal-gain fit (kMI ~ N^p) ---
    N4=41; N6=61; N9=79;
    pexp = log(kMI6/kMI9)/log(N6/N9);
    kMI4 = kMI9*(N4/N9)^pexp;
    aIs4 = aIs9*(kMI4/kMI9)^2;                 % Isppa ~ pressure^2 ~ kMI^2

    % --- helper to evaluate a config at a voltage ---
    %   pdscale = actual pulse duration / measured pulse duration
    %   (1 for imaging and 1900-cyc pushes; 1500/1900 for the shorter push).
    row = @(V,kMI,aIs,PD,pdscale,PRF) [V, kMI*V, aIs*V^2, aIs*V^2*PD*pdscale*PRF*1000];
    % max allowed V = first limit reached
    maxV = @(kMI,aIs,PD,pdscale,PRF) min([1.9/kMI, sqrt(190/aIs), ...
                                      sqrt(720/(aIs*PD*pdscale*PRF*1000))]);
    s15 = 1500/1900;

    R = {}; nm = {};
    nm{end+1}='Focused imaging, 30 V';            R{end+1}=row(30,kMIi,aIsi,PDi,1,PRF_img); %#ok<*AGROW>
    V=maxV(kMI9,aIs9,PD9,1,PRF_burst);
    nm{end+1}='Push 79el, 1900cyc, max V';        R{end+1}=row(V,kMI9,aIs9,PD9,1,PRF_burst);
    V=maxV(kMI6,aIs6,PD6,1,PRF_burst);
    nm{end+1}='Push 61el, 1900cyc, max V';        R{end+1}=row(V,kMI6,aIs6,PD6,1,PRF_burst);
    V=maxV(kMI6,aIs6,PD6,s15,PRF_burst);
    nm{end+1}='Push 61el, 1500cyc, max V';        R{end+1}=row(V,kMI6,aIs6,PD6,s15,PRF_burst);
    nm{end+1}='Push 41el, 1500cyc, 30 V (est)';   R{end+1}=row(30,kMI4,aIs4,PD6,s15,PRF_burst);

    % --- print ---
    fprintf('\n%-30s | %5s | %5s | %8s | %9s\n','Configuration','V(V)','MI','Isppa','Ispta');
    fprintf('%-30s | %5s | %5s | %8s | %9s\n','','','','[W/cm2]','[mW/cm2]');
    fprintf('%s\n', repmat('-',1,72));
    for i=1:numel(nm)
        f = @(x,l) sprintf('%s%s', num2str(x,'%.1f'), tern(x>l));
        fprintf('%-30s | %5.1f | %s | %s | %s\n', nm{i}, R{i}(1), ...
            flag(R{i}(2),LIM(1),'%.2f'), flag(R{i}(3),LIM(2),'%.1f'), flag(R{i}(4),LIM(3),'%.1f'));
    end
    fprintf('%s\n', repmat('-',1,72));
    fprintf('%-30s | %5s | %5.2f | %8.1f | %9.1f\n','FDA LIMIT','', LIM(1),LIM(2),LIM(3));
    fprintf('\n(* = exceeds limit;  pushes are Isppa-limited;  Ispta uses the burst duty)\n');
end

% ===================================================================
function [kMI,aIs,PDref] = readsweep(dir,pat,Vs,Home,rho,c)
%READSWEEP  Linear coefficients from a voltage sweep: MI=kMI*V, Isppa=aIs*V^2.
    Vs = Vs(:); MI = nan(numel(Vs),1); Is = nan(numel(Vs),1); dep = NaN; PDref = NaN;
    fawf = measureFawf(fullfile(dir,sprintf(pat,Vs(1)))); sens = sensVperMPa(fawf);
    for i=1:numel(Vs)
        [w,cfg] = readHWS(fullfile(dir,sprintf(pat,Vs(i)))); yb = w(1).y-median(w(1).y); dt = w(1).dt;
        if isnan(dep)
            tk = regexp(cfg.notes,'Pos\s*=\s*\(\s*([-\d.]+)\s*,\s*([-\d.]+)\s*,\s*([-\d.]+)','tokens','once');
            dep = norm([str2double(tk{1}) str2double(tk{2}) str2double(tk{3})]-Home);
        end
        der = db2mag(0.3*dep/10*fawf); Ider = (1/der)^2;
        Vpk = -min(yb); p = yb/sens*1e6; I = p.^2/(rho*c); cR = rescale(cumsum(I),0,1);
        i10 = find(cR>0.1,1); i90 = find(cR>0.9,1); PD = 1.25*(i90-i10)*dt; PII = trapz(I(i10:i90))*dt;
        MI(i) = Vpk/sens/der/sqrt(fawf); Is(i) = PII*Ider/PD/1e4; PDref = PD;
    end
    m = (MI./Vs) >= 0.9*max(MI./Vs);          % unsaturated points
    kMI = sum(MI(m).*Vs(m))/sum(Vs(m).^2);    % through-origin linear fit
    aIs = mean(Is(m)./Vs(m).^2);
end

function fawf = measureFawf(file)
    w = readHWS(file); y = w(1).y-median(w(1).y); dt = w(1).dt;
    Y = abs(fft(y.*hann(numel(y)))); f = (0:numel(Y)-1)'/(numel(Y)*dt); hb = 1:floor(numel(Y)/2);
    bf = f(hb)>1.3e6 & f(hb)<3.5e6; pk = max(Y(hb).*bf); idx = find(Y(hb)>=pk/2 & bf);
    fawf = (f(idx(1))+f(idx(end)))/2/1e6;
end

function s = sensVperMPa(ft)
    persistent HT AT
    if isempty(HT)
        S = load(miRoot('calibration','HydrophoneTables.mat'), ...
                 'AmplifierTable','HydrophoneTable'); HT = S.HydrophoneTable; AT = S.AmplifierTable;
    end
    ih = find(HT.FREQ_MHz>=ft,1); sh = db2mag(HT{ih,"SENS_DB"})*1e6; ch = HT{ih,"CAP_PF"}*1e-12;
    ia = find(AT.FREQ_MHZ>=ft,1); ga = db2mag(AT{ia,"GAIN_DB"});     ca = AT{ia,"CAP_PF"}*1e-12;
    s = ga*sh*ch/(ca+ch)*1e6;
end

function s = flag(x,lim,fmt)
    s = sprintf(fmt,x); if x>lim, s = [s '*']; end
    s = sprintf('%8s', s);
end
function s = tern(b); if b, s='*'; else, s=''; end; end
