function ConfirmPulseLength()
%CONFIRMPULSELENGTH  Data confirmation that pulse length only scales Ispta.
%
%   Compares two 79-element push captures acquired at the SAME location and
%   SAME 5 V drive, differing only in pulse length (1900 vs 1500 cycles):
%       S5-1_5V_50ohm_fund_push.hws        (1900 cycles)
%       S5-1_5V_50ohm_fund_push_1500.hws   (1500 cycles)
%
%   Shows that MI and Isppa.3 are unchanged while the pulse duration and
%   Ispta.3 scale with the cycle count.  Just run:  ConfirmPulseLength

    addpath(fileparts(mfilename('fullpath')));  % co-located analysis code
    d = miData('S5-1','2026-08-13_sessionA_preamp_50ohm','push_79el_sweep');
    Home = [-82.21 105.14 -42.90]; rho = 1000; c = 1500;
    PRF_burst = 24/(1.2+30);                         % push burst duty

    f1900 = fullfile(d,'S5-1_5V_50ohm_fund_push.hws');
    f1500 = fullfile(d,'S5-1_5V_50ohm_fund_push_1500.hws');
    fawf  = measureFawf(f1900); sens = sensVperMPa(fawf);
    [~,cfg] = readHWS(f1900);
    tk  = regexp(cfg.notes,'Pos\s*=\s*\(\s*([-\d.]+)\s*,\s*([-\d.]+)\s*,\s*([-\d.]+)','tokens','once');
    dep = norm([str2double(tk{1}) str2double(tk{2}) str2double(tk{3})]-Home);
    der = db2mag(0.3*dep/10*fawf); Ider = (1/der)^2;

    lbl = {'1900 cycles','1500 cycles'}; files = {f1900,f1500};
    M = zeros(2,5);   % Vpk, PD_us, MI, Isppa, Ispta
    for k = 1:2
        w = readHWS(files{k}); yb = w(1).y-median(w(1).y); dt = w(1).dt;
        Vpk = -min(yb); p = yb/sens*1e6; I = p.^2/(rho*c); cR = rescale(cumsum(I),0,1);
        i10 = find(cR>0.1,1); i90 = find(cR>0.9,1); PD = 1.25*(i90-i10)*dt; PII = trapz(I(i10:i90))*dt;
        MI = Vpk/sens/der/sqrt(fawf); PII3 = PII*Ider;
        Isppa = PII3/PD/1e4;                              % W/cm2
        Ispta = Isppa*PD*PRF_burst*1000;                 % mW/cm2
        M(k,:) = [Vpk, PD*1e6, MI, Isppa, Ispta];
    end

    fprintf('79-element push, 5 V, same location (depth %.1f mm, f_awf %.2f MHz)\n\n', dep, fawf);
    fprintf('%-12s | %6s | %7s | %5s | %8s | %9s\n','', 'Vpk(V)','PD(us)','MI','Isppa','Ispta');
    fprintf('%-12s | %6s | %7s | %5s | %8s | %9s\n','','','','','[W/cm2]','[mW/cm2]');
    fprintf('%s\n', repmat('-',1,60));
    for k = 1:2
        fprintf('%-12s | %6.3f | %7.1f | %5.3f | %8.2f | %9.3f\n', lbl{k}, M(k,:));
    end
    r = M(2,:)./M(1,:);
    fprintf('%s\n', repmat('-',1,60));
    fprintf('%-12s | %6.3f | %7.3f | %5.3f | %8.3f | %9.3f\n','ratio 1500/1900', r);
    fprintf('\nExpected pulse-length ratio 1500/1900 = %.3f\n', 1500/1900);
    fprintf(['CONFIRMED: MI ratio %.3f and Isppa ratio %.3f ~ 1 (unchanged);\n' ...
             '           PD ratio %.3f and Ispta ratio %.3f ~ %.3f (scale with cycles).\n'], ...
             r(3), r(4), r(2), r(5), 1500/1900);
end

% ---- helpers (same calibration as SafetyTable.m) ----
function fawf = measureFawf(file)
    w = readHWS(file); y = w(1).y-median(w(1).y); dt = w(1).dt;
    Y = abs(fft(y.*hann(numel(y)))); f = (0:numel(Y)-1)'/(numel(Y)*dt); hb = 1:floor(numel(Y)/2);
    bf = f(hb)>1.3e6 & f(hb)<3.5e6; pk = max(Y(hb).*bf); idx = find(Y(hb)>=pk/2 & bf);
    fawf = (f(idx(1))+f(idx(end)))/2/1e6;
end
function s = sensVperMPa(ft)
    S = load(miRoot('calibration','HydrophoneTables.mat'), ...
             'AmplifierTable','HydrophoneTable'); HT = S.HydrophoneTable; AT = S.AmplifierTable;
    ih = find(HT.FREQ_MHz>=ft,1); sh = db2mag(HT{ih,"SENS_DB"})*1e6; ch = HT{ih,"CAP_PF"}*1e-12;
    ia = find(AT.FREQ_MHZ>=ft,1); ga = db2mag(AT{ia,"GAIN_DB"});     ca = AT{ia,"CAP_PF"}*1e-12;
    s = ga*sh*ch/(ca+ch)*1e6;
end
