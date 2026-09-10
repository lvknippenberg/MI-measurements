function FreqDomainSensitivity()
%FREQDOMAINSENSITIVITY  Fixed-f0 vs frequency-domain hydrophone sensitivity.
%
%   As TX voltage rises the S5-1 push waveform develops stronger harmonics
%   (2nd harmonic at ~4.5 MHz). The default safety chain converts voltage to
%   pressure with the hydrophone sensitivity taken at the fundamental only
%   (p(t) = y(t)/S(f0)). This script instead deconvolves in the frequency
%   domain, applying the *frequency-dependent* loaded sensitivity
%       S(f) = M_oc(f) * C_h(f) / (C_h(f) + C_load)          [no-preamp path]
%   across a 1-12 MHz passband, so each harmonic is weighted by the sensitivity
%   at its own frequency. It reports pr / MI / I_sppa.3 both ways.
%
%   This mirrors the LaTeX reference chain ("Ultrasound safety indices"):
%     V_DAQ -> (remove offset, /G) -> FFT -> P(f)=V(f)/S(f) -> IFFT -> p(t)
%     -> pr->MI ; I(t)=p^2/(rho c) -> PII -> I_sppa=PII/tau_p ; I_spta=I_sppa*tau*PRF
%
%   CAVEAT (phase): a *complete* complex deconvolution needs the hydrophone
%   PHASE vs frequency. The Onda HGL-0400 calibration is magnitude-only, so
%   S(f) here is real (zero-phase) and the measured phase is retained. This
%   corrects the harmonic *amplitudes* but not the hydrophone's phase
%   distortion; treat the frequency-domain column as a magnitude-only refinement
%   until a phase-calibrated sensitivity is available.

    addpath(fileparts(mfilename('fullpath')));
    repo = miRoot('calibration');           % ...\Mechanical index
    d = miData('S5-1','2026-08-17_sessionB_preamp_vs_nopreamp','79el_repeats');

    S  = load(fullfile(repo,'HydrophoneTables.mat'));  HT = S.HydrophoneTable;
    rho = 1000; c = 1500; Cload = 123.3; f0 = 2.23; depth = 3.63;   % cm

    McF = @(fq) db2mag(interp1(HT.FREQ_MHz, HT.SENS_DB, fq,'linear','extrap'))*1e6;   % V/Pa
    ChF = @(fq) interp1(HT.FREQ_MHz, HT.CAP_PF,  fq,'linear','extrap');               % pF
    Sf  = @(fq) McF(fq).*ChF(fq)./(ChF(fq)+Cload);                                    % loaded V/Pa
    Sfix = Sf(f0);
    der  = db2mag(0.3*depth*f0);  Ider = (1/der)^2;                                   % derating

    V = [20 25 30 35 40 45 50];
    fprintf(['\nS(f): S(2.25)=%.3e  S(4.5)=%.3e  S(6.75)=%.3e V/Pa\n\n'], Sf(2.25),Sf(4.5),Sf(6.75));
    fprintf(' TX | H2/H1 | pr_fix pr_fd | dPr%% | MI_fix MI_fd | Isppa_fix Isppa_fd | dI%%\n');
    for iv = 1:numel(V)
        f = fullfile(d, sprintf('S5-1_%dV_NoPreamp_opt_79el.hws', V(iv)));
        raw = double(h5read(f,'/wfm_group0/vectors/vector0/data'));
        sc  = h5read(f,'/wfm_group0/axes/axis1/scale_coef');
        wf  = readHWS(f);  dt = wf(1).dt;                          % s/sample
        y   = sc(1)+sc(2)*raw;  y = y-median(y);  N = numel(y);

        % harmonic ratio (2nd/1st) from a windowed spectrum
        Yw = abs(fft(y.*hann(N)));  fa = (0:N-1)'/(N*dt);
        H2H1 = max(Yw(fa>4e6 & fa<5e6)) / max(Yw(fa>1.8e6 & fa<2.7e6));

        % fixed: time-domain scaling
        p_fix = y / Sfix;

        % frequency-domain magnitude deconvolution over 1-12 MHz
        fs = 1/dt; fsig = (0:N-1)'/(N*dt); fsig(fsig>fs/2) = fsig(fsig>fs/2)-fs;
        faM = abs(fsig)/1e6;  Smag = Sf(faM);  m = faM>=1 & faM<=12;
        Vf = fft(y);  P = zeros(N,1);  P(m) = Vf(m)./Smag(m);  p_fd = real(ifft(P));

        prf = -min(p_fix)/1e6;  prd = -min(p_fd)/1e6;
        MIf = (prf/der)/sqrt(f0);  MId = (prd/der)/sqrt(f0);
        Isf = isppa(p_fix,rho,c,dt,Ider);  Isd = isppa(p_fd,rho,c,dt,Ider);
        fprintf('%3d | %.2f  | %5.2f %5.2f | %+4.1f | %5.2f %5.2f | %8.1f %8.1f | %+4.1f\n', ...
            V(iv), H2H1, prf, prd, 100*(prd-prf)/prf, MIf, MId, Isf, Isd, 100*(Isd-Isf)/Isf);
    end
    fprintf(['\nFrequency-domain correction raises pr/MI/I_sppa by up to ~8%% ' ...
             '(conservative direction),\ngrowing with harmonic content. See phase CAVEAT in the header.\n']);
end

function Is = isppa(p, rho, c, dt, Ider)
    I  = p.^2/(rho*c);
    cR = rescale(cumsum(I),0,1);
    i10 = find(cR>0.1,1);  i90 = find(cR>0.9,1);
    PD  = (i90-i10)*dt*1.25;
    PII = trapz(I(i10:i90))*dt;
    Is  = PII*Ider/PD/1e4;                      % W/cm^2, derated
end
