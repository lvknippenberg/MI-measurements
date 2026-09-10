function [prNeg_MPa, p, info] = voltageToPressureBroadband(t, v, tables, opts)
%VOLTAGETOPRESSUREBROADBAND Frequency-domain deconvolution to pressure.
%   [prNeg_MPa, p] = voltageToPressureBroadband(t, v, tables) converts a
%   measured hydrophone voltage waveform v(t) [V] into the acoustic pressure
%   waveform p(t) [Pa] using the FULL complex transfer function of the
%   hydrophone + preamp chain (magnitude AND amplifier phase), rather than a
%   single-frequency magnitude sensitivity. Returns the peak-rarefactional
%   (peak-negative) pressure in MPa.
%
%   System transfer function [V/Pa]:
%       H(f) = M_hyd(f) * G_amp(f) e^{i*phi_amp(f)} * C_h(f)/(C_h(f)+C_amp(f))
%   with M_hyd the hydrophone open-circuit sensitivity (magnitude; phase
%   assumed ~linear/flat in-band, so it does not affect the peak), G_amp and
%   phi_amp the preamp gain magnitude and phase, and the capacitive divider.
%
%   'tables' is the struct loaded from HydrophoneTables.mat (fields
%   HydrophoneTable, AmplifierTable). 'opts' (optional) fields:
%       flo, fhi   passband [MHz], default 1 and 18 (hydrophone cal is
%                  valid 1-20 MHz; band-limiting avoids out-of-band noise
%                  blow-up from dividing by a small sensitivity).
%       taper      cosine taper width [MHz] at each band edge, default 1.
%
%   The deconvolution is linear, so a scope-impedance load correction can be
%   applied either to v beforehand or to prNeg afterward.

    if nargin < 4, opts = struct; end
    flo   = getf(opts,'flo',1);
    fhi   = getf(opts,'fhi',18);
    taper = getf(opts,'taper',1);

    v = v(:); t = t(:);
    N = numel(v); dt = t(2)-t(1); fs = 1/dt;

    % signed frequency axis (so phase can be made odd -> real ifft)
    f = (0:N-1)'/(N*dt);
    fsigned = f; fsigned(f > fs/2) = fsigned(f > fs/2) - fs;
    fa = abs(fsigned)/1e6;                                  % |f| in MHz

    HT = tables.HydrophoneTable; AT = tables.AmplifierTable;
    Mh  = interp1(HT.FREQ_MHz, HT.SENS_VPERPA,       fa, 'linear', NaN);      % V/Pa
    Ch  = interp1(HT.FREQ_MHz, HT.CAP_PF,            fa, 'linear', 'extrap')*1e-12;
    Ga  = interp1(AT.FREQ_MHZ, db2mag(AT.GAIN_DB),  fa, 'linear', NaN);       % lin
    Pha = interp1(AT.FREQ_MHZ, AT.PHASE_DEG,        fa, 'linear', NaN)*pi/180;
    Ca  = interp1(AT.FREQ_MHZ, AT.CAP_PF,           fa, 'linear', 'extrap')*1e-12;

    div  = Ch./(Ch+Ca);
    Hmag = Mh .* Ga .* div;                 % even in f
    phi  = Pha .* sign(fsigned);            % odd in f  -> Hermitian H
    H    = Hmag .* exp(1i*phi);             % V/Pa

    % passband window (raised-cosine tapered), even in |f|
    W = bandWindow(fa, flo, fhi, taper);

    valid = isfinite(H) & (Hmag > 0) & (W > 0);
    V = fft(v);
    P = zeros(N,1);
    P(valid) = V(valid) ./ H(valid) .* W(valid);          % Pa spectrum
    p = real(ifft(P));                                    % Pa

    prNeg_MPa = -min(p)/1e6;
    if nargout > 2
        info = struct('flo',flo,'fhi',fhi,'taper',taper, ...
                      'sens_at_peakfreq', Hmag(2));
    end
end

% -------------------------------------------------------------------
function W = bandWindow(fa, flo, fhi, taper)
    W = zeros(size(fa));
    W(fa>=flo & fa<=fhi) = 1;
    lo = fa>=(flo-taper) & fa<flo;
    W(lo) = 0.5*(1 - cos(pi*(fa(lo)-(flo-taper))/taper));
    hi = fa>fhi & fa<=(fhi+taper);
    W(hi) = 0.5*(1 + cos(pi*(fa(hi)-fhi)/taper));
end

function v = getf(s,f,d)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
