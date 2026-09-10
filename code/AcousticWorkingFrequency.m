function [fawf, info] = AcousticWorkingFrequency(y, dt, bandMHz)
%ACOUSTICWORKINGFREQUENCY  IEC 62127 acoustic working frequency [MHz].
%
%   fawf = AcousticWorkingFrequency(y, dt, bandMHz) returns the arithmetic
%   centre of the -6 dB band of the FUNDAMENTAL lobe of the Hann-windowed
%   spectrum of the waveform y sampled at interval dt, searched inside
%   bandMHz = [flo fhi].
%
%   fawf = AcousticWorkingFrequency(hwsFile, [], bandMHz) does the same for a
%   capture file, reading it with readHWS and removing the median baseline.
%
%   Use a LOW-drive capture: at high transmit voltage the harmonics grow and
%   bias the centre upward.
%
%   info returns the -6 dB edges f1, f2 and the peak frequency.

    if (ischar(y) || isstring(y)) && (nargin < 2 || isempty(dt))
        w  = readHWS(char(y));
        dt = w(1).dt;
        y  = w(1).y - median(w(1).y);
    end
    y = y(:);
    n = numel(y);
    win = 0.5 - 0.5*cos(2*pi*(0:n-1)'/(n-1));       % hann, symmetric (no toolbox)
    Y = abs(fft(y.*win));
    f = (0:n-1)'/(n*dt);
    hb = 1:floor(n/2);
    f = f(hb);  Y = Y(hb);

    in = f > bandMHz(1)*1e6 & f < bandMHz(2)*1e6;
    [pk, ipk] = max(Y.*in);
    idx  = find(Y >= pk/2 & in);                    % -6 dB = half amplitude
    fawf = (f(idx(1)) + f(idx(end)))/2 / 1e6;

    info = struct('f1_MHz',f(idx(1))/1e6, 'f2_MHz',f(idx(end))/1e6, ...
                  'fpeak_MHz',f(ipk)/1e6, 'bandMHz',bandMHz);
end
