function CF = spatialAvgFactor(f_MHz, Fnum, aperture_um, c)
%SPATIALAVGFACTOR Finite-aperture (spatial-averaging) correction factor.
%   CF = spatialAvgFactor(f_MHz, Fnum, aperture_um, c) returns CF >= 1 such
%   that   p_true(peak) = CF * p_measured   for a circular hydrophone of
%   diameter aperture_um measuring the FOCAL PLANE of a focused aperture.
%
%   A finite tip averages the focal pressure over its face and therefore
%   UNDER-reads the true spatial-peak. This routine models the focal field
%   as the Airy/jinc pattern of a circular focusing aperture,
%       p(r)/p(0) = 2*J1(u)/u ,   u = pi*r / (lambda*Fnum),
%   and returns CF = 1 / <p(r)/p(0)> averaged over the hydrophone face.
%
%   Inputs:  f_MHz       acoustic working frequency [MHz]
%            Fnum        f-number of the transmit beam at the measurement plane
%            aperture_um hydrophone geometric aperture [um] (e.g. 400)
%            c           speed of sound [m/s] (default 1500)
%
%   CAVEATS: (1) first-order estimate only; (2) a linear array gives a
%   sinc (rectangular) rather than jinc focal pattern, and its elevation
%   focus differs from azimuth, so treat CF as an upper-ish bound; (3) valid
%   at/near the focal plane. For a rigorous value use IEC 62127-1 or a
%   directly measured beam width. The correction INCREASES pressure/MI.
%
%   Example:
%       CF = spatialAvgFactor(5.0, 2, 400);   % ~1.04

    if nargin < 4 || isempty(c), c = 1500; end
    lambda = c/(f_MHz*1e6)*1e3;           % mm
    a  = (aperture_um/1e3)/2;             % hydrophone radius [mm]
    kr = pi/(lambda*Fnum);

    % area-weighted average of the normalized focal amplitude over the face
    Msa = integral(@(r) jinc(kr.*r).*(2*r)/a^2, 0, a);
    CF  = 1/Msa;
end

function v = jinc(u)
    v = ones(size(u));
    nz = u ~= 0;
    v(nz) = 2*besselj(1,u(nz))./u(nz);
end
