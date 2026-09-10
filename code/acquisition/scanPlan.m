function plan = scanPlan(f0, focusDepth, options)
%SCANPLAN Non-uniform 3-orthogonal-line hydrophone scan sized to the beam.
%   plan = scanPlan(f0, focusDepth) builds axial (z), azimuth (x) and
%   elevation (y) sample lines through the focus, with fine sampling around
%   the peak and coarse sampling in the plateau/wings, sized from the
%   diffraction scales of the beam (lambda, lateral FWHM, depth of focus).
%
%   Required inputs:
%       f0          - centre frequency [Hz]           (e.g. 3e6)
%       focusDepth  - geometric focal depth [mm]      (e.g. 40)
%
%   Name-value options:
%       Fnum          - f-number (focusDepth/aperture). If omitted, it is
%                       computed from 'aperture'; if that is also omitted a
%                       default of 2 is assumed (with a warning).
%       aperture      - active aperture [mm] (used only if Fnum not given)
%       elevationFnum - separate f-number for the elevation (lens) focus.
%                       Phased/linear arrays usually focus weaker in
%                       elevation, so its FWHM is larger. Default = Fnum.
%       c             - sound speed [m/s]. Default 1480 (water).
%       latRangeFWHM  - lateral half-range in units of FWHM. Default 2.5.
%       latFineFrac   - fine lateral step as a fraction of FWHM. Default 0.25.
%       latCoarseFrac - coarse lateral step as a fraction of FWHM. Default 0.5.
%       axRangeDOF    - axial half-range in units of DOF. Default 0.9.
%       axProxExtra   - extra proximal axial range (fraction of DOF) to catch
%                       the nonlinear proximal shift of the peak. Default 0.3.
%       axFineStep    - fine axial step [mm] near the max. Default 0.5.
%       axCoarseStep  - coarse axial step [mm] in the plateau. Default 1.5.
%       gridStep      - round all coordinates to this [mm] (stage resolution).
%                       Default 0.01.
%       lines         - which orthogonal lines to emit, any subset of
%                       ["axial" "azimuth" "elevation"]. Default all three.
%                       Use ["axial" "azimuth"] for a minimal Tier-1
%                       validation cross.
%       anchorOnly    - if true, return a single point at the centre - the
%                       Tier-0 worst-case anchor. Default false.
%       center        - [x y z] mm to centre the scan on. Default
%                       [0 0 focusDepth]. Pass the OBSERVED peak (found by a
%                       coarse search) so the fine cross sits on the true
%                       maximum, not the displaced geometric focus.
%       makePlot      - scatter3 the plan. Default false.
%
%   Measurement tiers (see below) trade points for reliance on simulation:
%       Tier 0 (per sequence):  scanPlan(f0,z,'anchorOnly',true)       1 pt
%       Tier 1 (per beam type):  scanPlan(f0,z,'Fnum',F, ...
%                                   'lines',["axial" "azimuth"])       ~20 pts
%       Tier 2 (full):           scanPlan(f0,z,'Fnum',F)               ~60 pts
%
%   Output struct 'plan':
%       .lambda_mm, .lat_FWHM_mm, .elev_FWHM_mm, .axial_DOF_mm, .Fnum
%       .x, .y, .z            - 1-D sample coordinates of each line [mm]
%       .axial, .azimuth, .elevation - N-by-3 [x y z] point lists per line
%       .points              - unique N-by-3 [x y z] of the whole scan [mm]
%       .tof_s               - time-of-flight to each point (|r|/c) [s],
%                              i.e. the trigger delay to place the capture
%                              window on the pulse arrival
%       .nPoints             - number of unique measurement locations
%
%   Coordinates: x = azimuth, y = elevation, z = depth. Lateral centre is
%   (x,y) = (0,0); the axial line runs along x=y=0; the lateral lines sit at
%   z = focusDepth. Feed .points to your positioner and .tof_s to the scope
%   trigger delay.
%
%   Example:
%       p = scanPlan(3e6, 40, 'Fnum', 2, 'makePlot', true);
%       % p.points -> ~50 locations to visit; p.tof_s -> per-point delay

    arguments
        f0            (1,1) double {mustBePositive}
        focusDepth    (1,1) double {mustBePositive}
        options.Fnum          (1,1) double = NaN
        options.aperture      (1,1) double = NaN
        options.elevationFnum (1,1) double = NaN
        options.c             (1,1) double {mustBePositive} = 1480
        options.latRangeFWHM  (1,1) double {mustBePositive} = 2.5
        options.latFineFrac   (1,1) double {mustBePositive} = 0.25
        options.latCoarseFrac (1,1) double {mustBePositive} = 0.5
        options.axRangeDOF    (1,1) double {mustBeNonnegative} = 0.9
        options.axProxExtra   (1,1) double {mustBeNonnegative} = 0.3
        options.axFineStep    (1,1) double {mustBePositive} = 0.5
        options.axCoarseStep  (1,1) double {mustBePositive} = 1.5
        options.gridStep      (1,1) double {mustBePositive} = 0.01
        options.lines         (1,:) string {mustBeMember(options.lines, ...
                                  ["axial","azimuth","elevation"])} = ...
                                  ["axial","azimuth","elevation"]
        options.anchorOnly    (1,1) logical = false
        options.center        (1,3) double = [NaN NaN NaN]
        options.makePlot      (1,1) logical = false
    end

    % --- f-number -------------------------------------------------
    Fnum = options.Fnum;
    if isnan(Fnum)
        if ~isnan(options.aperture)
            Fnum = focusDepth / options.aperture;
        else
            Fnum = 2;
            warning('scanPlan:defaultFnum', ...
                'No Fnum or aperture given; assuming Fnum = 2.');
        end
    end
    elevFnum = options.elevationFnum;
    if isnan(elevFnum); elevFnum = Fnum; end

    % --- diffraction scales --------------------------------------
    lambda  = options.c / f0 * 1e3;          % mm
    latFWHM = 1.0 * Fnum     * lambda;        % -6 dB lateral beamwidth (mm)
    elvFWHM = 1.0 * elevFnum * lambda;        % elevation beamwidth (mm)
    DOF     = 7.1 * Fnum^2   * lambda;        % -6 dB depth of focus (mm)

    % --- lateral lines (symmetric about 0) -----------------------
    xrel = symAxis(1.5*latFWHM, options.latFineFrac*latFWHM, ...
                   options.latRangeFWHM*latFWHM, options.latCoarseFrac*latFWHM);
    yrel = symAxis(1.5*elvFWHM, options.latFineFrac*elvFWHM, ...
                   options.latRangeFWHM*elvFWHM, options.latCoarseFrac*elvFWHM);

    % --- axial line (asymmetric: extra reach on the proximal side) ---
    distHalf = options.axRangeDOF*DOF;
    proxHalf = options.axRangeDOF*DOF + options.axProxExtra*DOF;
    zdist =  halfAxis(0.3*DOF, options.axFineStep, distHalf, options.axCoarseStep);
    zprox =  halfAxis(0.3*DOF, options.axFineStep, proxHalf, options.axCoarseStep);
    zrel  = unique([-fliplr(zprox), zdist]);

    % --- scan centre: the observed peak if given, else the focus ---
    %   Locate the peak with a coarse search first, then pass its
    %   coordinates as 'center' so the fine cross sits on the true maximum
    %   rather than the (displaced) geometric focus.
    ctr = options.center;
    if any(isnan(ctr)); ctr = [0, 0, focusDepth]; end

    % --- absolute coordinates, snapped to the stage grid ---------
    snap = @(v) round(v / options.gridStep) * options.gridStep;
    xabs = unique(snap(ctr(1) + xrel));
    yabs = unique(snap(ctr(2) + yrel));
    zabs = unique(snap(ctr(3) + zrel));

    % --- assemble the requested lines ----------------------------
    %   'anchorOnly' -> just the peak point (Tier 0). 'lines' selects any
    %   subset, e.g. ["axial","azimuth"] for a Tier-1 validation cross.
    wantAx = any(options.lines == "axial")     && ~options.anchorOnly;
    wantAz = any(options.lines == "azimuth")   && ~options.anchorOnly;
    wantEl = any(options.lines == "elevation") && ~options.anchorOnly;
    nAx = numel(zabs); nAz = numel(xabs); nEl = numel(yabs);
    axial     = [ctr(1)*ones(nAx,1), ctr(2)*ones(nAx,1), zabs(:)];
    azimuth   = [xabs(:),            ctr(2)*ones(nAz,1), ctr(3)*ones(nAz,1)];
    elevation = [ctr(1)*ones(nEl,1), yabs(:),            ctr(3)*ones(nEl,1)];

    parts = zeros(0,3);
    if wantAx; parts = [parts; axial];     else; axial     = zeros(0,3); end
    if wantAz; parts = [parts; azimuth];   else; azimuth   = zeros(0,3); end
    if wantEl; parts = [parts; elevation]; else; elevation = zeros(0,3); end
    if options.anchorOnly || isempty(parts)
        parts = snap(ctr);                 % single worst-case anchor
    end

    pts = unique(parts, 'rows');
    tof = vecnorm(pts, 2, 2) * 1e-3 / options.c;   % s

    % --- pack ----------------------------------------------------
    plan = struct();
    plan.Fnum          = Fnum;
    plan.lambda_mm     = lambda;
    plan.lat_FWHM_mm   = latFWHM;
    plan.elev_FWHM_mm  = elvFWHM;
    plan.axial_DOF_mm  = DOF;
    plan.center = ctr;
    plan.x = xabs(:); plan.y = yabs(:); plan.z = zabs(:);
    plan.axial = axial; plan.azimuth = azimuth; plan.elevation = elevation;
    plan.points = pts;
    plan.tof_s  = tof;
    plan.nPoints = size(pts, 1);

    % --- report --------------------------------------------------
    fprintf('scanPlan @ %.2f MHz, focus %.1f mm, F#=%.2f (elev F#=%.2f), c=%g m/s\n', ...
        f0/1e6, focusDepth, Fnum, elevFnum, options.c);
    fprintf('  centre (x,y,z) = (%.2f, %.2f, %.2f) mm\n', ctr(1), ctr(2), ctr(3));
    fprintf('  lambda=%.3f mm | lateral FWHM=%.2f mm | elevation FWHM=%.2f mm | DOF=%.1f mm\n', ...
        lambda, latFWHM, elvFWHM, DOF);
    if options.anchorOnly
        fprintf('  anchor only: 1 pt at centre (%.2f, %.2f, %.2f) mm\n', ctr(1), ctr(2), ctr(3));
    else
        if wantAx; fprintf('  axial:     %2d pts over z=[%.1f, %.1f] mm\n', ...
                nAx, zabs(1), zabs(end)); end
        if wantAz; fprintf('  azimuth:   %2d pts over x=[%.2f, %.2f] mm (step >=%.2f mm)\n', ...
                nAz, xabs(1), xabs(end), options.latFineFrac*latFWHM); end
        if wantEl; fprintf('  elevation: %2d pts over y=[%.2f, %.2f] mm (step >=%.2f mm)\n', ...
                nEl, yabs(1), yabs(end), options.latFineFrac*elvFWHM); end
    end
    fprintf('  total unique locations: %d | trigger delay range %.1f-%.1f us\n', ...
        plan.nPoints, min(tof)*1e6, max(tof)*1e6);

    if options.makePlot
        figure('Name', 'scanPlan');
        scatter3(pts(:,1), pts(:,3), pts(:,2), 28, tof*1e6, 'filled');
        xlabel('x azimuth [mm]'); ylabel('z depth [mm]'); zlabel('y elevation [mm]');
        cb = colorbar; cb.Label.String = 'trigger delay [\mus]';
        title('Hydrophone scan plan'); axis equal; grid on; view(35, 20);
    end
end

% ===================================================================
function p = halfAxis(fineHalf, fineStep, fullHalf, coarseStep)
%HALFAXIS Non-negative samples 0..fullHalf: fine to fineHalf, then coarse.
    if fullHalf <= fineHalf
        p = 0:fineStep:fullHalf;
    else
        core = 0:fineStep:fineHalf;
        wing = (fineHalf + coarseStep):coarseStep:fullHalf;
        p = [core, wing];
    end
    p = unique(p);
end

% ===================================================================
function a = symAxis(fineHalf, fineStep, fullHalf, coarseStep)
%SYMAXIS Symmetric non-uniform axis about 0.
    pos = halfAxis(fineHalf, fineStep, fullHalf, coarseStep);
    a = unique([-fliplr(pos(2:end)), pos]);
end
