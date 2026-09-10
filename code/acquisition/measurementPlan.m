function plan = measurementPlan(beams, options)
%MEASUREMENTPLAN Full hydrophone measurement protocol for several beam types.
%   plan = measurementPlan(beams) returns, for each beam definition, the
%   complete "what to measure with which parameters" recipe:
%       1. SEARCH  - a coarse grid to LOCATE the pressure peak (measure
%                    around the observed maximum, not the geometric focus).
%       2. FINE    - a function handle that, given the found peak [x y z],
%                    produces the fine measurement (anchor / validation cross
%                    / full scan) centred on that peak.
%       3. SWEEPS  - the 1-D parameter sweeps (transmit voltage, pulse
%                    length) done at the single peak point, not re-scanned
%                    over space.
%
%   'beams' is a struct array; each element may define:
%       .name        - label (string)                         [required]
%       .f0          - centre frequency [Hz]                  [required]
%       .focalDepth  - nominal focal / search-centre depth [mm][required]
%       .type        - "focused" | "widebeam" | "diverging"   (default focused)
%       .Fnum        - f-number                               (or .aperture)
%       .aperture    - active aperture [mm] (used if Fnum absent)
%       .tier        - fine-measurement tier once the peak is found:
%                        "anchor"   -> 1 pt   (trust the simulation shape)
%                        "validate" -> axial+azimuth cross (~30 pts)
%                        "full"     -> 3-line scan (~60 pts)   (default validate)
%       .voltages    - transmit-voltage sweep levels (at the peak point)
%       .cycles      - pulse-length sweep in cycles (at the peak point)
%
%   Name-value options (applied to every beam): c, gridStep, elevationFnum,
%   plus searchStepFrac / searchRangeFWHM / searchAxialSpan to size the
%   coarse search (see below). Anything else you want to override per fine
%   scan can be set later on the returned function handle.
%
%   Output 'plan' is a struct array with fields: .name .params .scales
%   .search (struct from scanPlan) .fine (function handle) .tier .sweeps.
%   A human-readable protocol is printed for each beam.
%
%   Typical use:
%       beams(1) = struct('name',"focused",  'f0',3e6,'focalDepth',40, ...
%                         'type',"focused",  'Fnum',2, 'tier',"validate", ...
%                         'voltages',10:10:50,'cycles',[2 4 8]);
%       beams(2) = struct('name',"diverging",'f0',3e6,'focalDepth',20, ...
%                         'type',"diverging",'Fnum',1, 'tier',"anchor", ...
%                         'voltages',10:10:50,'cycles',[]);
%       P = measurementPlan(beams);
%       % ... run P(1).search, find the peak, then:
%       fine1 = P(1).fine([xpk ypk zpk]);   % fine cross on the observed peak

    arguments
        beams (1,:) struct
        options.c              (1,1) double = 1480
        options.gridStep       (1,1) double = 0.01
        options.elevationFnum  (1,1) double = NaN
        options.searchRangeFWHM (1,1) double = 3     % lateral search half-range
        options.searchStepFrac  (1,1) double = 0.75  % coarse lateral step / FWHM
        options.searchAxialSpan (1,1) double = 1.5    % axial search half-range / DOF
        options.searchAxialStep (1,1) double = 2      % coarse axial step [mm]
    end

    plan = struct('name', {}, 'params', {}, 'scales', {}, 'search', {}, ...
                  'fine', {}, 'fineArgs', {}, 'tier', {}, 'sweeps', {});

    for i = 1:numel(beams)
        b = beams(i);
        b = fillDefaults(b);

        % common scanPlan args for this beam
        common = {'c', options.c, 'gridStep', options.gridStep};
        if ~isnan(options.elevationFnum)
            common = [common, {'elevationFnum', options.elevationFnum}]; %#ok<AGROW>
        end
        if ~isempty(b.Fnum) && ~isnan(b.Fnum)
            common = [common, {'Fnum', b.Fnum}];       %#ok<AGROW>
        elseif ~isempty(b.aperture) && ~isnan(b.aperture)
            common = [common, {'aperture', b.aperture}]; %#ok<AGROW>
        end

        % ---- 1) SEARCH: coarse grid to find the peak -------------
        %   Diverging/widebeam beams have no focal peak near focalDepth, so
        %   search a wide axial span from the near field down past it.
        if b.type == "diverging" || b.type == "widebeam"
            axSpan   = max(options.searchAxialSpan * 20, b.focalDepth); % mm each side
            searchCtr = [0 0 max(b.focalDepth/2, 10)];  % start shallow
            axStep    = options.searchAxialStep;
        else
            searchCtr = [0 0 b.focalDepth];
            axSpan    = options.searchAxialSpan;         % in DOF units below
            axStep    = options.searchAxialStep;
        end

        fprintf('\n===== BEAM %d: %s (%s) =====\n', i, b.name, b.type);
        fprintf('[1] SEARCH for the peak (coarse, then recentre):\n');
        if b.type == "diverging" || b.type == "widebeam"
            search = scanPlan(b.f0, b.focalDepth, common{:}, ...
                'lines', ["axial","azimuth"], 'center', searchCtr, ...
                'axRangeDOF', 0, 'axProxExtra', 0, ...
                'axFineStep', axStep, 'axCoarseStep', axStep, ...
                'latRangeFWHM', options.searchRangeFWHM, ...
                'latFineFrac', options.searchStepFrac, ...
                'latCoarseFrac', options.searchStepFrac);
            % widen axial manually for divergent case (flat field)
            search = widenAxial(search, searchCtr, axSpan, axStep, options.c, b);
            fprintf(['  (wide/diverging: axial search widened to z=[%.0f, %.0f] mm, ' ...
                'step %.1f mm; peak may lie in the near field)\n'], ...
                search.z(1), search.z(end), axStep);
        else
            search = scanPlan(b.f0, b.focalDepth, common{:}, ...
                'lines', ["axial","azimuth"], 'center', searchCtr, ...
                'axRangeDOF', axSpan, 'axProxExtra', 0.3, ...
                'axFineStep', axStep, 'axCoarseStep', axStep, ...
                'latRangeFWHM', options.searchRangeFWHM, ...
                'latFineFrac', options.searchStepFrac, ...
                'latCoarseFrac', options.searchStepFrac);
        end

        % ---- 2) FINE: recipe centred on the FOUND peak -----------
        tier = b.tier;
        switch tier
            case "anchor";   fineArgs = {'anchorOnly', true};
            case "validate"; fineArgs = {'lines', ["axial","azimuth"], ...
                                         'latFineFrac', 0.5, 'axFineStep', 1};
            otherwise;       fineArgs = {};   % full
        end
        % Freeze the concrete scanPlan arguments for this beam (everything
        % except 'center') and build the handle from an explicit copy via a
        % helper, so each beam's .fine is self-contained and inspectable via
        % .fineArgs rather than capturing loop workspace variables.
        baseArgs = [{b.f0, b.focalDepth}, common, fineArgs];
        fcn = makeFineFcn(baseArgs);
        fprintf('[2] FINE measurement (tier=%s) -> call plan(%d).fine(peakXYZ)\n', tier, i);

        % ---- 3) SWEEPS at the single peak point ------------------
        sweeps = struct('voltages', b.voltages, 'cycles', b.cycles, ...
                        'location', 'observed peak (1 point)');
        fprintf('[3] SWEEPS at the peak point only:\n');
        if ~isempty(b.voltages)
            fprintf('      transmit voltage: [%s]  -> pr-vs-V curve (nonlinear scaling)\n', ...
                num2str(b.voltages));
        end
        if ~isempty(b.cycles)
            fprintf('      pulse length:     [%s] cycles -> PII / Ispp scaling\n', ...
                num2str(b.cycles));
        end

        % Assign fields individually: storing a cell (baseArgs) or a handle
        % via struct(...) would trip the struct() cell-expansion trap.
        plan(i).name    = b.name;
        plan(i).params  = struct('f0', b.f0, 'focalDepth', b.focalDepth, ...
                                 'Fnum', b.Fnum, 'aperture', b.aperture, 'type', b.type);
        plan(i).scales  = struct('lambda_mm', search.lambda_mm, ...
                                 'lat_FWHM_mm', search.lat_FWHM_mm, ...
                                 'axial_DOF_mm', search.axial_DOF_mm);
        plan(i).search  = search;
        plan(i).fine    = fcn;
        plan(i).fineArgs = baseArgs;   % concrete frozen args (inspectable)
        plan(i).tier    = tier;
        plan(i).sweeps  = sweeps;
    end

    fprintf('\nWorkflow per beam: run .search -> locate peak -> .fine(peak) -> .sweeps at peak.\n');
    fprintf('Voltage/pulse-length are NOT re-scanned in space; they scale from the peak point.\n');
end

% ===================================================================
function h = makeFineFcn(baseArgs)
%MAKEFINEFCN Build the fine-scan handle from an explicit argument copy.
%   Creating the handle inside its own function guarantees it binds this
%   beam's frozen 'baseArgs' (passed by value) and nothing from the caller's
%   loop workspace.
    h = @(peakXYZ) scanPlan(baseArgs{:}, 'center', peakXYZ);
end

% ===================================================================
function b = fillDefaults(b)
    if ~isfield(b,'type') || isempty(b.type);   b.type = "focused"; end
    if ~isfield(b,'tier') || isempty(b.tier);   b.tier = "validate"; end
    if ~isfield(b,'Fnum');       b.Fnum = NaN;     end
    if ~isfield(b,'aperture');   b.aperture = NaN; end
    if ~isfield(b,'voltages');   b.voltages = [];  end
    if ~isfield(b,'cycles');     b.cycles = [];    end
    b.type = string(b.type);
    b.tier = string(b.tier);
    b.name = string(b.name);
end

% ===================================================================
function s = widenAxial(s, ctr, axSpanMM, axStep, c, b) %#ok<INUSD>
%WIDENAXIAL Replace the axial line with a wide, uniform near-field sweep for
%   diverging/wide beams whose peak can lie anywhere from ~near field to the
%   nominal depth. Rebuilds points and tof.
    zline = (max(1, ctr(3)-axSpanMM)):axStep:(ctr(3)+axSpanMM);
    axial = [ctr(1)*ones(numel(zline),1), ctr(2)*ones(numel(zline),1), zline(:)];
    s.axial = axial;
    s.z = zline(:);
    s.points = unique([axial; s.azimuth], 'rows');
    s.tof_s  = vecnorm(s.points, 2, 2) * 1e-3 / c;
    s.nPoints = size(s.points, 1);
end
