function C = CaptureConditions(cfg, opts)
%CAPTURECONDITIONS  The acquisition conditions a capture must carry, resolved.
%
%   C = CaptureConditions(cfg) reads the free-text note that the operator wrote
%   into the .hws file (cfg is the second output of readHWS) and returns the
%   conditions the analysis needs. Anything the note does not carry is asked
%   for, so a capture with an incomplete note can still be processed instead of
%   failing or, worse, being processed against a silent assumption.
%
%   Required -- the numbers cannot be computed without them
%   -------------------------------------------------------
%     Home, Pos           the derating depth is |Pos - Home|
%     Scope_impedance     50 or 1e6; a factor 2 on MI and 4 on the intensities
%
%   Optional -- a documented default applies
%   ----------------------------------------
%     Pre-amp             "No" means the AH-2010 was out of the chain.
%                         ABSENT MEANS IT WAS IN, which is the convention the
%                         earlier sessions were recorded under.
%     Probe, TX_voltage, Pulse, PushCycles, PushElements
%                         provenance only; never enter the calculation.
%
%   What happens when a required field is missing
%   ---------------------------------------------
%   In an interactive session you are prompted for it. Under `matlab -batch`
%   there is nobody to ask, so it raises an error naming the field and the
%   option that supplies it -- a script never silently guesses.
%
%   Name-value options (any of these skip the note and the prompt)
%   -------------------------------------------------------------
%     Depth_cm             set the depth directly, ignoring Home/Pos
%     Home_mm              supply Home when the note records only Pos
%     ScopeImpedance_Ohm   50 or 1e6
%     WithPreamp           true / false
%     Interactive          force prompting on or off
%
%   Returned fields
%   ---------------
%     depth_cm, chain ("preamp"|"nopreamp"), scopeImpedance_Ohm,
%     probe, txVoltage_V, pulse, pushCycles, pushElements,
%     source  -- where each of the four resolved values came from:
%                "note", "argument", "prompt" or "default"
%
%   See also: readHWS, DepthFromNotes, DemoSingleCapture.

    arguments
        cfg struct
        opts.Depth_cm           double = []
        opts.Home_mm            double = []
        opts.ScopeImpedance_Ohm double = []
        opts.WithPreamp         = []
        opts.Interactive        = []
    end

    interactive = opts.Interactive;
    if isempty(interactive)
        interactive = ~batchStartupOptionUsed && usejava('jvm');
    end

    note = '';
    if isfield(cfg,'notes'), note = cfg.notes; end
    C = struct('source', struct());

    % ---------------------------------------------------------- provenance
    C.probe        = noteText(note, 'Probe');
    C.pulse        = noteText(note, 'Pulse');
    C.txVoltage_V  = noteNum(note,  'TX_voltage');
    C.pushCycles   = noteNum(note,  'PushCycles');
    C.pushElements = noteNum(note,  'PushElements');

    % ---------------------------------------------------- depth (required)
    if ~isempty(opts.Depth_cm)
        C.depth_cm = opts.Depth_cm;
        C.source.depth_cm = "argument";
    else
        pos  = noteCoord(note, 'Pos');
        home = opts.Home_mm(:)';
        if isempty(home), home = noteCoord(note, 'Home'); end
        if ~isempty(pos) && ~isempty(home)
            C.depth_cm = norm(pos - home)/10;
            if isempty(opts.Home_mm)
                C.source.depth_cm = "note";
            else
                C.source.depth_cm = "argument";   % Home came from the caller
            end
        else
            missing = 'Pos and Home';
            if ~isempty(pos),  missing = 'Home'; end
            if ~isempty(home), missing = 'Pos';  end
            C.depth_cm = askNumber(interactive, ...
                sprintf(['The capture note has no %s, so the derating depth is unknown.\n' ...
                         '  Depth = distance from the transducer face to the hydrophone tip.'], missing), ...
                'depth [cm]', 'Depth_cm', @(v) v >= 0 && v < 50);
            C.source.depth_cm = "prompt";
        end
    end

    % ------------------------------------------------ impedance (required)
    if ~isempty(opts.ScopeImpedance_Ohm)
        C.scopeImpedance_Ohm = opts.ScopeImpedance_Ohm;
        C.source.scopeImpedance_Ohm = "argument";
    else
        tk = regexp(note, 'Scope_impedance\s*=\s*([0-9.eE+]+)', 'tokens', 'once');
        if ~isempty(tk)
            C.scopeImpedance_Ohm = str2double(tk{1});
            C.source.scopeImpedance_Ohm = "note";
        else
            C.scopeImpedance_Ohm = askChoice(interactive, ...
                ['The capture note does not record the scope input impedance.', newline, ...
                 '  This matters: reading a pre-amplifier on 1 MOhm instead of 50 Ohm', newline, ...
                 '  overstates MI by 2x and the intensities by 4x.'], ...
                {'50 Ohm (terminated)', '1 MOhm (high-Z)'}, [50 1e6], 'ScopeImpedance_Ohm');
            C.source.scopeImpedance_Ohm = "prompt";
        end
    end

    % ------------------------------------------------- pre-amp (optional)
    if ~isempty(opts.WithPreamp)
        withPreamp = logical(opts.WithPreamp);
        C.source.chain = "argument";
    else
        tk = regexp(note, 'Pre-?amp\s*=\s*(\w+)', 'tokens', 'once');
        if ~isempty(tk)
            withPreamp = ~strcmpi(tk{1}, 'No');
            C.source.chain = "note";
        else
            withPreamp = true;              % the documented convention
            C.source.chain = "default";
        end
    end
    if withPreamp, C.chain = "preamp"; else, C.chain = "nopreamp"; end
end

% =====================================================================
%  note parsing
% =====================================================================
function v = noteCoord(note, key)
    v = [];
    if isempty(note), return, end
    tk = regexp(note, [key '\s*=\s*\(\s*([-\d.]+)\s*,\s*([-\d.]+)\s*,\s*([-\d.]+)'], ...
                'tokens', 'once');
    if ~isempty(tk), v = str2double(tk)'; v = v(:)'; end
end

function s = noteText(note, key)
    s = "";
    if isempty(note), return, end
    tk = regexp(note, [key '\s*=\s*([^\r\n]+)'], 'tokens', 'once');
    if ~isempty(tk), s = string(strtrim(tk{1})); end
end

function v = noteNum(note, key)
    v = NaN;
    if isempty(note), return, end
    tk = regexp(note, [key '\s*=\s*([-\d.]+)'], 'tokens', 'once');
    if ~isempty(tk), v = str2double(tk{1}); end
end

% =====================================================================
%  prompting
% =====================================================================
function v = askChoice(interactive, why, labels, values, optName)
%ASKCHOICE  Offer a short menu; error instead if nobody can answer.
    requireInteractive(interactive, why, optName);
    fprintf('\n%s\n\n', why);
    for k = 1:numel(labels)
        fprintf('    [%d] %s\n', k, labels{k});
    end
    while true
        r = input(sprintf('  choose 1-%d: ', numel(labels)));
        if isnumeric(r) && isscalar(r) && any(r == 1:numel(labels))
            v = values(r);
            fprintf('  -> using %s\n\n', labels{r});
            return
        end
        fprintf('  please enter a number between 1 and %d.\n', numel(labels));
    end
end

function v = askNumber(interactive, why, label, optName, isValid)
%ASKNUMBER  Ask for one number; error instead if nobody can answer.
    requireInteractive(interactive, why, optName);
    fprintf('\n%s\n\n', why);
    while true
        v = input(sprintf('  %s: ', label));
        if isnumeric(v) && isscalar(v) && isfinite(v) && isValid(v)
            fprintf('  -> using %g\n\n', v);
            return
        end
        fprintf('  that is not a usable value; try again.\n');
    end
end

function requireInteractive(interactive, why, optName)
    if interactive, return, end
    error('CaptureConditions:missingField', ...
        ['%s\n\nThis session cannot prompt (running under -batch or without a JVM).\n' ...
         'Pass %s= explicitly, or add the field to the capture note.'], why, optName);
end
