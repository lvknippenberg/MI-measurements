function [wfm, cfg] = readHWS(file, expectedPulseLen)
%READHWS Read an NI Hierarchical Waveform Storage (.hws) scope capture.
%   [wfm, cfg] = readHWS(file) returns the waveform(s) and the scope
%   configuration that was saved inside the file.
%
%   [wfm, cfg] = readHWS(file, expectedPulseLen) additionally warns if the
%   captured record is shorter than expectedPulseLen seconds (the pulse
%   would be truncated). A warning is always issued when the file was
%   acquired with RIS (equivalent-time) sampling, which is only valid for
%   repetitive signals and NOT for single-shot push/transient waveforms.
%
%   wfm is a struct array, one element per acquired trace/channel:
%       .name  - channel name (e.g. 'Chan 0')
%       .t     - time vector in seconds (column vector)
%       .y     - scaled samples in the y-axis units, i.e. volts (column vector)
%       .v     - alias of .y (kept for backward compatibility)
%       .units - y-axis units string (e.g. 'V')
%       .dt    - sample interval (s)
%       .fs    - sample rate (Hz)
%       .t0    - time of first sample (s)
%
%   cfg is a struct with the acquisition settings the scope stored:
%       .sample_rate        - Hz
%       .dt                 - sample interval (s)
%       .record_length      - samples per record
%       .number_of_records  - segments acquired (1 = single record)
%       .acquisition_type   - raw enum (0 = normal)
%       .reference_pos_pct  - trigger reference position (% of record)
%       .channels(k)        - per-channel struct:
%             .index, .name, .enabled,
%             .bandwidth_Hz, .probe_attenuation,
%             .input_impedance_ohm  (decoded), .input_impedance_raw,
%             .coupling             (decoded 'AC'/'DC'/'GND'), .coupling_raw,
%             .vertical_range_V, .vertical_offset_V
%       .writer             - software that wrote the file
%
%   .hws files are HDF5 underneath. Samples are stored as raw ADC codes and
%   converted to physical units via the axis scale_coef polynomial:
%       value = c0 + c1*raw   (+ higher-order terms if present)
%
%   Example:
%       [w, cfg] = readHWS('FunctionGeneratorSine.hws');
%       plot(w(1).t*1e6, w(1).y); xlabel('\mus'); ylabel(w(1).units);
%       fprintf('%.0f MS/s, %d samples, %d records, %g MHz BW\n', ...
%           cfg.sample_rate/1e6, cfg.record_length, ...
%           cfg.number_of_records, cfg.channels(1).bandwidth_Hz/1e6);

    % ---------------------------------------------------------------
    %  Waveform(s)
    % ---------------------------------------------------------------
    quiet = quietHDF5TimestampWarning();  %#ok<NASGU>  see that function for why
    info  = h5info(file, '/wfm_group0/traces');
    nTr  = numel(info.Groups);
    wfm  = struct('name', {}, 't', {}, 'y', {}, 'v', {}, 'units', {}, ...
                  'dt', {}, 'fs', {}, 't0', {});

    for k = 1:nTr
        tracePath = info.Groups(k).Name;                       % .../trace0
        name  = h5readatt(file, tracePath, 'name');
        xIdx  = h5readatt(file, tracePath, 'x_axis');          % axis index
        yIdx  = h5readatt(file, tracePath, 'y_axis');
        xAxis = sprintf('/wfm_group0/axes/axis%d', xIdx);
        yAxis = sprintf('/wfm_group0/axes/axis%d', yIdx);

        % --- y data: raw samples -> physical units ---
        raw = double(readTraceData(file, tracePath, yAxis));
        sc  = h5read(file, [yAxis '/scale_coef']);             % [c0 c1 ...]
        y   = zeros(size(raw));
        for p = 1:numel(sc)
            y = y + sc(p) * raw.^(p-1);                        % polynomial
        end
        units = h5readatt(file, yAxis, 'units');

        % --- x (time) axis: implicit start + increment ---
        t0 = h5readatt(file, xAxis, 'start');
        dt = h5readatt(file, xAxis, 'increment');
        t  = t0 + (0:numel(y)-1)' * dt;

        wfm(k) = struct('name', char(name), 't', t(:), ...
                        'y', y(:), 'v', y(:), 'units', char(units), ...
                        'dt', dt, 'fs', 1/dt, 't0', t0);
    end

    % ---------------------------------------------------------------
    %  Configuration (best-effort; empty if the file has no cfg group)
    % ---------------------------------------------------------------
    if nargout > 1 || nargin > 1
        cfg = readConfig(file);

        % --- sanity warnings -------------------------------------
        if isequal(cfg.is_RIS, true)
            warning('readHWS:RIS', ['File acquired with RIS (equivalent-time) ' ...
                'sampling at %.0f MS/s. Valid only for repetitive signals, ' ...
                'NOT single-shot push/transient waveforms.'], cfg.sample_rate/1e6);
        end
        if nargin > 1 && ~isempty(expectedPulseLen) && ~isempty(cfg.record_duration) ...
                && cfg.record_duration < expectedPulseLen
            warning('readHWS:ShortRecord', ['Record duration %.3g s is shorter ' ...
                'than the expected pulse length %.3g s - the pulse will be ' ...
                'truncated. Increase record length or lower the sample rate.'], ...
                cfg.record_duration, expectedPulseLen);
        end
    end
end

% ===================================================================
function raw = readTraceData(file, tracePath, yAxis)
%READTRACEDATA Locate the raw sample vector for a trace.
%   Samples usually live in the y-axis' data_vector reference, which points
%   at /wfm_group0/vectors/vectorN/data. Falls back to vector0 for the
%   common single-channel case.
    try
        vecRef = h5readatt(file, [yAxis '/data_vector'], 'vector');
        raw = h5read(file, sprintf('/wfm_group0/vectors/vector%d/data', vecRef));
    catch
        raw = h5read(file, '/wfm_group0/vectors/vector0/data');
    end
end

% ===================================================================
function cfg = readConfig(file)
%READCONFIG Pull acquisition settings out of the /cfg_scope0 group.
    cfg = struct('sample_rate', [], 'dt', [], 'record_length', [], ...
                 'record_duration', [], 'number_of_records', [], ...
                 'acquisition_type', [], 'sample_mode', [], ...
                 'reference_pos_pct', [], 'is_RIS', [], ...
                 'digitizer', '', 'notes', '', ...
                 'channels', struct([]), 'writer', '');

    try
        cfg.writer = char(h5readatt(file, '/wfm_group0', 'writer'));
    catch
    end

    % --- digitizer model ---
    try
        cfg.digitizer = char(h5read(file, '/cfg_scope0/data/glb/public/names'));
    catch
    end

    % --- user notes (e.g. push location) ---
    for np = {'/wfm_group0/id', '/cfg_scope0/data/file/public/id'}
        for na = {'note', 'notes'}
            try
                s = char(h5readatt(file, np{1}, na{1}));
                if ~isempty(strtrim(s)); cfg.notes = s; end
            catch
            end
        end
        if ~isempty(cfg.notes); break; end
    end

    % --- horizontal / timebase ---
    hor = '/cfg_scope0/data/hor/public0';
    cfg.sample_rate       = attn(file, hor, 'sample_rate');
    cfg.record_length     = attn(file, hor, 'record_length');
    cfg.number_of_records = attn(file, hor, 'number_of_records');
    cfg.acquisition_type  = attn(file, hor, 'acquisition_type');
    cfg.sample_mode       = attn(file, hor, 'sample_mode');
    cfg.reference_pos_pct = attn(file, hor, 'reference_position');
    if ~isempty(cfg.sample_rate) && cfg.sample_rate > 0
        cfg.dt = 1 / cfg.sample_rate;
        if ~isempty(cfg.record_length)
            cfg.record_duration = cfg.record_length / cfg.sample_rate;
        end
    end

    % --- RIS (equivalent-time) sampling detection ---
    %   RIS_on is a per-timebase flag vector; time-div_index selects the
    %   active timebase. If that entry is set, the scope reconstructed the
    %   record from many triggers (repetitive-signal only).
    try
        risVec = double(h5read(file, '/cfg_scope0/data/hor/private/RIS_on'));
        tdi    = attn(file, '/cfg_scope0/data/hor/private', 'time-div_index');
        if ~isempty(tdi) && tdi+1 >= 1 && tdi+1 <= numel(risVec)
            cfg.is_RIS = risVec(tdi+1) == 1;
        end
    catch
    end

    % --- per-channel vertical settings ---
    nCh = attn(file, '/cfg_scope0/data/chan/private', 'nChannels');
    if isempty(nCh); nCh = 0; end
    ch = struct('index', {}, 'name', {}, 'enabled', {}, ...
                'bandwidth_Hz', {}, 'probe_attenuation', {}, ...
                'input_impedance_ohm', {}, 'input_impedance_raw', {}, ...
                'coupling', {}, 'coupling_raw', {}, ...
                'vertical_range_V', {}, 'vertical_offset_V', {});
    for c = 0:nCh-1
        p = sprintf('/cfg_scope0/data/chan/publicUI%d', c);
        zraw = attn(file, p, 'input_impedance');
        craw = attn(file, p, 'vertical_coupling');
        ch(end+1) = struct( ...                                    %#ok<AGROW>
            'index',               c, ...
            'name',                sprintf('Chan %d', c), ...
            'enabled',             logical(nz(attn(file, p, 'channel_enabled'))), ...
            'bandwidth_Hz',        attn(file, p, 'bandwidth'), ...
            'probe_attenuation',   attn(file, p, 'probe_attenuation'), ...
            'input_impedance_ohm', decodeImpedance(zraw), ...
            'input_impedance_raw', zraw, ...
            'coupling',            decodeCoupling(craw), ...
            'coupling_raw',        craw, ...
            'vertical_range_V',    attn(file, p, 'vertical_range'), ...
            'vertical_offset_V',   attn(file, p, 'vertical_offset'));
    end
    cfg.channels = ch;
end

% ---- small helpers -------------------------------------------------
function v = attn(file, node, name)
%ATTN  Read a numeric attribute, returning [] if it is absent.
    try
        v = double(h5readatt(file, node, name));
    catch
        v = [];
    end
end

function x = nz(v)
    if isempty(v); x = 0; else; x = v; end
end

function s = decodeImpedance(raw)
%DECODEIMPEDANCE  Scope-SFP input_impedance enum -> ohms.
%   Known mapping: 0 -> 1e6 (1 MOhm), 1 -> 50. Returns NaN if unknown.
    if isempty(raw);      s = NaN;
    elseif raw == 0;      s = 1e6;
    elseif raw == 1;      s = 50;
    else;                 s = NaN;
    end
end

function s = decodeCoupling(raw)
%DECODECOUPLING  Scope-SFP vertical_coupling enum -> string.
%   Known mapping: 0 -> AC, 1 -> DC, 2 -> GND.
    switch nz(raw)
        case 0; s = 'AC';
        case 1; s = 'DC';
        case 2; s = 'GND';
        otherwise; s = 'unknown';
    end
end
