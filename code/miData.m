function p = miData(varargin)
%MIDATA  Path into the bundled example captures (data/).
%
%   miData()                          -> ...\MI-measurements\data
%   miData('S5-1','2026-08-13_sessionA_preamp_50ohm','push_79el_sweep')
%
%   Set the environment variable MI_DATA_ROOT to point somewhere else (e.g. a
%   network share holding the full, unabridged measurement set); the bundled
%   data/ folder is used when it is unset.

    root = getenv('MI_DATA_ROOT');
    if isempty(root), root = miRoot('data'); end
    p = fullfile(root, varargin{:});
end
