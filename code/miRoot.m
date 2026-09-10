function p = miRoot(varargin)
%MIROOT  Absolute path to the repository root, or to something inside it.
%
%   miRoot()                 -> ...\MI-measurements
%   miRoot('docs','report')  -> ...\MI-measurements\docs\report
%
%   Every path in this repository is built from here, so nothing has to be
%   edited after a clone.

    p = fileparts(fileparts(mfilename('fullpath')));   % code\miRoot.m -> root
    if ~isempty(varargin), p = fullfile(p, varargin{:}); end
end
