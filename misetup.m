function root = misetup()
%MISETUP  Put this repository on the MATLAB path. Run once per session.
%
%   misetup            adds code/, code/analysis/, code/acquisition/ to the path
%   root = misetup()   also returns the repository root folder
%
%   Everything in this repository resolves its own paths relative to that root
%   (see miRoot, miData, miCalibration), so the repository can be cloned
%   anywhere and needs no editing.
%
%   Requires MATLAB (developed on R2025b). Base MATLAB only -- `hann` comes
%   from Signal Processing Toolbox but is only used by the analysis scripts,
%   not by readHWS/SafetyIndices.

    root = fileparts(mfilename('fullpath'));
    addpath(fullfile(root,'code'), ...
            fullfile(root,'code','analysis'), ...
            fullfile(root,'code','acquisition'), ...
            fullfile(root,'demo'));
    fprintf('MI-measurements on the path: %s\n', root);
    if nargout == 0, clear root; end
end
