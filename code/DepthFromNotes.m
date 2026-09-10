function [depth_cm, pos, home] = DepthFromNotes(cfg, Home_mm)
%DEPTHFROMNOTES  Measurement depth [cm] from the note stored in the .hws file.
%
%   The capture note written at acquisition time carries the positioner
%   coordinates, e.g.
%       Home = (55.8,105.84,-17.3) mm      <- transducer centre
%       Pos  = (46.1,70.84,-17.3) mm       <- hydrophone tip
%   and the depth used for derating is |Pos - Home|.
%
%   [depth_cm, pos, home] = DepthFromNotes(cfg) reads both from the note.
%   DepthFromNotes(cfg, Home_mm) overrides Home (older captures wrote only Pos).
%
%   cfg is the second output of readHWS.

    pos  = coords(cfg, 'Pos');
    if nargin > 1 && ~isempty(Home_mm)
        home = Home_mm(:)';
    else
        home = coords(cfg, 'Home');
    end
    if isempty(pos) || isempty(home)
        error('DepthFromNotes:missing', ...
              'Capture note has no Pos/Home coordinates; pass Home_mm or set the depth by hand.');
    end
    depth_cm = norm(pos - home)/10;
end

function v = coords(cfg, key)
    v = [];
    if ~isfield(cfg,'notes') || isempty(cfg.notes), return, end
    tk = regexp(cfg.notes, [key '\s*=\s*\(\s*([-\d.]+)\s*,\s*([-\d.]+)\s*,\s*([-\d.]+)'], ...
                'tokens','once');
    if ~isempty(tk), v = str2double(tk)'; v = v(:)'; end
end
