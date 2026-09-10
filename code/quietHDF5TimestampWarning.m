function c = quietHDF5TimestampWarning()
%QUIETHDF5TIMESTAMPWARNING  Silence one specific, unavoidable HDF5 warning.
%
%   c = quietHDF5TimestampWarning() turns off the warning below and restores
%   the caller's warning state when c goes out of scope. Call it immediately
%   before an h5info on a .hws file:
%
%       quiet = quietHDF5TimestampWarning();  %#ok<NASGU>
%       info  = h5info(file);
%
%   What it silences
%   ----------------
%       There is no MATLAB integer type that corresponds to an HDF5 integer
%       type with 16 bytes. The data will be treated as an uninterpreted
%       uint8 array.
%       (MATLAB:imagesci:hdf5dataset:integerButNoSuchMatlabClass)
%
%   Why it is safe
%   --------------
%   Every .hws capture carries three NI 128-bit timestamp attributes:
%       /wfm_group0/id           timestamp
%       /wfm_group0/axes/axis0   ref_time
%       /wfm_group0/axes/axis1   ref_time
%   MATLAB has no 128-bit integer type, so h5info warns whenever it
%   ENUMERATES them -- once per attribute, per file. Nothing in this
%   repository reads them: the samples, the scale coefficients, the time axis
%   and the whole scope configuration are unaffected, and every number the
%   analysis produces is identical with or without the warning.
%
%   It cannot be avoided by requesting a different data type. The 16-byte
%   integer is a property of the file, written by the NI Soft Front Panel, and
%   h5info reports on it before any type could be chosen.
%
%   Suppression is therefore by IDENTIFIER and scoped to the single call.
%   Do NOT reach for warning('off','all') instead: readHWS raises two warnings
%   that matter -- RIS (equivalent-time) sampling, which is invalid for a
%   single-shot push, and a record shorter than the pulse, which silently
%   truncates PII and both intensities.
%
%   See also: readHWS.

    id  = 'MATLAB:imagesci:hdf5dataset:integerButNoSuchMatlabClass';
    old = warning('off', id);
    c   = onCleanup(@() warning(old));
end
