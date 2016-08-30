%
% Name
%   mms_latest_file
%
% Purpose
%   Search for the latest version of an MMS data file.
%
%   At the SDC, data files are initially saved in a dropbox location before
%   being mored to their final resting place in the directory structure. This
%   program provide a mechanism for searching in one folder or the other, or
%   both.
%
% Calling Sequence
%   FILENAME = mms_latest_file(DIR, SC, INSTR, MODE, LEVEL, TSTART);
%     Search for a file in directory DIR with spacecraft identifier SC,
%     instrument identifier INSTR, telemetry mode MODE, data level LEVEL,
%     and data start time TSTART. The file name with the highest version
%     number is returned.
%
%   [__] = mms_latest_file(..., 'ParamName', ParamValue);
%     Also supply any of the parameter name-value pairs listed below.
%
% Parameters
%   FOLDER:         in, required, type=char
%   SC:             in, required, type=char
%   INSTR:          in, required, type=char
%   MODE:           in, required, type=char
%   LEVEL:          in, required, type=char
%   TSTART:         in, required, type=char, default='%Y%M%d' or '%Y%M%d%H%m%S'
%   'OptDesc':      in, required, type=char, default=''
%                   Optional descriptor for the file name.
%   'Version':      in, required, type=char, default='*'
%                   File version number, formatted as 'X.Y.Z', where X, Y, and
%                     Z are integers.
%   'RootDir':      in, required, type=logical/char, default='*'
%                   If boolean and true, then `DIR` is taken to be the root of
%                     an SDC-like directory structure. If char, then "RootDir"
%                     is itself taken to be the root of an SDC-like directory
%                     structure. In the latter case, results from both ROOTDIR
%                     and DIR are returned.
%
% MATLAB release(s) MATLAB 8.2.0.701 (R2013b)
% Required Products None
%
% History:
%   2016-04-01      Written by Matthew Argall
%
function filename = mms_latest_file(folder, sc, instr, mode, level, tstart, varargin)

	% Search for files
	filename = mms_file_search(folder, sc, instr, mode, level, tstart, varargin{:});

	% Return the latest version
	if ~isempty(filename)
		filename = mms_latest_version(filename);
	end
end