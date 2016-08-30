%
% Name
%   mms_latest_zversion
%
% Purpose
%   Return the next newest Z-version of a given file. The z-version returned is one
%   larger than the latest existing file.
%
% Calling Sequence
%   VZ = mms_latest_zversion(FOLDER, FILENAME);
%     Search for files with names FILES in directory FOLDER.
%     File z-versions are examined and the next highest available
%     z-version is returned in VZ.
%
% Parameters
%   FILES:          in, required, type=char/cell
%   'RootDir':      in, optional, type=logical/char, default='*'
%                   If boolean and true, then `FOLDER` is taken to be the root of
%                     an SDC-like directory structure. If char, then "RootDir"
%                     is itself taken to be the root of an SDC-like directory
%                     structure. In the latter case, results from both ROOTDIR
%                     and FOLDER are returned.
%
% Output
%   NEWEST_FILE:    out, required, type=char
%   VERSION:        out, optional, type=char
%
% MATLAB release(s) MATLAB 8.2.0.701 (R2013b)
% Required Products None
%
% History:
%   2016-04-01      Written by Matthew Argall
%
function vz = mms_latest_zversion(folder, filename, varargin)

	% Only one file allowed.
	if ~ischar(filename) && ~isrow(filename)
		error('MMS:badDatatype', 'FILENAME must be a row vector of class char.')
	end
	
	% Remove the directory on FILENAME
	[~, fname, ext] = fileparts(filename);
	fname           = [fname ext];
	
	% Dissect the file name and version number
	[sc, instr, mode, level, tstart, version, optdesc] = mms_dissect_filename(filename);
	[vx, vy, vz] = mms_parse_version(version);

	% Replace the z-version number with wild card and search
	vout     = [vx '.' vy '.*'];
	fXYmatch = mms_file_search(folder, sc, instr, mode, level, tstart, ...
	                           'OptDesc', optdesc, ...
	                           'Version', vout,    ...
	                           varargin{:} );

	% If no match, the z-version is zero
	if isempty(fXYmatch)
		vz = 0;
	
	% Match
	else
		% Dissect the version numbers
		[~, ~, ~, ~, ~, version] = mms_dissect_filename(fXYmatch);
		[vx, vy, vz] = mms_parse_version(version);

		% Find maximum and add one for newest
		if iscell(fXYmatch)
			vz = cellfun(@str2num, vz);
		else
			vz = str2num(vz);
		end
		vz = max(vz) + 1;
	end
end