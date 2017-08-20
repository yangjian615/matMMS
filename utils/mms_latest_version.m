%
% Name
%   mms_latest_version
%
% Purpose
%   Find the newest version of a given set of files. File names are assumed to be
%   identical, with the exception of their version numbers.
%
% Calling Sequence
%   NEWEST_FILE = mms_latest_version(FILES);
%     Filter CDF file names given in FILES for the newest version and
%     return them in NEWEST_FILE.
%
%   [..., VERSION] = mms_latest_version(__);
%     Also return the version of the newest file, formatted as 'X.Y.Z', where X,
%     Y, and Z are integers.
%
% Parameters
%   FILES:          in, required, type=char/cell
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
function [newest_file, version] = mms_latest_version(files)

	% Single file: Already done
	if ischar(files) && isrow(files)
		newest_file = files;
		return
		
	% Cell array of file names
	elseif iscell(files)
		newest_file = files;
	
	% Bad datatype
	else
		error('MMS:badDatatype', 'FILES must be a single or cell array of char row vectors.')
	end
	
	% Number of files given
	nFiles = length(files);
	
	% File names end in vX.Y.Z.cdf, where X, Y, and Z are integers
	%   - Extract X, Y, and Z
	v = regexp(newest_file, 'v([0-9]+)\.([0-9]+)\.([0-9]+)\.cdf$', 'tokens');
	v = vertcat( v{:} );
	v = reshape( cellfun(@str2num, [ v{:} ]), 3, nFiles );

	% Step through the X-, Y-, and Z-verison numbers
	for ii = 1 : 3
		% Find highest version
		idx = find( v(ii,:) == max( v(ii,:) ) );
		
		% Filter out
		v = v(:,idx);
		newest_file = newest_file(:,idx);
	end

	% Change from cell to string
	%   - If there is more than one file with the same version, this
	%     selects the first.
	newest_file = newest_file{1};

	% Return the latest version
	if nargout() == 2
		version = sprintf('%i.%i.%i', v(:,1));
	end
end