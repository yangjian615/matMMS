%
% Name
%   mms_comp_version
%
% Purpose
%   Compare the versions of two files in MMS.
%
% Calling Sequence
%   RESULT = mms_comp_version(FILE1, FILE2)
%       Compares the version of FILE1 to that of FILE2. The RESULT is -1 if
%       FILE1 is older, 0 if FILE1 is the same version, and 1 if FILE1 is
%       newer.
%
%   RESULT = mms_comp_version(__, REGEX)
%       Specify the regular expression REGEX used to extract version
%       numbers. Used with the 'tokens' flag set in regexp().
%
% Parameters
%   FILE1           in, required, type = char
%   FILE2           in, required, type = char
%   REGEX           in, optional, type = char, default = '([0-9]+)\.([0-9]+)\.([0-9]+)'
%
% Returns:
%   RESULT          out, required, type = integer
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-02      Written by Matthew Argall
%
function [result] = mms_comp_version(file1, file2, regex)

	% String for extracting version elements
	if nargin < 3
		regex = '([0-9]+)\.([0-9]+)\.([0-9]+)';
	end
	
	% Extract the version from each file
	vfile1 = regexp(file1, regex, 'tokens');
	vfile2 = regexp(file2, regex, 'tokens');
	
	% Get rid of the nested cells
	vfile1 = vfile1{1};
	vfile2 = vfile2{1};
	
	% Compare	
	ii     = 1;
	result = 0;
	while result == 0 && ii <= length(vfile1)
		
		% Newer version
		if vfile1{ii} >= vfile2{ii}
			result = 1;
			
		% Older version
		elseif vfile1{ii} <= vfile2{ii}
			result = -1;
		end
		
		% Move to the next part
		ii = ii + 1;
	end
end