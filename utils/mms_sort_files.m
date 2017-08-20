%
% Name
%   mms_sort_files
%
% Purpose
%   Sort files by time.
%
% Calling Sequence
%   FSORTED = mms_parse_time(FILES);
%     Sort MMS file names given in FILES by their start time.
%
% Parameters
%   FILES:          in, required, type=char/cell
%
% Returns
%   FSORTED         out, required, type=char/cell
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2016-07-19      Written by Matthew Argall
%
function fsorted = mms_sort_files(files)

	% Parse the file names to get the start time
	[~, ~, ~, ~, fstart] = mms_dissect_filename(files);

	% Parse the times into tt2000 values
	[~, tt2000] = mms_parse_time(fstart);
	
	% Sort files by time
	[~, isort] = sort(tt2000);
	fsorted    = files(isort);
end