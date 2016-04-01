%
% Name
%   mms_parse_time
%
% Purpose
%   Parse MMS file times.
%
% Calling Sequence
%   TVEC = mms_parse_time(TIMES);
%     TIMES is a string or cell array of strings containing the tstart
%     times from MMS file names. Their formatting is either 'yyyymmdd'
%     or 'yyyymmddHHMMSS'. These times are parsed and turned into a
%     Nx6 date vector, TVEC. If 'HHMMSS' are not provided, '000000' is
%     assumed.
%
%   [..., TT2000] = mms_parse_time(__);
%     Also return the CDF TT2000 times corresponding to TIMES.
%
% Parameters
%   TIMES:          in, required, type=char/cell
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-12-10      Written by Matthew Argall
%
function [tvec, tt2000] = mms_parse_time(times)
	% If input is a string, convert it to a cell
	if ischar(times) && isrow(times)
		times = { times };
	end

	% Find short times
	tf_short = regexp(times, '^[0-9]{4}[0-9]{2}[0-9]{2}$');
	tf_short = ~cellfun(@isempty, tf_short);
	
	% Find long times
	tf_long = regexp(times, '^[0-9]{4}[0-9]{2}[0-9]{2}[0-9]{2}[0-9]{2}[0-9]{2}$');
	tf_long = ~cellfun(@isempty, tf_long);

	% Make sure all were short or long
	assert( min( tf_short | tf_long ) == 1, 'Inputs must be formatted as "yyyymmdd" or "yyyymmddHHMMSS".' );
	
	% Convert to datevec
	%   - Datevec does not handle empty strings
	tvec = zeros( length(times), 9 );
	if sum(tf_short) > 0
		tvec(tf_short, 1:6) = datevec( times(tf_short), 'yyyymmdd' );
	end
	if sum(tf_long) > 0
		tvec(tf_long,  1:6) = datevec( times(tf_long),  'yyyymmddHHMMSS' );
	end

	% Convert to tt2000?
	if nargout() == 2
		tt2000 = spdfcomputett2000(tvec);
	end
	
	% Trim off the last 3 elements
	tvec = tvec(:, 1:6);
end