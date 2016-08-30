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
%     Nx6 cell array of strings, TVEC. If 'HHMMSS' are not provided,
%    '000000' is assumed.
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
%   2016-04-02      Return cell array of strings.
%
function [tvec, tt2000] = mms_parse_time(times)

%------------------------------------%
% Check Inputs                       %
%------------------------------------%
	
	% If input is a string, convert it to a cell
	if ischar(times) && isrow(times)
		times = { times };
	end
	
	% Number of times given
	nPts = length(times);

	% Find short and long times
	tf_short = ~cellfun(@isempty, regexp(times, '^[0-9]{4}[0-9]{2}[0-9]{2}$') );
	tf_long  = ~cellfun(@isempty, regexp(times, '^[0-9]{4}[0-9]{2}[0-9]{2}[0-9]{2}[0-9]{2}[0-9]{2}$') );

	% Make sure all were short or long
	assert( min( tf_short | tf_long ) == 1, 'Inputs must be formatted as "yyyymmdd" or "yyyymmddHHMMSS".' );

%------------------------------------%
% Convert to Date-Time Vectors       %
%------------------------------------%
	
	% Allocate memory
	tvec = cell( nPts, 6 );
	
	% Short format
	if sum(tf_short) > 0
		% Parse string
		tstr = regexp( times(tf_short), '^([0-9]{4})([0-9]{2})([0-9]{2})$', 'tokens' );
		tstr = vertcat( tstr{:} );
		tstr = vertcat( tstr{:} );

		% Separate time components
		tvec(tf_short, 1) = tstr(:,1);      % Year
		tvec(tf_short, 2) = tstr(:,2);      % Month
		tvec(tf_short, 3) = tstr(:,3);      % Day
		tvec(tf_short, 4) = { '00' };       % Hour
		tvec(tf_short, 5) = { '00' };       % Minute
		tvec(tf_short, 6) = { '00' };       % Second
	end
	
	% Long format
	if sum(tf_long) > 0
		% Parse string
		tstr = regexp( times(tf_long), '^([0-9]{4})([0-9]{2})([0-9]{2})([0-9]{2})([0-9]{2})([0-9]{2})$', 'tokens' );
		tstr = vertcat( tstr{:} );
		
		% Separate time components
		tvec{tf_long}    = tstr{:}{1};
		tvec{tf_long, 2} = tstr{:}{2};
		tvec{tf_long, 3} = tstr{:}{3};
		tvec{tf_long, 4} = tstr{:}{4};
		tvec{tf_long, 5} = tstr{:}{5};
		tvec{tf_long, 6} = tstr{:}{6};
	end

%------------------------------------%
% Convert to TT2000                  %
%------------------------------------%

	% Convert to tt2000?
	%   - SPDFComputeTT2000 expects an Nx9 array
	%   - Convert TVEC to number, then append three zeros to the end of each row.
	if nargout() == 2
		tt2000 = spdfcomputett2000( [ cellfun(@str2num, tvec) repmat([0,0,0], nPts, 1) ] );
	end
end