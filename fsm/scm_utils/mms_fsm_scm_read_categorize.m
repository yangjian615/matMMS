%
% Name
%   mms_fsm_scm_read_categorize
%
% Purpose
%   Read SCM data and sort it into fast and slow survey portions.
%
% Calling Sequence
%   DATA = mms_fsm_scm_read_categorize(T, B)
%     Categorize SCM magnetic field data by mode (slow/fast/brst).
%     Use FGM information obtained from mms_fsm_bkgd_fgm_read to apply the deck 32/64
%     flag to SCM data.
%
% Parameters
%   T               in, required, type = uint64
%   B               in, required, type = 3xN single
%
% Returns
%   DATA            out, required, type=struct
%                   Fields are:
%                       't'      - TT2000 epoch times
%                       'b'      - Magnetic field in 123 coordinates
%                       'flag'   - Flag indicating operational mode:
%                                    1: unset = Slow    set = Fast
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-06-01      Written by Matthew Argall
%   2016-09-07      Added the FGM parameter to categorize by deck or range. - MRA
%   2016-10-01      Renamed from mms_fsm_bkgd_scm_read_categorize to mms_fsm_scm_read_categorize. - MRA
%   2016-10-21      Removed the FGM parameter (SCM is not affected by DAC 32/64). - MRA
%
function data = mms_fsm_scm_read_categorize(t, b)

%------------------------------------%
% Isolate Orbit                      %
%------------------------------------%

	%
	% After Sept 1, 2015 (?), SCM sampling rate in slow survey
	% was changed to be 32 S/s, the same as in fast survey
	%
	
	% Compute the sampling rate
	dt = double( diff(t) ) * 1e-9;
	sr = round(1 ./ dt);
	
	% First and last srvy data point
	is0 = find( sr == 16, 1, 'first');
	if1 = find( sr == 32, 1, 'last');
	
	% Brst mode
	if isempty(is0) && isempty(if1)
		is0 = 1;        % So that is0:if1 returns the entire orbit
		is1 = 0;        % So that is0:is1 returns an empty array.
		if0 = 1;
		if1 = length(t);
	
	% Srvy mode (slow == fast)
	elseif isempty(is0)
		is0 = 1;        % So that is0:if1 returns the entire orbit
		is1 = 0;        % So that is0:is1 returns an empty array.
		if0 = 1;
		if1 = length(t);
	
	% Srvy mode
	else
		% Fast and Slow boundary
		if0 = find( sr == 32, 1, 'first' );
		is1 = if0 - 1;
	end
	
	% Trim the data
	sr = [sr(is0:if1-1) sr(end)];
	t  = t(is0:if1);
	b  = b(:, is0:if1);
	
	% Adjust the indices
	is0 = is0 - is0 + 1;
	is1 = is1 - is0 + 1;
	if0 = if0 - is0 + 1;
	if1 = if1 - is0 + 1;

%------------------------------------%
% Create the Flag                    %
%------------------------------------%
	nSlow = length( t(is0:is1) );
	nFast = length( t(if0:if1) );
	flag  = zeros( 1, nSlow + nFast, 'uint8' );

%------------------------------------%
% Slow / Fast / Brst (bit 1)         %
%------------------------------------%
	
	% Create a flag
	%   0 = Slow  1 = Fast Survey
	flag(if0:if1) = bitset( flag(if0:if1), 1 );

%------------------------------------%
% Collect Data                       %
%------------------------------------%
	
	% Data structure
	data = struct( 't',    t,     ...
	               'b',    b,     ...
	               'sr',   sr,    ...
	               'flag', flag );
end