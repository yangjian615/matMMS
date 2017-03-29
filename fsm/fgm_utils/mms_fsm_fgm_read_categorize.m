%
% Name
%   mms_fsm_fgm_read_categorize
%
% Purpose
%   Read FGM data.
%
% Calling Sequence
%   DATA = mms_fsm_fgm_read_categorize( T_PACKET, T, B_OMB, RATE, RANGE, PFMODE )
%     Categorize FGM data by using the Deck 32/64 flag PFMODE, data sampling rate
%     RATE, and the hi/lo-range status flag RANGE. PFMODE is provided on packet time-
%     scale T_PACKET and is expanded to the full time resolution T of the other data.
%     Data is returned in the structure DATA.
%
% Parameters
%   T_PACKET        in, required, type = int64 (cdf_time_tt2000)
%   T               in, required, type = int64 (cdf_time_tt2000)
%   B_OMB           in, required, type = 4xN single
%   RATE            in, required, type = 1xN uint8
%   RANGE           in, required, type = 1xN uint8
%   PFMODE          in, required, type = 1xN uint8
%
% Returns
%   DATA            out, required, type=struct
%                   Fields are:
%                       't'      - TT2000 epoch times
%                       'b'      - Magnetic field in DMPA coordinates
%                       'flag'   - Bit flag indicating operational mode:
%                                    1: unset = Slow           set = Fast/Brst
%                                    2: unset = Lo-Range       set = Hi-Range
%                                    3: unset = Deck64         set = Deck32     (dfg)
%                                    3: unset = ADCA           set = ADCB       (afg)
%                                    4: unset = Post-Perigee   set = Pre-Perigee
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-06-01      Written by Matthew Argall
%   2016-09-13      Lo/HiRange occurs only in slow survey while Pre/Post-Perigee occur
%                     only in slow survey, lo range. They are given their own bit to
%                     make all states unique. - MRA
%   2016-10-01      Removed unused TRANGE parameter. Improved documentation. Renamed
%                     from mms_fsm_bkgd_fgm_read_categorize to mms_fsm_fgm_read_categorize - MRA
%   2017-01-30      Switch bit setting for DEC23/64 (now DEC64=unset and DEC32=set) to be
%                     consistent with parent PFMODE setting. Lo/Hi-Range are the same bit.  - MRA
%
function fgm = mms_fsm_fgm_read_categorize( t_packet, t, b_omb, rate, range, pfmode )

%------------------------------------%
% Interpolate PFMODE                 %
%------------------------------------%
	% Convert to seconds
	t_ref        = min( [t_packet(1) t(1)] );
	t_packet_ssm = MrCDF_epoch2sse(t_packet, t_ref);
	t_ssm        = MrCDF_epoch2sse(t, t_ref);

	% Locate the data points with the packet times
	inds   = MrValue_Locate(t_packet_ssm, t_ssm, 'RoundUp', true);
	pfmode = pfmode(inds);

	% Clear the time
	clear t_packet_ssm t_ssm

%------------------------------------%
% Isolate Orbit                      %
%------------------------------------%

	% Find the first & last slow survey data point
	rate = round(rate);
	is0  = find( rate ==  8, 1, 'first' );

	% Brst
	if isempty(is0)
		is0 = 1;        % So that is0:if1 returns the entire orbit
		is1 = 0;        % So that is0:is1 returns an empty array.
		if0 = 1;
		if1 = length(t);
	
	else
		if0 = find( rate(is0:end) == 16, 1, 'first' );
		if1 = find( rate(if0:end) == 16, 1, 'last' );
		if1 = if1 + if0 - 1;
		if0 = if0 + is0 - 1;
		is1 = if0 - 1;
	end

	% Trim the data
	t      = t( is0:if1 );
	b_omb  = b_omb(:, is0:if1);
	range  = range( is0:if1 );
	rate   = rate( is0:if1 );
	pfmode = pfmode( is0:if1 );
	
	% Adjust the indices
	is0 = is0 - is0 + 1;
	is1 = is1 - is0 + 1;
	if0 = if0 - is0 + 1;
	if1 = if1 - is0 + 1;

%------------------------------------%
% Allocate Memory                    %
%------------------------------------%
	nSlow = length( t(1:is1) );
	nFast = length( t(if0:end) );
	nPts  = nSlow + nFast;
	flag  = zeros(1, nPts, 'uint8');

%------------------------------------%
% Fast & Slow (bit 1)                %
%------------------------------------%
	% Set bit:
	%   - 0 = Slow
	%   - 1 = Fast/brst
	flag(if0:end ) = bitset( flag(if0:end), 1 );

%------------------------------------%
% Hi & Lo range (bits 2 & 3)         %
%------------------------------------%

	% Hi/Lo range
	%   |---- lo ----|-- hi --|---- lo ----|
	iHi0 = find(range == 1, 1, 'first');
	iHi1 = find(range == 1, 1, 'last');
	
	% Set bit:
	%   - 0 = Lo-Range
	%   - 2 = Hi-Range
%	flag( 1:iHi0-1     ) = bitset( flag(1:iHi0-1),     2 );
	flag( iHi0:iHi1    ) = bitset( flag(iHi0:iHi1),    2 );
%	flag( iHi1+1:nSlow ) = bitset( flag(iHi1+1:nSlow), 2 );

%------------------------------------%
% Pre- & Post-Perigee (bits 4 & 5)   %
%------------------------------------%
	
	% Set bit:
	%   - 0 = Post-Perigee
	%   - 4 = Pre-Perigee
	flag( 1:iHi0-1     ) = bitset( flag(1:iHi0-1),     4, 0 );
	flag( iHi1+1:nSlow ) = bitset( flag(iHi1+1:nSlow), 4, 1 );

%------------------------------------%
% DEC 32 & 64 (bit 3)                %
%------------------------------------%

	% Find first point
	i64_0 = find(pfmode == 0, 1, 'first');
	i32_0 = find(pfmode == 1, 1, 'first');

	% DEC32 and DEC64 should change at the outbound range change (hi to lo).
	%   - Range change con occur in the middle of a packet.
	%   - Mag team accounts for this in data products above L1A.
	%   - Adjust DEC32\64 change to match that of hi/lo range change.
	if ~isempty(i64_0) && ~isempty(i32_0)
		if i64_0 > i32_0
			dt    = double( t(i64_0) - t(iHi1+1) ) * 1e-9;
			i32_1 = i64_0 - 1;
			i64_1 = nPts;
		else
			dt    = double( t(i32_0) - t(iHi1+1) ) * 1e-9;
			i64_1 = i32_0 - 1;
			i32_1 = nPts;
		end
		
		% Indicate delay
		if dt ~= 0
			% Print an informational message
			mrfprintf( 'logtext', 'Adjusting DEC32/64 to Range Change:' );
			mrfprintf( 'logtext', '   dt = %d min %7.3f sec', round(dt / 60.0), sign(dt) * mod( abs(dt), 60.0) );

			% Make sure there were no other weird range changes
			%   - Should be +/- 1 packet length (8s)
			assert(dt < 8.5, 'Deck 32\64 change not when expected.')
			
			% Adjust DEC32\64 time
			if i64_0 > i32_0
				i32_1 = iHi1;           % [1, iHi]
				i64_0 = iHi1 + 1;       % [iHi+1, end]
				i64_1 = nPts;
			else
				i64_1 = iHi1;           % [1, iHi]
				i32_0 = iHi1 + 1;       % [iHi+1, end]
				i32_1 = nPts;
			end
		end
	
	% Only DEC32 available
	elseif isempty(i64_0)
		i32_1 = length(t);
		i64_1 = 0;              % So that i64_0:i64_1 returns the empty array
	
	% Only DEC64 available
	elseif isempty(i32_0)
		i32_1 = 0;              % So that i32_0:i32_1 returns the empty array
		i64_1 = length(t);
	
	% Unknown DEC
	else
		error('Neither Deck 32 nor Deck 64 are present.');
	end
	
	% Set bit:
	%   - 0 = DEC64
	%   - 3 = DEC32
	flag(i32_0:i32_1) = bitset( flag(i32_0:i32_1), 3 );

%------------------------------------%
% Output                             %
%------------------------------------%
	
	% Create the output structure
	fgm = struct( 't',    t,     ...
	              'sr',   rate,  ...
	              'b',    b_omb, ...
	              'flag', flag );
end