%
% Name
%   mms_fsm_bkgd_fgm_categorize
%
% Purpose
%   Read FGM data.
%
% Calling Sequence
%   DATA = mms_fsm_bkgd_fgm_read(FILES, TRANGE)
%     Read and categorize FGM data based on the following factors: 1) slow/fast srvy,
%     2) hi/lo-range 3) before/after perigee, 4) deck 32/64. L1A file names F_L1A are
%     required to read deck 32/64 flags, while all other information comes from the
%     L2Pre files with names given by F_L2PRE. Data is read over the time interval
%     TRANGE, which must be formatted as 'yyyy-mm-ddTHH:MM:SS'. Categorized data is
%     returned in the structure DATA.
%
% Parameters
%   F_L1A           in, required, type = char/cell
%   F_L2PRE         in, required, type = char/cell
%   TRANGE          in, required, type = 1x2 cell
%
% Returns
%   DATA            out, required, type=struct
%                   Fields are:
%                       't'      - TT2000 epoch times
%                       'b'      - Magnetic field in DMPA coordinates
%                       'flag'   - Bit flag indicating operational mode:
%                                    0 = Slow    1 = Fast
%                                    0 = Lo-     2 = Hi-Range
%                                    0 = Before  4 = After perigee (for slow survey)
%                                    0 = Deck32  8 = Deck64
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-06-01      Written by Matthew Argall
%
function fgm = mms_fsm_bkgd_fgm_read( l1a_files, l2pre_files, trange )

%------------------------------------%
% Read Data                          %
%------------------------------------%

	% Read the L1A files
	[pfmode, t_packet] = mms_fsm_bkgd_fgm_read_l1a(l1a_files, trange);
	
	% Read the L2Pre files
	[b_bcs, range, rate, t] = mms_fsm_bkgd_fgm_read_l2pre(l2pre_files, trange);

%------------------------------------%
% Rotate BCS -> OMB                  %
%------------------------------------%
	zMPA = [0; 0; 1];
	mrfprintf('logwarn', 'Using zMPA = [0; 0; 1].');

	% Create transformation matrix
	bcs2smpa = mms_fg_xbcs2smpa(zMPA);
	omb2smpa = mms_fg_xomb2smpa();
	bcs2omb  = bcs2smpa * omb2smpa';

	% Rotate data
	b_omb = mrvector_rotate( bcs2omb, b_bcs(1:3,:) );
	
	clear b_bcs

%------------------------------------%
% Categorize Data                    %
%------------------------------------%
	
	% Categorize the data
	fgm = mms_fsm_bkgd_fgm_read_categorize( trange, t_packet, t, b_omb, rate, range, pfmode );
	
	% Clear all of the data
	clear t_packet t b_omb rate range pfmode

%------------------------------------%
% Isolate Modes                      %
%------------------------------------%
	
	% Flag
	%   - Bit 1 set:     Fast and Brst
	%   - Bit 1 not set: Slow 
	%   - iFast serves as iBrst when iSlow is empty
	iFast = find( bitget(fgm.flag, 1, 'uint8') );
	iSlow = find( bitget(fgm.flag, 1, 'uint8') == 0 );

%------------------------------------%
% High-Pass Filter                   %
%------------------------------------%
	% Filter parameters
	f0 = 0.5;           % bass-cut frequency

	% Slow
	if ~isempty(iSlow)
		% Create filter
		%   - filefilt requires doubles
		sr    = double( fgm.sr( iSlow(1) ) );
		fN    = sr / 2.0;
		[b,a] = butter(9, f0/fN, 'high');

		% Apply filter
		fgm.b(1,iSlow) = filtfilt(b, a, fgm.b(1,iSlow));
		fgm.b(2,iSlow) = filtfilt(b, a, fgm.b(2,iSlow));
		fgm.b(3,iSlow) = filtfilt(b, a, fgm.b(3,iSlow));
	end

	% Fast
	if ~isempty(iFast)
		% Create filter
		sr    = double( fgm.sr( iFast(1) ) );
		fN    = sr / 2.0;
		[b,a] = butter(9, f0/fN, 'high');
		
		% Apply filter
		fgm.b(1,iFast) = filtfilt(b, a, fgm.b(1,iFast));
		fgm.b(2,iFast) = filtfilt(b, a, fgm.b(2,iFast));
		fgm.b(3,iFast) = filtfilt(b, a, fgm.b(3,iFast));
	end
end


%
% Name
%   mms_fsm_bkgd_fgm_read_l1a
%
% Purpose
%   Read FGM data.
%
% Calling Sequence
%   [PFMODE, T] = mms_fsm_bkgd_fgm_read_l1a(FILES, TRANGE)
%     Read data files with names FILES over the course of an orbit, outlined by
%     the time range TRANGE, formatted as 'yyyy-mm-ddTHH:MM:SS'. Return the flag
%     PFFMODE that indicates whether deck 32 or 64 is in use, and also its time
%     tags, which are at the beginning of each packet.
%
% Parameters
%   FILES           in, required, type = char/cell
%   TRANGE          in, required, type = 1x2 cell
%
% Returns
%   PFMODE          out, required, type=float
%   T               out, required, type=int64
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-06-01      Written by Matthew Argall
%
function [pfmode, t] = mms_fsm_bkgd_fgm_read_l1a(files, trange)

%
% Each orbit should look like this
%
%                         | ------------ ROI  ------------ |
%    | ------ slow ------ | ------------ fast ------------ |
%    >--lo--|--hi--|--lo----------------------------------->
%    <--- 32/64 ---| -- 64/32 ----------------------------->
%
% We want to grab a little extra on either end
% just in case TS_RANGE or TF_RANGE are not precise.
%

%------------------------------------%
% Metadata                           %
%------------------------------------%

	% Make FILES a cell array for consistency.
	if ischar(files)
		files = { files };
	end
	
	% Parse the file name
	[sc, instr, mode, level] = mms_dissect_filename(files{1});

%------------------------------------%
% Extract Data for This Orbit        %
%------------------------------------%

	%
	% Srvy:
	%   Read an entire orbit. The input start time, taken from the
	%   SITL BSS FOM structure, is often off by several seconds.
	%   Make adjustments and use the sampling rate to distinguish
	%   between slow and fast survey.
	%
	%   It will be difficult to sort the files by time, since slow
	%   and fast survey files can have the same date. So, instead
	%   read the files in the order they are given, then sort the
	%   time afterward.
	%
	% Brst:
	%   Read the entire file.
	%

	% Variable names
	pf_name = [sc '_' instr '_' 'pfmode'];
	
	% Read data
	[pfmode, t] = MrCDF_nRead( files, pf_name );

%------------------------------------%
% Sort the Data                      %
%------------------------------------%
	
	if nargin == 2 && ~isempty(trange)
		% Extend the time interval by 5 minutes on either side.
		tt2000_range = MrCDF_Epoch_Parse( trange, 'CDF_TIME_TT2000');
		t0           = tt2000_range(1) - int64(5*60*1e9);
		t1           = tt2000_range(2) + int64(5*60*1e9);
	
		% Sort the data
		[t, isort] = sort(t);
		pfmode     = pfmode(isort);
	
		% Restrict times
		iKeep  = find( t >= t0 & t <= t1 );
		t      = t(iKeep);
		pfmode = pfmode(iKeep);
	end
end


%
% Name
%   mms_fsm_bkgd_fgm_read_l2pre
%
% Purpose
%   Read FGM L2Pre data.
%
% Calling Sequence
%   [B_OMB, RANGE, RATE, T] = mms_fsm_bkgd_fgm_read(FILES, TRANGE)
%     Read data files with names FILES over the course of an orbit, outlined by
%     the time range TRANGE, formatted as 'yyyy-mm-ddTHH:MM:SS'. Return data in
%     the structure DATA. Return the vector magnetic field on OMB coordinates,
%     the hi- and lo-range flag, the sampling rate, and the time tags.
%
% Parameters
%   FILES           in, required, type = char/cell
%   TRANGE          in, required, type = 1x2 cell
%
% Returns
%   B_OMB           out, required, type=3xN single
%   RANGE           out, required, type=int8
%   RATE            out, required, type=single
%   T               out, required, type=int64
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-06-01      Written by Matthew Argall
%
function [b_bcs, range, rate, t] = mms_fsm_bkgd_fgm_read_l2pre(files, trange)

%
% Each orbit should look like this
%
%                         | ------------ ROI  ------------ |
%    | ------ slow ------ | ------------ fast ------------ |
%    >--lo--|--hi--|--lo----------------------------------->
%    <--- 32/64 ---| -- 64/32 ----------------------------->
%
% We want to grab a little extra on either end
% just in case TS_RANGE or TF_RANGE are not precise.
%

%------------------------------------%
% Metadata                           %
%------------------------------------%
	% Make FILES a cell array for consistency.
	if ischar(files)
		files = { files };
	end

	% Parse the file name
	[sc, instr, mode, level] = mms_dissect_filename(files{1});

%------------------------------------%
% Extract Data for This Orbit        %
%------------------------------------%

	%
	% Srvy:
	%   Read an entire orbit. The input start time, taken from the
	%   SITL BSS FOM structure, is often off by several seconds.
	%   Make adjustments and use the sampling rate to distinguish
	%   between slow and fast survey.
	%
	% Brst:
	%   Read the entire file.
	%
	
	% Extend the time interval by 5 minutes on either side.
	trange_ext = {'', ''};
	if nargin == 2 && ~isempty(trange)
		tt2000_range  = MrCDF_Epoch_Parse( trange, 'CDF_TIME_TT2000');
		tt2000_ext    = zeros(1, 2, 'int64');
		tt2000_ext(1) = tt2000_range(1) - int64(5*60*1e9);
		tt2000_ext(2) = tt2000_range(2) + int64(5*60*1e9);
		trange_ext    = MrCDF_Epoch_Encode(tt2000_ext);
	end

	% Variable names
	b_bcs_name    = [sc '_' instr '_' mode '_' level '_bcs'];
	range_name    = [sc '_' instr '_' mode '_' level '_hirange'];
	rate_name     = [sc '_' instr '_' mode '_' level '_rate'];
	l1a_mode_name = [sc '_' instr '_' mode '_' level '_l1a_mode'];
	
	% Read data
	[b_bcs, t] = MrCDF_nRead(files, b_bcs_name,    'sTime', trange_ext{1}, 'eTime', trange_ext{2});
	range      = MrCDF_nRead(files, range_name,    'sTime', trange_ext{1}, 'eTime', trange_ext{2}); % 1=hi 0=lo
	rate       = MrCDF_nRead(files, rate_name,     'sTime', trange_ext{1}, 'eTime', trange_ext{2});
	mode       = MrCDF_nRead(files, l1a_mode_name, 'sTime', trange_ext{1}, 'eTime', trange_ext{2});
end


%
% Name
%   mms_fsm_bkgd_fgm_categorize
%
% Purpose
%   Read FGM data.
%
% Calling Sequence
%   DATA = mms_fsm_bkgd_fgm_read(FILES, TRANGE)
%     Read and categorize FGM data based on the following factors: 1) slow/fast srvy,
%     2) hi/lo-range 3) before/after perigee, 4) deck 32/64. L1A file names F_L1A are
%     required to read deck 32/64 flags, while all other information comes from the
%     L2Pre files with names given by F_L2PRE. Data is read over the time interval
%     TRANGE, which must be formatted as 'yyyy-mm-ddTHH:MM:SS'. Categorized data is
%     returned in the structure DATA.
%
% Parameters
%   F_L1A           in, required, type = char/cell
%   F_L2PRE         in, required, type = char/cell
%   TRANGE          in, required, type = 1x2 cell
%
% Returns
%   DATA            out, required, type=struct
%                   Fields are:
%                       't'      - TT2000 epoch times
%                       'b'      - Magnetic field in DMPA coordinates
%                       'flag'   - Bit flag indicating operational mode:
%                                    0 = Slow    1 = Fast
%                                    0 = Lo-     2 = Hi-Range
%                                    0 = Before  4 = After perigee (for slow survey)
%                                    0 = Deck32  8 = Deck64
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-06-01      Written by Matthew Argall
%
function fgm = mms_fsm_bkgd_fgm_read_categorize( trange, t_packet, t, b_omb, rate, range, pfmode )

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
	%   - 0 = Slow  1 = Fast/brst
	flag(if0:end ) = bitset( flag(if0:end), 1 );

%------------------------------------%
% Hi & Lo range (bits 2 & 3)         %
%------------------------------------%

	% Hi/Lo range
	%   |---- lo ----|-- hi --|---- lo ----|
	iHi0 = find(range == 1, 1, 'first');
	iHi1 = find(range == 1, 1, 'last');
	
	% Set bit:
	%   - 0 = Lo-   2 = Hi-Range
	%   - 0 = Pre-  4 = Post-Perigee
%	flag( 1:iHi0-1     ) = bitset( flag(1:iHi0-1),     0 );
	flag( iHi0:iHi1    ) = bitset( flag(iHi0:iHi1),    2 );
	flag( iHi1+1:nSlow ) = bitset( flag(iHi1+1:nSlow), 3 );

%------------------------------------%
% Deck 32 & 64 (bit 4)               %
%------------------------------------%

	% Find first point
	i64_0 = find(pfmode == 0, 1, 'first');
	i32_0 = find(pfmode == 1, 1, 'first');

	% Deck32 and Deck64 should change at the outbound range change (hi to lo).
	%   - Range change con occur in the middle of a packet.
	%   - Mag team accounts for this in data products above L1A.
	%   - Adjust Deck32\64 change to match that of hi/lo range change.
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
			mrfprintf( 'logtext', 'Adjusting Deck32/64 to Range Change:' );
			mrfprintf( 'logtext', '   dt = %d min %7.3f sec', round(dt / 60.0), sign(dt) * mod( abs(dt), 60.0) );

			% Make sure there were no other weird range changes
			%   - Should be +/- 1 packet length (8s)
			assert(dt < 8.5, 'Deck 32\64 change not when expected.')
			
			% Adjust Deck32\64 time
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
	
	% Only Deck 32 available
	elseif isempty(i64_0)
		i32_1 = length(t);
		i64_1 = 0;              % So that i64_0:i64_1 returns the empty array
	
	% Only Deck 64 available
	elseif isempty(i32_0)
		i32_1 = 0;              % So that i32_0:i32_1 returns the empty array
		i64_1 = length(t);
	
	% Unknown Deck
	else
		error('Neither Deck 32 nor Deck 64 are present.');
	end
	
	% Set bit:
	%   - 0 = Deck32   8 = Deck64
	flag(i64_0:i64_1) = bitset( flag(i64_0:i64_1), 4 );

%------------------------------------%
% Output                             %
%------------------------------------%
	
	% Create the output structure
	fgm = struct( 't',    t,     ...
	              'sr',   rate,  ...
	              'b',    b_omb, ...
	              'flag', flag );
end
