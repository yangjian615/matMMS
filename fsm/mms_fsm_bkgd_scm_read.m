%
% Name
%   mms_fsm_bkgd_scm_read
%
% Purpose
%   Read SCM data and sort it into fast and slow survey portions.
%
% Calling Sequence
%   DATA = mms_fsm_bkgd_scm_read(FILES)
%     Read data files with names FILES over the course of an orbit, outlined by
%     the time range TRANGE, formatted as 'yyyy-mm-ddTHH:MM:SS'. Return data in
%     the structure DATA.
%
%   DATA = mms_fsm_bkgd_scm_read(FILES, TRANGE)
%     Read data within the time range TRANGE, formatted as 'yyyy-mm-ddTHH:MM:SS'.
%     Return data in the structure DATA.
%
% Parameters
%   FILES           in, required, type = char/cell
%   TRANGE          in, optional, type = 1x2 cell
%
% Returns
%   DATA            out, required, type=struct
%                   Fields are:
%                       't'      - TT2000 epoch times
%                       'b'      - Magnetic field in 123 coordinates
%                       'flag'   - Flag indicating operational mode:
%                                    0 = Slow    1 = Fast
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-06-01      Written by Matthew Argall
%
function scm = mms_fsm_bkgd_scm_read(file, trange)

%------------------------------------%
% Read Data                          %
%------------------------------------%
	
	% Read the L2Pre files
	[b_123, t] = mms_fsm_bkgd_scm_read_l1b(file, trange);
	
	% Categorize the data
	scm = mms_fsm_bkgd_scm_read_categorize(t, b_123);
	
	% Clear the data
	clear t b_123

%------------------------------------%
% Interpolate Over NaNs              %
%------------------------------------%
	% Find NaNs
	tf_nan = isnan(scm.b(1,:));
	nNaN   = sum(tf_nan);
	if nNaN > 0
		% Issue warning
		mrfprintf('logwarn', '%d NaN values found. Inteprolating.', nNaN);
		
		% Convert time to doubles (required by interp1)
		t_ssm = MrCDF_epoch2ssm(scm.t);
		
		% Interpolate
		scm.b(1,tf_nan) = interp1( t_ssm(~tf_nan), scm.b(1,~tf_nan), t_ssm(tf_nan) );
		scm.b(2,tf_nan) = interp1( t_ssm(~tf_nan), scm.b(2,~tf_nan), t_ssm(tf_nan) );
		scm.b(3,tf_nan) = interp1( t_ssm(~tf_nan), scm.b(3,~tf_nan), t_ssm(tf_nan) );
		
		% Clear data
		clear t_ssm
	end

%------------------------------------%
% Isolate Modes                      %
%------------------------------------%
	
	% Flag
	%   - Bit 1 set:     Fast and Brst
	%   - Bit 1 not set: Slow 
	%   - iFast serves as iBrst when iSlow is empty
	iFast = find( bitget(scm.flag, 1, 'uint8') );
	iSlow = find( bitget(scm.flag, 1, 'uint8') == 0 );

%------------------------------------%
% High-Pass Filter                   %
%------------------------------------%
	% Filter parameters
	f0 = 0.5;           % bass-cut frequency

	% Slow
	if ~isempty(iSlow)
		% Create filter
		fN    = scm.sr( iSlow(1) ) / 2.0;
		[b,a] = butter(9, f0/fN, 'high');
		
		% Apply filter
		scm.b(1,iSlow) = filtfilt(b, a, double( scm.b(1,iSlow) ));
		scm.b(2,iSlow) = filtfilt(b, a, double( scm.b(2,iSlow) ));
		scm.b(3,iSlow) = filtfilt(b, a, double( scm.b(3,iSlow) ));
	end

	% Fast & Brst
	if ~isempty(iFast)
		% Create filter
		fN    = scm.sr( iFast(1) ) / 2.0;
		[b,a] = butter(9, f0/fN, 'high');
		
		% Apply filter
		scm.b(1,iFast) = filtfilt(b, a, double( scm.b(1,iFast) ));
		scm.b(2,iFast) = filtfilt(b, a, double( scm.b(2,iFast) ));
		scm.b(3,iFast) = filtfilt(b, a, double( scm.b(3,iFast) ));
	end
end


%
% Name
%   mms_fsm_bkgd_scm_read
%
% Purpose
%   Read SCM data and sort it into fast and slow survey portions.
%
% Calling Sequence
%   DATA = mms_fsm_bkgd_scm_read(FILES)
%     Read data files with names FILES over the course of an orbit, outlined by
%     the time range TRANGE, formatted as 'yyyy-mm-ddTHH:MM:SS'. Return data in
%     the structure DATA.
%
%   DATA = mms_fsm_bkgd_scm_read(FILES, TRANGE)
%     Read data within the time range TRANGE, formatted as 'yyyy-mm-ddTHH:MM:SS'.
%     Return data in the structure DATA.
%
% Parameters
%   FILES           in, required, type = char/cell
%   TRANGE          in, optional, type = 1x2 cell
%
% Returns
%   DATA            out, required, type=struct
%                   Fields are:
%                       't'      - TT2000 epoch times
%                       'b'      - Magnetic field in 123 coordinates
%                       'flag'   - Flag indicating operational mode:
%                                    0 = Slow    1 = Fast
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-06-01      Written by Matthew Argall
%
function [b_123, t] = mms_fsm_bkgd_scm_read_l1b(files, trange)

%
% Each orbit should look like this
%
%    | ---- slow ---- | -------- fast -------- |
%
% Or
%
%    | ---------------- srvy ----------------- |
%
% We want to grab a little extra on either end
% just in case TS_RANGE or TF_RANGE are not precise.
%


%------------------------------------%
% Metadata                           %
%------------------------------------%
	
	% Make sure FILES is a cell array for consistency
	if ischar(files)
		files = { files };
	end

	% Parse the file name
	[sc, instr, mode, level] = mms_dissect_filename(files{1});

%------------------------------------%
% Read Data                          %
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

	switch mode
		case 'srvy'
			optdesc = 'scsrvy';
		case 'brst'
			optdesc = 'scb';
		case 'slow'
			optdesc = 'scs';
		case 'fast'
			optdesc = 'scf';
		otherwise
			error(['Unrecognized telemetry mode: "' mode '".']);
	end

	% Variable names
	b_123_name = [sc '_' instr '_acb_scm123_' optdesc '_' mode '_' level];
	
	% Read data
	[b_123, t] = MrCDF_nRead(files, b_123_name, 'sTime', trange_ext{1}, 'eTime', trange_ext{2});
end


%
% Name
%   mms_fsm_bkgd_scm_read
%
% Purpose
%   Read SCM data and sort it into fast and slow survey portions.
%
% Calling Sequence
%   DATA = mms_fsm_bkgd_scm_read(FILES)
%     Read data files with names FILES over the course of an orbit, outlined by
%     the time range TRANGE, formatted as 'yyyy-mm-ddTHH:MM:SS'. Return data in
%     the structure DATA.
%
%   DATA = mms_fsm_bkgd_scm_read(FILES, TRANGE)
%     Read data within the time range TRANGE, formatted as 'yyyy-mm-ddTHH:MM:SS'.
%     Return data in the structure DATA.
%
% Parameters
%   FILES           in, required, type = char/cell
%   TRANGE          in, optional, type = 1x2 cell
%
% Returns
%   DATA            out, required, type=struct
%                   Fields are:
%                       't'      - TT2000 epoch times
%                       'b'      - Magnetic field in 123 coordinates
%                       'flag'   - Flag indicating operational mode:
%                                    0 = Slow    1 = Fast
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-06-01      Written by Matthew Argall
%
function data = mms_fsm_bkgd_scm_read_categorize(t, b_123)

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
	sr    = [sr(is0:if1-1) sr(end)];
	t     = t(is0:if1);
	b_123 = b_123(:, is0:if1);
	
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
	
	% Create a flag
	%   0 = Slow  1 = Fast Survey
	if is1 == 0
		flag(:) = uint8( 1 );
	else
		flag(is0:is1) = uint8( 0 );
		flag(if0:if1) = uint8( 1 );
	end

%------------------------------------%
% Collect Data                       %
%------------------------------------%
	
	% Data structure
	data = struct( 't',    t,     ...
	               'b',    b_123, ...
	               'sr',   sr,    ...
	               'flag', flag );
end