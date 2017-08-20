%
% Name
%   mms_fsm_scm_read_l1b
%
% Purpose
%   Read SCM data and sort it into fast and slow survey portions.
%
% Calling Sequence
%   DATA = mms_fsm_scm_read_l1b(FILES)
%     Read data files with names FILES over the course of an orbit, outlined by
%     the time range TRANGE, formatted as 'yyyy-mm-ddTHH:MM:SS'. Return data in
%     the structure DATA.
%
%   DATA = mms_fsm_scm_read_l1b(FILES, TRANGE)
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
%   2015-10-01      Renamed from mms_fsm_bkgd_scm_read_l1b to mms_fsm_scm_read_l1b. - MRA
%
function [b_123, t] = mms_fsm_scm_read_l1b(files, trange)

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
	if nargin == 2 && ~isempty( trange{1} )
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

%------------------------------------%
% Brst Data                          %
%------------------------------------%
	if nargin < 2 || isempty( trange{1} )
		[b_123, t] = MrCDF_nRead(files, b_123_name);
		
%------------------------------------%
% Srvy Data                          %
%------------------------------------%
	else
		
		% Allocate memory
		b_123 = zeros(3, 0, 'single');
		t     = [];
		nToss = 0;
		
		% Loop over each file
		for ii = 1 : length(files)
			% Read data
			[b_temp, t_temp] = MrCDF_Read(files{ii}, b_123_name, 'sTime', trange_ext{1}, 'eTime', trange_ext{2});
			
			
			% File start time for reference
			[~, ~, ~, ~, fstart] = mms_dissect_filename( files{ii} );
			[~, fstart]          = mms_parse_time(fstart);

			% Exclude time outside of the current day
			ssm_temp = MrCDF_epoch2ssm(t_temp, fstart);
			iKeep    = find( ssm_temp >= 0.0 & ssm_temp < 86400.0 );
			nToss    = nToss + (length(t_temp) - length(iKeep));

			% Save the data
			t    =  [t     t_temp(iKeep)   ];
			b_123 = [b_123 b_temp(:,iKeep) ];
		end
		
		% Warn about lost data
		mrfprintf('logwarn', '%i overlapping SCM data points.', nToss);
	end
end
