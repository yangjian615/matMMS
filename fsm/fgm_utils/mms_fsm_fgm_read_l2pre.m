%
% Name
%   mms_fsm_fgm_read_l2pre
%
% Purpose
%   Read FGM L2Pre data.
%
% Calling Sequence
%   [B_OMB, RANGE, RATE, T] = mms_fsm_fgm_read_l2pre(FILES, TRANGE)
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
%   2015-10-01      Renamed from mms_fsm_bkgd_fgm_read_l2pre to mms_fsm_fgm_read_l2pre. - MRA
%
function [b_bcs, range, rate, t] = mms_fsm_fgm_read_l2pre(files, trange)

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

	% Variable names
	b_bcs_name    = [sc '_' instr '_' mode '_' level '_bcs'];
	range_name    = [sc '_' instr '_' mode '_' level '_hirange'];
	rate_name     = [sc '_' instr '_' mode '_' level '_rate'];
	l1a_mode_name = [sc '_' instr '_' mode '_' level '_l1a_mode'];

%------------------------------------%
% Burst Data                         %
%------------------------------------%
	if nargin < 2 || isempty(trange{1})
		[b_bcs, t] = MrCDF_nRead(files, b_bcs_name);
		range      = MrCDF_nRead(files, range_name); % 1=hi 0=lo
		rate       = MrCDF_nRead(files, rate_name);
		mode       = MrCDF_nRead(files, l1a_mode_name);

%------------------------------------%
% Survey Data                        %
%------------------------------------%
	else 
		% Declare variables
		t     = [];
		b_bcs = zeros(4, 0, 'single');
		range = [];
		rate  = [];
		mode  = [];
		nToss = 0;
		
		% Extend the time interval by 5 minutes on either side.
		tt2000_range  = MrCDF_Epoch_Parse( trange, 'CDF_TIME_TT2000');
		tt2000_ext    = zeros(1, 2, 'int64');
		tt2000_ext(1) = tt2000_range(1) - int64(5*60*1e9);
		tt2000_ext(2) = tt2000_range(2) + int64(5*60*1e9);
		trange_ext    = MrCDF_Epoch_Encode(tt2000_ext);
		
		% Step through each file
		for ii = 1 : length(files)
			% Read data
			[b_temp, t_temp] = MrCDF_Read(files{ii}, b_bcs_name,    'sTime', trange_ext{1}, 'eTime', trange_ext{2});
			range_temp       = MrCDF_Read(files{ii}, range_name,    'sTime', trange_ext{1}, 'eTime', trange_ext{2}); % 1=hi 0=lo
			rate_temp        = MrCDF_Read(files{ii}, rate_name,     'sTime', trange_ext{1}, 'eTime', trange_ext{2});
			mode_temp        = MrCDF_Read(files{ii}, l1a_mode_name, 'sTime', trange_ext{1}, 'eTime', trange_ext{2});
			
			% File start time for reference
			[~, ~, ~, ~, fstart] = mms_dissect_filename( files{ii} );
			[~, fstart]          = mms_parse_time(fstart);

			% Exclude time outside of the current day
			ssm_temp = MrCDF_epoch2ssm(t_temp, fstart);
			iKeep    = find( ssm_temp >= 0.0 & ssm_temp < 86400.0 );
			nToss    = nToss + (length(t_temp) - length(iKeep));

			% Save the data
			t    =  [t         t_temp(iKeep)   ];
			b_bcs = [b_bcs     b_temp(:,iKeep) ];
			range = [range range_temp(iKeep)   ];
			rate  = [rate   rate_temp(iKeep)   ];
			mode  = [mode   mode_temp(iKeep)   ];
		end
		
		% Warn about lost data
		mrfprintf('logwarn', '%i overlapping FGM data points.', nToss);
	end
end