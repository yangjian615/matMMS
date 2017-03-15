%
% Name
%   mms_fsm_fgm_read_l2temp
%
% Purpose
%   Read FGM L2 intermediate data.
%
% Calling Sequence
%   [B_OMB, HIRANGE, RATE, T] = mms_fsm_fgm_read_l2temp(FILES, TRANGE)
%     Read data files with names FILES over the course of an orbit, outlined by
%     the time range TRANGE, formatted as 'yyyy-mm-ddTHH:MM:SS'. Return data in
%     the structure DATA. Return the vector magnetic field in OMB coordinates,
%     the hi-range flag, the sampling rate, and the time tags.
%
% Parameters
%   FILES           in, required, type = char/cell
%   TRANGE          in, required, type = 1x2 cell
%
% Returns
%   B_OMB           out, required, type=3xN single
%   HIRANGE         out, required, type=int8
%   RATE            out, required, type=single
%   T               out, required, type=int64
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2017-03-14      Written by Matthew Argall
%
function [b_omb, hirange, rate, t] = mms_fsm_fgm_read_l2temp(files, trange)

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
	b_name        = [sc '_' instr '_b_omb_'    mode '_' level];
	hirange_name  = [sc '_' instr '_hirange_'  mode '_' level];
	rate_name     = [sc '_' instr '_rate_'     mode '_' level];
	l1a_mode_name = [sc '_' instr '_l1a_mode_' mode '_' level];
	
%------------------------------------%
% Burst Data                         %
%------------------------------------%
	if nargin < 2 || isempty(trange{1})
		[b_omb, t] = MrCDF_nRead(files, b_name);
		hirange    = MrCDF_nRead(files, hirange_name); % 1=hi 0=lo
		rate       = MrCDF_nRead(files, rate_name);
		mode       = MrCDF_nRead(files, l1a_mode_name);
		
%------------------------------------%
% Survey Data                        %
%------------------------------------%
	else 
		% Declare variables
		t       = [];
		b_omb   = zeros(4, 0, 'single');
		hirange = [];
		rate    = [];
		mode    = [];
		nToss   = 0;
		
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
			hirange_temp     = MrCDF_Read(files{ii}, hirange_name,  'sTime', trange_ext{1}, 'eTime', trange_ext{2}); % 1=hi 0=lo
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
			t       = [t        t_temp(iKeep)       ];
			b_omb   = [b_omb    b_temp(:,iKeep)     ];
			hirange = [hirange  hirange_temp(iKeep) ];
			rate    = [rate     rate_temp(iKeep)    ];
			mode    = [mode     mode_temp(iKeep)    ];
		end
		
		% Warn about lost data
		mrfprintf('logwarn', '%i overlapping FGM data points.', nToss);
	end
end