%
% Name
%   mms_fsm_fgm_read_l1a
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
%   2015-10-01      Renamed from mms_fsm_bkgd_fgm_read_l1a to mms_fsm_fgm_read_l1a. - MRA
%
function [pfmode, t] = mms_fsm_fgm_read_l1a(files, trange)

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

	% Variable names
	pf_name = [sc '_' instr '_' 'pfmode'];

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

%------------------------------------%
% Burst Data                         %
%------------------------------------%
	% Burst data files do not suffer from overlapping packet times
	%   - TODO: Is this true? What about adjacent burst files?
	if nargin() == 1 || isempty(trange{1})
	
		% Read data
		[pfmode, t] = MrCDF_nRead( files, pf_name );

%------------------------------------%
% Survey Data                        %
%------------------------------------%
	else
		% Declare variables
		t      = [];
		t_ssm  = [];
		pfmode = [];
		nToss  = 0;

	%------------------------------------%
	% Read Each File                     %
	%------------------------------------%
		for ii = 1 : length(files)
			% Read data
			[pf_temp, t_temp]    = MrCDF_nRead( files{ii}, pf_name );
			[~, ~, ~, ~, fstart] = mms_dissect_filename( files{ii} );
			[~, fstart]          = mms_parse_time(fstart);

			% Exclude time outside of the current day
			ssm_temp = MrCDF_epoch2ssm(t_temp, fstart);
			iKeep    = find( ssm_temp >= 0.0 & ssm_temp < 86400.0 );
			nToss    = nToss + (length(t_temp) - length(iKeep));
		
			% Store the data
			t_ssm  = [ t_ssm  ssm_temp(iKeep) ];
			t      = [ t      t_temp(iKeep) ];
			pfmode = [ pfmode pf_temp(iKeep) ];
		end

	%------------------------------------%
	% Overlapping Packets                %
	%------------------------------------%
	
		% The sample interval is 16s (8s) for slow (fast) survey.
		%   - There are overlapping points
		iToss         = find( diff(t_ssm) < 4.0 );
		t_ssm(iToss)  = [];
		t(iToss)      = [];
		pfmode(iToss) = [];
	
		% Are there more weird points remaining?
		iOdd = find( diff(t_ssm) < 8.0 );
	
		% Indicate how many points were removed
		mrfprintf('logwarn', '%i overlapping points removed.', nToss);
		mrfprintf('logwarn', '%i points with sample rate < 4.0s', length(iToss));
		mrfprintf('logwarn', '%i points with sample rate < 8.0s', length(iOdd));

	%------------------------------------%
	% Sort & Trim the Data               %
	%------------------------------------%
	
		% Sort the data
		[t, iSort] = sort(t);;
		pfmode     = pfmode(iSort);

		% Extend the time interval by 5 minutes on either side.
		%   - The orbit times are not exact.
		tt2000_range = MrCDF_Epoch_Parse( trange, 'CDF_TIME_TT2000');
		t0           = tt2000_range(1) - int64(5*60*1e9);
		t1           = tt2000_range(2) + int64(5*60*1e9);

		% Restrict times
		iKeep  = find( t >= t0 & t <= t1 );
		t      = t(iKeep);
		pfmode = pfmode(iKeep);
	end
end