%
% Name
%   mms_fsm_bkgd_process
%
% Purpose
%   Generate calibration data for the fluxgate-searchcoil merged (FSM) data
%   product.
%
% Calling Sequence
%   [] = mms_fsm_bkgd_process(SC, MODE, TSTART)
%     Read burst mode FGM and SCM magnetometer data, compute noise floor
%     parameters, and write the results to FILE_FGM and FILE_SCM. TSTART
%     must be formatted as yyyymmdd or yyyymmddHHMMSS
%
%   [] = mms_fsm_bkgd_process(SC, MODE, TSTART, TEND)
%     Read slow, fast, or srvy mode FGM and SCM magnetometer data,
%     compute noise floor parameters, and write the results to 
%     FILE_FGM and FILE_SCM. TSTART and TEND are intended to outline
%     an orbit of data, starting with slow survy, and ending after
%     fast survey.
%
%   [] = mms_fsm_bkgd_process(__, 'ParamName', ParamValue)
%     Any parameter name-value pairs given below may be used with any of
%     the above
%
% Parameters
%   SC              in, required, type = char
%   MODE            in, required, type = char
%   TSTART          in, required, type = char
%   TEND            in, optional, type = char, default = ''
%   'Duration'      in, optional, type = double, default = 20.0
%                   Number of seconds per calibration interval.
%   'FGMInstr'      in, optional, type = double, default = 'dfg'
%                   FGM instrument to use. Options are {'afg' | 'dfg'}
%   'NoLog'         in, optional, type = boolean, default = false
%                   If true, no log file will be created.
%
% Returns
%   FGM_FILE        out, required, type=string
%   SCM_FILE        out, required, type=string
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2016-07-15      Written by Matthew Argall
%
function [] = mms_fsm_bkgd_process(sc, mode, tstart, tend, varargin)

	% Declare global attributes
	global cal_path_root data_path_root dropbox_root hk_root log_path_root unh_data_root
	
	% Initialize
	mms_fsm_init();
	
	% Defaults
	tf_log    = true;
	fgm_instr = 'dfg';
	
	nOptArgs = length(varargin);
	for ii = 1 : 2 : nOptArgs
		switch varargin{ii}
			case 'NoLog'
				tf_log = ~varargin{ii+1};
			case 'FGMInstr'
				fgm_instr = varargin{ii+1};
			otherwise
				% Do nothing. Pass others on to *_sdc()
		end
	end

%------------------------------------%
% Check Inputs                       %
%------------------------------------%
	
	% Constants
	instr = 'fsm';
	level = 'l2plus';

	assert( min( ismember(sc,   {'mms1', 'mms2', 'mms3', 'mms4'}) ) == 1, 'Invalid spacecraft.' );
	assert( min( ismember(mode, {'slow', 'fast', 'srvy', 'brst'}) ) == 1, 'Invalid mode.' );
	
	% Defaults
	if nargin < 4
		tend = {};
	end
	
	% Turn strings into cells for consistency
	if ischar(sc)
		sc = { sc };
	end
	if ischar(mode)
		mode = { mode };
	end
	
	% Reformat the times
	start_vec = mms_parse_time(tstart);
	end_vec   = mms_parse_time(tend);

	% Convert to ISO
	tstart = datestr( cellfun(@str2num, start_vec), 'yyyy-mm-ddTHH:MM:SS');
	tend   = datestr( cellfun(@str2num, end_vec),   'yyyy-mm-ddTHH:MM:SS');

%------------------------------------%
% Create Log File                    %
%------------------------------------%
	% Name of log file
	if tf_log 
		% Calculate the current time
		%   - Append a random number in case processes start at same time
		date = datestr( now(), 'yyyymmdd' );
		time = datestr( now(), 'HHMMSS' );
		rnum = sprintf( '%06i', fix( rand(1,1,'single')*1e6 ) );
		
		% Name of log file
		logDir = fullfile(log_path_root, 'batch_logs');
		if exist(logDir, 'dir') ~= 7
			mkdir(logDir);
		end
		fLog = fullfile(logDir, ['mms_fsm_bkgd_' date '_' time '_' rnum '.log']);
	else
		fLog = 'stderr';
	end
	
	% Create log file
	oLog = MrLogFile(fLog);

%------------------------------------%
% Step Through Each File Type        %
%------------------------------------%
	% Statistics
	nalloc   = 100;
	files    = cell(2, nalloc);
	status   = zeros(nalloc, 'uint8');
	telapsed = zeros(nalloc);
	count    = 0;
	
	% Start run
	t0 = now();
	oLog.AddText( 'Processing FSM background calibrations.' );
	oLog.AddText( ['Starting at: ' datestr(t0, 'yyyy-mm-ddTHH:MM:SS') ] );
	oLog.AddText( '' );

	% Loop
	for ii = 1 : length(sc)
	for jj = 1 : length(mode)

	%------------------------------------%
	% Log Info                           %
	%------------------------------------%
		oLog.AddText( '------------------------------------------------------' );
		oLog.AddText( ['PROCESSING ' sc{ii} ' ' mode{jj} ' ' tstart ' - ' tend] );
		oLog.AddText( '' );

	%------------------------------------%
	% Time Interval                      %
	%------------------------------------%
		%
		% Brst:
		%   L1A files should provide the most complete dataset. Higher level
		%   data that has not been processed yet will cause errors to tell us
		%   that processing needs to be completed. Time intervals are taken
		%   directly from the file names.
		%
		% Srvy
		%   We will be processing on a per-orbit basis, so the time interval
		%   will be defined by orbit times, not by file times. We use the STIL
		%   ROI times to generate orbit intervals.
		%

		% BRST files
		if strcmp(mode{jj}, 'brst')
			% FGM
			fgm_files = mms_find_file( sc{ii}, fgm_instr, mode{jj}, 'l1a', ...
			                           'Dropbox', dropbox_root,            ...
			                           'SDCroot', data_path_root,          ...
			                           'TStart' , tstart,                  ...
			                           'TEnd',    tend );
			
			% SCM
			scm_files = mms_find_file( sc{ii}, 'scm', mode{jj}, 'l1a', ...
			                           'Dropbox', dropbox_root,        ...
			                           'SDCroot', data_path_root,      ...
			                           'OptDesc', 'scb',               ...
			                           'TStart',  tstart,              ...
			                           'TEnd',    tend );
			
			% File start times
			[~, ~, ~, ~, fgm_fstart] = mms_dissect_filename(fgm_files);
			[~, ~, ~, ~, scm_fstart] = mms_dissect_filename(scm_files);
			
			if ischar(fgm_fstart)
				fgm_fstart = { fgm_fstart };
			end
			if ischar(scm_fstart)
				scm_fstart = { scm_fstart };
			end

			% Take only unique values
			times = unique( [ fgm_fstart scm_fstart ] );
			times = MrTimeParser(times, '%Y%M%d%H%m%S', '%Y-%M-%dT%H:%m:%S');

		% SRVY files
		else
			% Times of each orbit
			tt2000_orbit = mms_bss_roi_get(tstart, tend, 'Orbit', true);
			times        = reshape( MrCDF_Epoch_Encode(tt2000_orbit), size(tt2000_orbit) );
		end

	%------------------------------------%
	% Process Data                       %
	%------------------------------------%
		for kk = 1 : size(times, 2)
			% Keep track of processing time
			f0 = now();
			oLog.AddText( ['Processing started at ' datestr(f0, 'yyyy-mm-ddTHH:MM:SS')] );
		
			% Attempt to create the files
			try
				if strcmp(mode{jj}, 'brst')
					[fstat, fgm_file, scm_file] = mms_fsm_bkgd_sdc(sc{ii}, mode{jj}, times{kk}, varargin{:});
				else
					[fstat, fgm_file, scm_file] = mms_fsm_bkgd_sdc(sc{ii}, mode{jj}, times{1,kk}, times{2,kk}, varargin{:});
				end
			
			% Capture any errors
			catch ME
				fstat    = 100;
				fstart   = MrTimeParser(times{1,kk}, '%Y-%M-%dT%H:%m:%S', '%Y%M%d%H%m%S');
				fgm_file = [sc{ii} '_' instr '_' mode{ii} '_' level '_cal-' fgm_instr '_' fstart];
				scm_file = [sc{ii} '_' instr '_' mode{ii} '_' level '_cal-scm'        '_' fstart];
				oLog.AddError(ME);
			end

		%------------------------------------%
		% Results                            %
		%------------------------------------%
			% Save results
			f1              = now();
			count           = count + 1;
			files(:, count) = {fgm_file, scm_file};
			status(count)   = fstat;
			telapsed(count) = (f1 - f0) * 86400.0;
			dt              = sprintf( '%dm %0.2fs', floor( telapsed(count) / 60.0 ), mod( telapsed(count), 60.0 ) );
			
			% Record to log file
			oLog.AddText( ['Processing finished at ' datestr(now(), 'yyyy-mm-ddTHH:MM:SS')] );
			oLog.AddText( ['    FGM Output: ' fgm_file] );
			oLog.AddText( ['    SCM Output: ' scm_file] );
			oLog.AddText( ['    Status:     ' num2str(fstat)] );
			oLog.AddText( ['    Proc Time:  ' dt] );
			oLog.AddText( '-----------------------------------------------------------' );
		end
	end
	end

%------------------------------------%
% Executive Summary                  %
%------------------------------------%
	% Trim results
	files  = files(:, 1:count);
	status = status(1:count);
	
	% Output information
	nWarn    = sum( find( status > 0 & status <  100 ) );
	nErr     = sum( find( status >= 100 ) );
	dt       = (f1 - t0) * 86400.0;
	dt_total = sprintf( '%dh %dm %0.2fs', floor(dt/3600.0), floor(mod(dt, 3600.0) / 60.0), mod(dt, 60.0) );

	% Summary
	oLog.AddText( '' );
	oLog.AddText( '////////////////////////////////////////////////////' );
	oLog.AddText( '////////////////////////////////////////////////////' );
	oLog.AddText( 'EXECUTIVE SUMMARY:' );
	oLog.AddText( ['    Number of Files:    ' num2str(count) ] );
	oLog.AddText( ['    Number of Warnings: ' num2str(nWarn) ] );
	oLog.AddText( ['    Number of Errors:   ' num2str(nErr)  ] );
	oLog.AddText( ['    Total Time Elapsed: ' dt_total  ] );
	oLog.AddText(  '    Status   Name' );
	
	% Add each file
	for ii = 1 : count
		oLog.AddText(  sprintf('    %3i      %s  %s', status(ii), files{:,ii}) );
	end
	
	oLog.AddText( '////////////////////////////////////////////////////' );
	oLog.AddText( '////////////////////////////////////////////////////' );
end