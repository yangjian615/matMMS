%
% Name
%   mms_fsm_bkgd_consolidate
%
% Purpose
%   Consolidate noise floor statistics from individual orbits (for survey data) or
%   burst intervals. Histogram statistics are summed, fit to a bi-gaussian distribution,
%   then the results are written to a CDF file.
%
% Calling Sequence
%   [STATUS, FOUT] = mms_fsm_bkgd_sdc(SC, MODE, OPTDESC, TSTART, TEND)
%     Read calibration data for spacecraft SC, data telemetry mode MODE, and
%     optional descriptor OPTDESC between the time of TSTART and TEND. The
%     times must be formatted as yyyymmdd. Processing success is reported
%     in STATUS and the resulting CDF file name returned in FOUT.
%
%   [__] = mms_fsm_bkgd_sdc(..., 'ParamName', ParamValue)
%     Any parameter name-value pairs given below may be used with any of
%     the above calling sequences.
%
% Parameters
%   SC              in, required, type = char
%   MODE            in, required, type = char
%   OPTDESC         in, required, type = char
%   TSTART          in, required, type = char
%   TEND            in, optional, type = char
%   'Interval'      in, optional, type = char, default = 'week'
%                   The time period over which to accumulate statistics. Options are:
%                     { 'day' | 'week' | 'month' }
%   'NoLog'         in, optional, type = logical, default = false
%                   If true, no log file will be created. If false, all messages will
%                     be directed to the log file.
%
% Returns
%   STATUS          out, required, type=integer
%   FOUT            out, required, type=string
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2016-09-23      Written by Matthew Argall
%
function [status, fout] = mms_fsm_bkgd_consolidate(sc, mode, optdesc, tstart, tend, varargin)

	% Declare global attributes
	global cal_path_root data_path_root dropbox_root hk_root log_path_root unh_data_root

	% Initialize
	mms_fsm_init();
	t0 = now();

%------------------------------------%
% Check Inputs                       %
%------------------------------------%
	% Defaults
	interval = 'week';
	tf_log   = true;
	
	% Optional parameters
	nOptArgs = length(varargin);
	for ii = 1 : 2 : nOptArgs
		switch varargin{ii}
			case 'Interval'
				interval = varargin{ii+1};
			case 'NoLog'
				tf_log = ~varargin{ii+1};
			otherwise
				error( ['Unrecognized optional parameter: "' varargin{ii} '".'] );
		end
	end

	% Restrictions
	assert( all( ismember(sc,       {'mms1', 'mms2', 'mms3', 'mms4'}) ),  'SC must be "mms1", "mms2", "mms3", or "mms4".' );
	assert( all( ismember(optdesc,  {'cal-afg', 'cal-dfg', 'cal-scm'}) ), 'OPTDESC must be {"cal-dfg" | "cal-afg" | "cal-fgm"}.');
	assert( ismember(mode,          {'srvy', 'brst'}),                    'MODE must be "srvy" or "brst".' );
	assert( ismember(interval,      {'day', 'week', 'month'}),            'INTERVAL must be {"day" | "week" | "month"}.');
	assert( MrTokens_IsMatch(tstart, '%Y%M%d'),                           'TSTART must be formatted as YYYYMMDD');
	assert( MrTokens_IsMatch(tend,   '%Y%M%d'),                           'TEND must be formatted as YYYYMMDD');
	assert( ~( strcmp(mode, 'srvy') && strcmp(interval, 'day') ),         'Cannot make daily plots of survey data.');

	% Constants
	instr   = 'fsm';
	level   = 'l2plus';
	outdesc = [optdesc '-' interval];

%------------------------------------%
% Create Log File                    %
%------------------------------------%
	%
	% Set "stdlog" so that all messages sent via mrfprintf land
	% in the log file.
	%

	if tf_log
		% Current time
		ahora = datestr(now(), 'yyyymmdd_HHMMSS');

		% Log file name
		logFile = [sc '_' instr '_' mode '_' level '_' outdesc '_' ahora '.log'];
		
		% Build log directory
		logDir = mms_create_path(log_path_root, sc, instr, mode, level, tstart, optdesc);
		if exist(logDir, 'dir') ~= 7
			mkdir( logDir );
		end
		
		% Set the log file
		logFile = fullfile(logDir, logFile);
	else
		logFile = 'stderr';
	end
	
	% Get the log file object
	oLog = mrstdlog(logFile);

%------------------------------------%
% Outline Data Interval              %
%------------------------------------%
	
	% Convert input time to TT2000
	[~, tstart_tt2000] = mms_parse_time(tstart);
	[~, tend_tt2000]   = mms_parse_time(tend);

%------------------------------------%
% Loop Through Products              %
%------------------------------------%
	% File start time
	fstart_tt2000 = tstart_tt2000;
	
	% Number of nanoseconds between files
	switch interval
		case 'day'
			tvec        = MrCDF_Epoch_Breakdown(fstart_tt2000);
			tvec(3)     = tvec(3) + 1;
			fend_tt2000 = MrCDF_Epoch_Compute( tvec, 'CDF_TIME_TT2000' );
		case 'week'
			tvec        = MrCDF_Epoch_Breakdown(fstart_tt2000);
			tvec(3)     = tvec(3) + 7;
			fend_tt2000 = MrCDF_Epoch_Compute( tvec, 'CDF_TIME_TT2000' );
		case 'month'
			tvec        = MrCDF_Epoch_Breakdown(fstart_tt2000);
			tvec(2)     = tvec(2) + 1;
			fend_tt2000 = MrCDF_Epoch_Compute( tvec, 'CDF_TIME_TT2000' );
		otherwise
			error( 'Invalid value for INTERVAL.')
	end
	
	% Number of files processed
	nalloc = 100;
	count  = 0;
	fout   = cell(1,nalloc);
	status = zeros(1,nalloc, 'uint8');

	% Step through spacecraft, instr, mode
	for ii = 1 : length(sc)
	for jj = 1 : length(optdesc)
	while fstart_tt2000 < tend_tt2000

	%------------------------------------%
	% Read Files                         %
	%------------------------------------%
		% Convert tt2000 times to strings
		fstart = MrCDF_Epoch_Encode(fstart_tt2000);
		fend   = MrCDF_Epoch_Encode(fend_tt2000);
		
		% Find files between start and end times
		fsm_files = mms_find_file( sc, instr, mode, level,      ...
		                           'Dropbox',   dropbox_root,   ...
		                           'OptDesc',   optdesc,        ...
		                           'SDCroot',   data_path_root, ...
		                           'TimeOrder', '%Y%M%d%H%m%S', ...
		                           'TStart',    fstart,         ...
		                           'TEnd',      fend );
		
		% Read the data
		data = mms_fsm_bkgd_read(fsm_files);

	%------------------------------------%
	% Fit Data                           %
	%------------------------------------%
		
		% Fit the distribution with a gaussian
		data.('psd_floor') = zeros( length(data.comp), length(data.flag), length(data.f), 'single' );
%		data.('psd_floor') = mms_fsm_bkgd_fit_bigauss( data.psd_hist, data.f, data.psd_bins );

	%------------------------------------%
	% Output to File                     %
	%------------------------------------%
		% Show parent files
		oLog.AddText('');
		oLog.AddText('-----------------------------------');
		oLog.AddText('| PARENT FILES                    |');
		oLog.AddText('-----------------------------------');
		oLog.AddText( fsm_files );
		oLog.AddText('');
		
		% Parents
		[~, names, ext] = cellfun( @fileparts, fsm_files, 'UniformOutput', false );
		parents         = strcat( names, ext );
		
		% Start time of file
		fstart = MrTimeParser( fstart, '%Y-%M-%dT%H:%m:%S%f', '%Y%M%d%H%m%S' );
		
		% Write data to file
		fout{count+1} = mms_fsm_bkgd_consol_write( sc, mode, outdesc, fstart, data, ...
		                                          'Parents', parents );

	%------------------------------------%
	% Record Output                      %
	%------------------------------------%
		% Processing time
		dt      = (now() - t0) * 86400.0;
		dt_text = sprintf('%dh %dm %0.2fs', floor(dt/3600.0), floor(mod(dt, 3600.0) / 60.0), mod(dt, 60));
		count   = count + 1;
	
		% Results
		mrfprintf('logtext', '');
		mrfprintf('logtext', '-----------------------------------');
		mrfprintf('logtext', '| RESULTS                         |');
		mrfprintf('logtext', '-----------------------------------');
		mrfprintf('logtext', fout{count} );
		mrfprintf('logtext', 'Processing time: %s', dt_text);
		mrfprintf('logtext', '');

	%------------------------------------%
	% Next Interval                      %
	%------------------------------------%
		fstart_tt2000 = fend_tt2000;
	
		% Number of nanoseconds between files
		switch interval
			case 'day'
				tvec        = MrCDF_Epoch_Breakdown(fstart_tt2000);
				tvec(3)     = tvec(3) + 1;
				fend_tt2000 = MrCDF_Epoch_Compute( tvec, 'CDF_TIME_TT2000' );
			case 'week'
				tvec        = MrCDF_Epoch_Breakdown(fstart_tt2000);
				tvec(3)     = tvec(3) + 7;
				fend_tt2000 = MrCDF_Epoch_Compute( tvec, 'CDF_TIME_TT2000' );
			case 'month'
				tvec        = MrCDF_Epoch_Breakdown(fstart_tt2000);
				tvec(2)     = tvec(2) + 1;
				fend_tt2000 = MrCDF_Epoch_Compute( tvec, 'CDF_TIME_TT2000' );
			otherwise
				error( 'Invalid value for INTERVAL.')
		end
	end
	end
	end

%------------------------------------%
% Executive Summary                  %
%------------------------------------%
	% Trim results
	fout   = fout(1:count);
	status = status(1:count);
	
	% Output information
	nWarn    = sum( find( status > 0 & status <  100 ) );
	nErr     = sum( find( status >= 100 ) );
	dt       = (now() - t0) * 86400.0;
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
		oLog.AddText(  sprintf('    %3i      %s  %s', status(ii), fout{ii}) );
	end
	
	oLog.AddText( '////////////////////////////////////////////////////' );
	oLog.AddText( '////////////////////////////////////////////////////' );

%------------------------------------%
% Clean Up                           %
%------------------------------------%

	% Close the log file by returning to stderr
	oLog = mrstdlog('stderr');
end