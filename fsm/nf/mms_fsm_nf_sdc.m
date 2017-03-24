%
% Name
%   mms_fsm_nf_sdc
%
% Purpose
%   Read noise floor results from individual orbits (srvy) or intervals (brst)
%   and compute the mean noise floor.
%
% Calling Sequence
%   FNAMES = mms_fsm_nf_create( SC, MODE, INTERVAL, TSTART )
%     Process MMS FSM data for spacecraft SC in telemetry mode MODE throughout a
%     time period INTERVAL beginning at TSTART. INTERVAL can have values of
%     {'day' | 'week' | 'month'} and TSTART must be formatted as 'YYYY-MM-DD'.
%     Data files are created and their names are returned in FNAMES.
%
%   [__] = mms_fsm_nf_create( __, 'ParamName', ParamValue )
%     Any of the name-value pairs listed below may be given.
%
% Inputs
%   SC              in, required, type=char
%   MODE            in, required, type=char
%   INTERVAL        in, required, type=char
%   TSTART          in, required, type=char
%   'FGMInstr'      in, optional, type=char, default='dfg'
%                   FGM instrument for which to calculate the noise floor.
%   'NoLog'         in, optional, type=boolean, default=false
%                   If true, no log file will be created.
%
% Returns
%   FNAMES          out, required, type=string
%
% MATLAB release(s) MATLAB 9.0.0.341360 (R2016a)
% Required Products None
%
% History:
%   2016-10-20      Written by Matthew Argall
%
%***************************************************************************
function [file_nf_fgm, file_nf_scm] = mms_fsm_nf_sdc(sc, mode, interval, tstart, varargin)

	% Declare global attributes
	global cal_path_root data_path_root dropbox_root hk_root log_path_root unh_data_root
	
	% Initialize
	mms_fsm_init();
	t0 = now();

%------------------------------------%
% Check Inputs                       %
%------------------------------------%
	
	% Defaults
	fgm_instr = 'dfg';
	tf_log    = true;
	
	% Optional arguments
	nOptArgs = length(varargin);
	for ii = 1 : 2 : nOptArgs
		switch varargin{ii}
			case 'FGMInstr'
				fgm_instr = varargin{ii+1};
			case 'NoLog'
				tf_log = ~varargin{ii+1};
			otherwise
				error( ['Optional argument not recognized: "' varargin{ii} '".'] );
		end
	end
	
	% Source file constants
	instr       = 'fsm';
	level       = 'l2plus';
	optdesc_fgm = ['cal-' fgm_instr];
	optdesc_scm = 'cal-scm';
	
	% Destination file constants
	outdesc     = 'nf';
	outdesc_fgm = [outdesc '-' fgm_instr '-' interval];
	outdesc_scm = [outdesc '-' 'scm'     '-' interval];
	
	% TSTART must be a date
	assert( MrTokens_IsMatch( tstart, '%Y-%M-%d' ), 'TSTART must be formatted as YYYY-MM-DD.' );

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
		logFile = [sc '_' instr '_' mode '_' level '_' outdesc '_' fstart '_' ahora '.log'];
		
		% Build log directory
		logDir = mms_create_path(log_path_root, sc, instr, mode, level, fstart, optdesc);
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
% Time Range                         %
%------------------------------------%
	% Convert TSTART to TT2000
	[~, ts_tt2000] = mms_parse_time( MrTimeParser(tstart, '%Y-%M-%d', '%Y%M%d') );
	fstart         = [ tstart 'T00:00:00' ];
	
	% Time interval in nanoseconds
	switch interval
		case 'day'
			dt = int64( 86400.0 * 1e9 );
		case 'week'
			dt = int64( 7.0 * 86400.0 * 1e9 );
		case 'month'
			dt = int64( 30.0 * 86400.0 * 1e9 );
		otherwise
			error( ['Invalid INTERVAL value: "' interval '".'] )
	end
	
	% Compute end date
	te_tt2000 = ts_tt2000 + dt;
	te_str    = MrCDF_Epoch_Encode(te_tt2000);
	tend      = te_str{1}(1:19);
	
%------------------------------------%
% Find FGM Files                     %
%------------------------------------%
	fgm_files = mms_find_file( sc, instr, mode, level,      ...
	                           'Dropbox',   dropbox_root,   ...
	                           'OptDesc',   optdesc_fgm,    ...
	                           'SDCroot',   data_path_root, ...
	                           'TimeOrder', '%Y%M%d%H%m%S', ...
	                           'TStart' ,   fstart,         ...
	                           'TEnd',      tend );
	
%------------------------------------%
% Find SCM Files                     %
%------------------------------------%
	scm_files = mms_find_file( sc, instr, mode, level,      ...
	                           'Dropbox',   dropbox_root,   ...
	                           'OptDesc',   optdesc_scm,    ...
	                           'SDCroot',   data_path_root, ...
	                           'TimeOrder', '%Y%M%d%H%m%S', ...
	                           'TStart' ,   fstart,         ...
	                           'TEnd',      tend );
	
%------------------------------------%
% Find NF Files                      %
%------------------------------------%
	% FGM
	fgm_nf_file = mms_find_file( sc, 'fsm', mode, 'l2plus',       ...
	                             'Dropbox',       dropbox_root,   ...
	                             'OptDesc',       ['nf-' fgm_instr '-' interval], ...
	                             'SDCroot',       data_path_root, ...
	                             'TimeOrder',     '%Y%M%d%H%m%S', ...
	                             'TStart',        fstart,         ...
	                             'TEnd',          tend,           ...
	                             'RelaxedTStart', true );
	
	% SCM
	scm_nf_file = mms_find_file( sc, 'fsm', mode, 'l2plus',       ...
	                             'Dropbox',       dropbox_root,   ...
	                             'OptDesc',       ['nf-scm-' interval], ...
	                             'SDCroot',       data_path_root, ...
	                             'TimeOrder',     '%Y%M%d%H%m%S', ...
	                             'TStart',        fstart,         ...
	                             'TEnd',          tend,           ...
	                             'RelaxedTStart', true );

%------------------------------------%
% Parent Files                       %
%------------------------------------%
	oLog.AddText('');
	oLog.AddText('-----------------------------------');
	oLog.AddText('| PARENT FILES                    |');
	oLog.AddText('-----------------------------------');
	oLog.AddText(fgm_files);
	oLog.AddText(fgm_nf_file);
	oLog.AddText(scm_files);
	oLog.AddText(scm_nf_file);
	oLog.AddText('');

%------------------------------------%
% Process the Data                   %
%------------------------------------%
	% Determine the noise floor
	nf_new_fgm = mms_fsm_nf_create(fgm_files);
	nf_new_scm = mms_fsm_nf_create(scm_files);

	% Append times
	nf_new_fgm.t  = ts_tt2000;
	nf_new_fgm.dt = dt;
	nf_new_scm.t  = ts_tt2000;
	nf_new_scm.dt = dt;

	% Read old NF files
%	nf_old_fgm = mms_fsm_nf_read(fgm_nf_file);
%	nf_old_scm = mms_fsm_nf_read(scm_nf_file);

	% Combine old and new
%	nf_fgm = mms_fsm_nf_add(nf_old_fgm, nf_new_fgm);
%	nf_scm = mms_fsm_nf_add(nf_old_scm, nf_new_scm);
	nf_fgm = nf_new_fgm;
	nf_scm = nf_new_scm;

%------------------------------------%
% Write to File                      %
%------------------------------------%
	% FGM Parent files
	[~, fnames, ext] = cellfun(@fileparts, fgm_files, 'UniformOutput', false);
	parents_fgm      = strcat(fnames, ext);
	
	% SCM Parent files
	[~, fnames, ext] = cellfun(@fileparts, scm_files, 'UniformOutput', false);
	parents_scm      = strcat(fnames, ext);

	% Output data
	fstart      = MrTimeParser(fstart, '%Y-%M-%dT%H:%m:%S', '%Y%M%d%H%m%S');
	file_nf_fgm = mms_fsm_nf_write(sc, mode, outdesc_fgm, fstart, nf_fgm, 'Parents', parents_fgm);
	file_nf_scm = mms_fsm_nf_write(sc, mode, outdesc_scm, fstart, nf_scm, 'Parents', parents_scm);

%------------------------------------%
% Record Output                      %
%------------------------------------%
	% Processing time
	dt      = (now() - t0) * 86400.0;
	dt_text = sprintf('%dh %dm %0.2fs', floor(dt/3600.0), floor(mod(dt, 3600.0) / 60.0), mod(dt, 60));
	
	% Results
	mrfprintf('logtext', '');
	mrfprintf('logtext', '-----------------------------------');
	mrfprintf('logtext', '| RESULTS                         |');
	mrfprintf('logtext', '-----------------------------------');
	mrfprintf('logtext', file_nf_fgm);
	mrfprintf('logtext', file_nf_scm);
	mrfprintf('logtext', 'Processing time: %s', dt_text);
	mrfprintf('logtext', '');
	

%------------------------------------%
% Clean Up                           %
%------------------------------------%

	% Close the log file by returning to stderr
	oLog = mrstdlog('stderr');
end
