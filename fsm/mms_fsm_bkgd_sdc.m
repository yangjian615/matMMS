%
% Name
%   mms_fsm_bkgd_sdc
%
% Purpose
%   Generate calibration data for the fluxgate-searchcoil merged (FSM) data
%   product.
%
% Calling Sequence
%   [FILE_FGM, FILE_SCM] = mms_fsm_bkgd_sdc(SC, MODE, TSTART)
%     Read burst mode FGM and SCM magnetometer data, compute noise floor
%     parameters, and write the results to FILE_FGM and FILE_SCM. TSTART
%     must be formatted as yyyy-mm-ddTHH:MM:SS
%
%   [FILE_FGM, FILE_SCM] = mms_fsm_bkgd_sdc(SC, MODE, TSTART, TEND)
%     Read slow, fast, or srvy mode FGM and SCM magnetometer data,
%     compute noise floor parameters, and write the results to 
%     FILE_FGM and FILE_SCM. TSTART and TEND are intended to outline
%     an orbit of data, starting with slow survy, and ending after
%     fast survey.
%
%   FG_QL = mms_fg_read_ql(__, 'ParamName', ParamValue)
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
%   'fgmInstr'      in, optional, type = double, default = 'dfg'
%                   FGM instrument to use. Options are {'afg' | 'dfg'}
%   'noLog'         in, optional, type = boolean, default = false
%                   If true, error messages are printed to the console instead of
%                       a log file.
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
function [file_fgm, file_scm] = mms_fsm_bkgd_sdc(sc, mode, tstart, tend, varargin)

	% Declare global attributes
	global cal_path_root data_path_root dropbox_root hk_root log_path_root unh_data_root

	% Initialize
	mms_fsm_init();
	t0 = now();

%------------------------------------%
% Check Inputs                       %
%------------------------------------%

	% Srvy or Brst?
	%   - In Brst, TEND will be an optional parameter.
	if nargin > 3 && strcmp(mode, 'brst')
		varargin = [ tend, varargin ];
	end

	% Defaults
	T         = 20.0;
	fgm_instr = 'dfg';
	tf_log    = true;
	
	% Optional parameters
	nOptArgs = length(varargin);
	for ii = 1 : 2 : nOptArgs
		switch varargin{ii}
			case 'Duration'
				T = varargin{ii+1};
			case 'fgmInstr'
				fgm_instr = varargin{ii+1};
			case 'noLog'
				tf_log = ~varargin{ii+1};
			otherwise
				error( ['Unrecognized optional parameter: "' varargin{ii} '".'] );
		end
	end

	% Restrictions
	assert( ismember(sc, {'mms1', 'mms2', 'mms3', 'mms4'}), 'SC must be "mms1", "mms2", "mms3", or "mms4".' )
	assert( ismember(mode, {'srvy', 'brst'}), 'MODE must be "srvy" or "brst".' )

	% Constants
	instr   = 'fsm';
	level   = 'l2plus';
	optdesc = [fgm_instr '-scm'];
	
	% Convert input times to file times
	fstart = MrTimeParser(tstart, '%Y-%M-%dT%H:%m:%S', '%Y%M%d%H%m%S');
	if ~strcmp(mode, 'brst')
		fend   = MrTimeParser(tend,   '%Y-%M-%dT%H:%m:%S', '%Y%M%d%H%m%S');
	end

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
		logFile = [sc '_' instr '_' mode '_' level '_' optdesc '_' fstart '_' ahora '.log'];
		
		% Build log directory
		logDir = mms_create_path(log_path_root, sc, instr, mode, level, fstart, optdesc);
		if exist(logDir, 'dir') ~= 7
			mkdir( logDir );
		end
		
		% Set the log file
		logFile = fullfile(logDir, logFile);
	else
		logFile = 'stdlog';
	end
	
	% Get the log file object
	oLog = mrstdlog(logFile);

%------------------------------------%
% Find time of ROI                   %
%------------------------------------%
	%
	% TODO: In phase 2, the orbit will be more than one day long!
	%       Calibration intervals, therefore, will also be multiple
	%       days. How to increment the orbit?
	%

	% Brst: Read the entire file
	% Srvy: Read all data within an orbit [slow, fast]
	%   - Convert tt2000 to date-time string
	%   - yyyy-mm-ddThh:mm:ss.fff[...]
	if strcmp(mode, 'brst')
		t_orbit = {};
	else
		tt2000_orbit = mms_bss_roi_get(tend, '', 'Orbit', true);
		t_orbit      = MrCDF_Epoch_Encode(tt2000_orbit);
	end
	
%------------------------------------%
% Find Files (BRST)                  %
%------------------------------------%
	
	% SCM Optional Descriptor
	switch mode
		case 'srvy'
			scm_optdesc = 'scsrvy';
		case 'brst'
			scm_optdesc = 'scb';
		otherwise
			error(['Invalid mode: "' mode '".']);
	end
	
	% BRST
	if strcmp(mode, 'brst')
		% FGM L2PRE
		f_l2pre_fgm = mms_latest_file( dropbox_root, sc, fgm_instr, mode, 'l2pre', fstart, ...
		                               'RootDir', data_path_root);
		
		% FGM L1B
		f_l1a_fgm   = mms_latest_file( dropbox_root, sc, fgm_instr, mode, 'l1a', fstart, ...
		                               'RootDir', data_path_root);
		
		% SCM L1B
		f_l1b_scm   = mms_latest_file( dropbox_root, sc, 'scm', mode, 'l1b', fstart, ...
		                               'OptDesc', scm_optdesc, ...
		                               'RootDir', data_path_root);
		
		% Make sure all files are found
		assert(~isempty(f_l2pre_fgm), ['No ' fgm_instr ' brst l2pre files found.']);
		assert(~isempty(f_l1a_fgm),   ['No ' fgm_instr ' brst l1a files found.']);
		assert(~isempty(f_l1b_scm),   ['No scm brst l1b files found.']);

%------------------------------------%
% Find SRVY Files                    %
%------------------------------------%
	else
		% FGM L2Pre
		f_l2pre_fgm = mms_find_file( sc, 'dfg', mode, 'l2pre', ...
		                             'TStart', t_orbit{1},     ...
		                             'TEnd',   t_orbit{2} );
	
		% SCM L1B
		%   - After Sept. ##, SLOW and FAST are the same.
		f_l1b_scm = mms_find_file( sc, 'scm', mode, 'l1b', ...
		                           'OptDesc', scm_optdesc, ...
		                           'TStart',  t_orbit{1},  ...
		                           'TEnd',    t_orbit{2} );

		% FGM L1A has 'brst', 'slow', 'fast' (not 'srvy')
		if ismember(mode, {'fast', 'srvy'})
			f_l1a_fgm = mms_find_file( sc, 'dfg', 'fast', 'l1a', ...
			                           'TStart', t_orbit{1},     ...
			                           'TEnd',   t_orbit{2} );
		end
		if ismember(mode, {'slow', 'srvy'})
			l1a_slow_fgm = mms_find_file( sc, 'dfg', 'slow', 'l1a', ...
			                             'TStart', t_orbit{1},      ...
			                             'TEnd',   t_orbit{2} );
			if isempty(f_l1a_fgm)
				f_l1a_fgm = l1a_slow_fgm;
			else
				f_l1a_fgm = [f_l1a_fgm, l1a_slow_fgm];
			end
		end
		
		% Make sure all files were found
		%   TODO: It might be better to use the second output COUNT to
		%         ensure the correct number of files.
		assert(~isempty(f_l2pre_fgm),  ['No ' fgm_instr ' srvy l2pre files found.']);
		assert(~isempty(f_l1a_fgm),    ['No ' fgm_instr ' fast l1a files found.']);
		assert(~isempty(l1a_slow_fgm), ['No ' fgm_instr ' slow l1a files found.']);
		assert(~isempty(f_l1b_scm),    ['No scm l1b files found.']);
	end

%------------------------------------%
% Parent Files                       %
%------------------------------------%
	oLog.AddText('');
	oLog.AddText('-----------------------------------');
	oLog.AddText('| PARENT FILES                    |');
	oLog.AddText('-----------------------------------');
	oLog.AddText(f_l1a_fgm);
	oLog.AddText(f_l2pre_fgm);
	oLog.AddText(f_l1b_scm);
	oLog.AddText('');

%------------------------------------%
% Read Data                          %
%------------------------------------%
	% Fast Survey
	fgm = mms_fsm_bkgd_fgm_read( f_l1a_fgm, f_l2pre_fgm, t_orbit );
	scm = mms_fsm_bkgd_scm_read( f_l1b_scm, t_orbit );

%------------------------------------%
% Compute Background                 %
%------------------------------------%
	fgm = mms_fsm_bkgd_compute_one(fgm, T);
	scm = mms_fsm_bkgd_compute_one(scm, T);

%------------------------------------%
% Write to File                      %
%------------------------------------%
	% Parent files
	if ischar(f_l1a_fgm)
		files = {f_l1a_fgm, f_l2pre_fgm, f_l1b_scm};
	else
		files = [f_l1a_fgm, f_l2pre_fgm, f_l1b_scm];
	end
	[~, fnames, ext] = cellfun(@fileparts, files, 'UniformOutput', false);
	parents          = strcat(fnames, ext);
	
	% Start time
	%   - Cal files are by orbit, not by day. Include hour, minute, seconds
	if ~strcmp(mode, 'brst')
		fstart = MrCDF_Epoch_Encode(fgm.t(1));
		fstart = MrTimeParser(fstart, '%Y-%M-%dT%H:%m:%S%f', '%Y%M%d%H%m%S');
	end

	% Output data
	file_fgm = mms_fsm_l2plus_cal_dfg_write(sc, mode, fstart, fgm, ...
	                                        'Parents', parents);
	file_scm = mms_fsm_l2plus_cal_scm_write(sc, mode, fstart, scm, ...
	                                        'Parents', parents);

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
	mrfprintf('logtext', file_fgm);
	mrfprintf('logtext', file_scm);
	mrfprintf('logtext', 'Processing time: %s', dt_text);
	mrfprintf('logtext', '');
	

%------------------------------------%
% Clean Up                           %
%------------------------------------%

	% Close the log file by returning to stderr
	oLog = mrstdlog('stderr');
end