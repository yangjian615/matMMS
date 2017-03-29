%
% Name
%   mms_fsm_l3_sdc
%
% Purpose
%   Merge MMS fluxgate and search coil magnetometer data in the frequency
%   domain.
%
% Calling Sequence
%   [B_MERGED, T_MERGED] = mms_fsm_l3_sdc(SC. TSTART, TEND);
%     Given the MMS spacecraft number SC (e.g. 'mms1'), and a time
%     interval, [TSTART, TEND), gather all of the required search coil and
%     fluxgate magnetometer data and merge them into a single dataset
%     B_MERGED with time stamps T_MERGED.
%
%   [..., B_FG, T_FG] = mms_fsm_l3_sdc(__);
%     Also return the calibrated FGM magnetic field B_FG and its time
%     stamps T_FG.
%
%   [..., B_SC, T_SC] = mms_fsm_l3_sdc(__);
%     Also return the *UN*calibrated SCM magnetic field B_SC and its time
%     stamps T_SC.
%
% Parameters
%   SC:             in, required, type=char
%   TSTART:         in, required, type=char
%   TEND:           in, required, type=char
%   'Duration':     in, optional, type=double, default=20.0
%                   The duration of each merging interval. Sets the
%                     frequency resolution of the final dataset.
%   'FMax':         in, optional, type=double, default=Nyquist frequency
%                   The maximum of the frequency range to merge.
%   'FMin':         in, optional, type=double, default=df
%                   The minimum ( > 0 ) of the frequency range to merge.
%   'FHeavySide':   in, optional, type=double, default=4.0
%                   In the absence of a filter function for joining SCM and FGm data,
%                       the default will be a heavy-side filter with this cut-off
%                       frequency. FGM will be used below and SCM will be use above
%                       this frequency.
%   'FGMCalDir':    in, optional, type=char, default=''
%                   Directory in which to find FGM calibration data.
%   'NoLog':        in, optional, type=logical, default=false;
%                   If set, status information will be written to stderr instead of
%                       a log file.
%   'OptDesc'       in, optional, type=char, default=''
%                   Optional descriptor for the output file name.
%   'SCMCalDir':    in, optional, type=char, default='/home/argall/data/lpp_scm/cal'
%                   Directory in which to find SCM calibration data.
%   'FGMModelDir':  in, optional, type=char, default='/home/argall/MATLAB/fischer/'
%                   Directory in which to find SCM calibration data.
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-12-08      Written by Matthew Argall
%
function status = mms_fsm_l3_sdc(sc, mode, tstart, tend, varargin)

	% Global path variables
	global cal_path_root data_path_root dropbox_root hk_root log_path_root unh_data_root

	% Establish data paths
	mms_fsm_init();
	t0     = now();
	status = 0;
	
	% Create a temporary error logging object
	oLog = MrLogFile('stderr');

%------------------------------------%
% Defaults                           %
%------------------------------------%

	% Srvy or Brst?
	%   - In Brst, TEND will be an optional parameter.
	if nargin > 3 && strcmp(mode, 'brst')
		varargin = [ tend, varargin ];
	end
	
	% Defaults
	duration      = 2.0;
	f_heavyside   = 4.0;
	fgm_instr     = 'dfg';
	fgm_cal_dir   = '';
	outoptdesc    = '';
	scm_cal_dir   = '/home/argall/data/lpp_scm/cal';
	fgm_model_dir = '/home/argall/MATLAB/fischer/';
	tf_log        = true;

	% Optional parameters
	nOptArgs = length(varargin);
	for ii = 1 : 2 : nOptArgs
		switch varargin{ii}
			case 'Duration'
				duration = varargin{ii+1};
			case 'FGMInstr'
				fgm_instr = varargin{ii+1};
			case 'FGMCalDir'
				fgm_cal_dir = varargin{ii+1};
			case 'FHeavySide'
				f_heavyside = varargin{ii+1};
			case 'NoLog'
				tf_log = ~varargin{ii+1};
			case 'OptDesc'
				outoptdesc = varargin{ii+1};
			case 'SCMCalDir'
				scm_cal_dir = varargin{ii+1};
			case 'FGMModelDir'
				fgm_model_dir = varargin{ii+1};
			otherwise
				error('MMS:FSM_L2Plus:SDC', 'Optional argument not recognized: "%s"', varargin{ii})
		end
	end

%------------------------------------%
% Restrictions                       %
%------------------------------------%

	% Check inputs types
	if ~ischar(sc)
		status = 100;
		oLog.AddError('MMS:FSM_L3_SDC:BadInput', 'SC must be a char array.', status);
		return
	end
	if ~ischar(mode)
		status = 100;
		oLog.AddError('MMS:FSM_L3_SDC:BadInput', 'MODE must be a char array.', status);
		return
	end
	if ~ischar(tstart)
		status = 100;
		oLog.AddError('MMS:FSM_L3_SDC:BadInput', 'TSTART must be a char array.', status);
		return
	end
	
	% SC
	if ~ismember(sc, {'mms1', 'mms2', 'mms3', 'mms4'})
		status = 100;
		oLog.AddError('MMS:FSM_L3_SDC:BadInput', 'Invalid spacecraft identifier: "%s".', sc)
		return
	end
	
	% MODE
	if ~ismember(mode, {'srvy', 'brst'})
		status = 100;
		oLog.AddError('MMS:FSM_L3_SDC:BadInput', sprintf('Invalid telemetry mode: "%s".', mode), status)
		return
	end
	
	% TSTART
	if ~MrTokens_IsMatch(tstart, '%Y%M%d%H%m%S')
		status = 100;
		oLog.AddError('MMS:FSM_L3_SDC:BadInput', 'TSTART must be formatted as YYYYMMDDhhmmss', status);
		return
	end
	
	% TEND
	if ~strcmp(mode, 'brst') && ~MrTokens_IsMatch(tend, '%Y%M%d%H%m%S')
		status = 100;
		oLog.AddError('MMS:FSM_L3_SDC:BadInput', 'TEND must be formatted as YYYYMMDDhhmmss', status);
		return
	end
	
	% Constants for output file
	outinstr   = 'fsm';
	outlevel   = 'l3';
	outmode    = mode;

%------------------------------------%
% Create Log File                    %
%------------------------------------%
	%
	% Set "stdlog" so that all messages sent via mrfprintf lands
	% in the log file.
	%
	
	if tf_log
		% Current time
		ahora = datestr(now(), 'yyyymmdd_HHMMSS');

		% Log file name
		logFile = [sc '_' outinstr '_' mode '_' outlevel '_' tstart '_' ahora '.log'];
		
		% Build log directory
		logDir = mms_create_path(log_path_root, sc, outinstr, mode, outlevel, tstart, outoptdesc);
		if exist(logDir, 'dir') ~= 7
			mkdir( logDir );
		end
		
		% Set the log file
		logFile = fullfile(logDir, logFile);
	else
		logFile = 'stderr';
	end
	
	% Create the log file object
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
		t_orbit = cell(1,2);
	else
		iso_end      = MrTimeParser(tend, '%Y%M%d%H%m%S', '%Y-%M-%dT%H:%m:%S');
		tt2000_orbit = mms_bss_roi_get(iso_end, '', 'Orbit', true);
		t_orbit      = MrCDF_Epoch_Encode(tt2000_orbit);
	end
	
%------------------------------------%
% Find Files (BRST)                  %
%------------------------------------%
	
	% SCM Optional Descriptor
	switch mode
		case 'brst'
			optdesc_scm = 'scb';
		case 'srvy'
			optdesc_scm = 'scsrvy';
		case 'fast'
			optdesc_scm = 'scf';
		case 'slow'
			optdesc_scm = 'scs';
		otherwise
			error( ['Invalid mode: "' mode '".'] );
	end
	
	% BRST
	if strcmp(mode, 'brst')
		% FGM L1B
		f_l1a_fgm   = mms_latest_file( dropbox_root, sc, fgm_instr, mode, 'l1a', tstart, ...
		                               'RootDir', data_path_root);
		
		% SCM L1B
		f_l1b_scm   = mms_latest_file( dropbox_root, sc, 'scm', mode, 'l1b', tstart, ...
		                               'OptDesc', optdesc_scm, ...
		                               'RootDir', data_path_root);

		% Make sure all files are found
		if isempty(f_l1a_fgm)
			status = 101;
			oLog.AddError( 'MMS:FSM_L3_SDC:NoFileFound', ['No ' fgm_instr ' brst l1a files found.'], status);
			return
		end
		if isempty(f_l1b_scm)
			status = 101;
			oLog.AddError( 'MMS:FSM_L3_SDC:NoFileFound', 'No scm brst l1b files found.', status );
			return
		end

%------------------------------------%
% Find SRVY Files                    %
%------------------------------------%
	else
		% SCM L1B
		%   - After Sept. ##, SLOW and FAST are the same.
		f_l1b_scm = mms_find_file( sc, 'scm', mode, 'l1b', ...
		                           'OptDesc', optdesc_scm, ...
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
		if isempty(f_l1a_fgm)
			status = 101;
			oLog.AddError( 'MMS:FSM_L3_SDC:NoFileFound', ['No ' fgm_instr ' fast l1a files found.'], status );
			return
		end
		if isempty(l1a_slow_fgm)
			status = 101;
			oLog.AddError( 'MMS:FSM_L3_SDC:NoFileFound', ['No ' fgm_instr ' slow l1a files found.'], status );
			return
		end
		if isempty(f_l1b_scm)
			status = 101;
			oLog.AddError( 'MMS:FSM_L3_SDC:NoFileFound', 'No scm l1b files found.', status );
			return
		end
	end

%------------------------------------%
% DSS: Find File                     %
%------------------------------------%

	% DSS files are formatted as YYYYMMDD, even in burst mode
	dss_file = mms_latest_file(hk_root, sc, 'fields', 'hk', 'l1b', tstart(1:8), ...
	                           'OptDesc', '101', ...
	                           'RootDir', hk_root);

%------------------------------------%
% FDOA: Find File                    %
%------------------------------------%
	defatt_file = mms_anc_search(dropbox_root, sc, 'defatt', tstart, 'RootDir', data_path_root);
	
	if isempty(defatt_file)
		status = 101;
		oLog.AddError( 'MMS:FSM_L3_SDC:NoFileFound', 'No DEFATT file found.', status );
		return
	end

%------------------------------------%
% Parent Files                       %
%------------------------------------%
	% Write parents to log file
	mrfprintf('logtext', '')
	mrfprintf('logtext', '---------------------------------')
	mrfprintf('logtext', '| Parent Files                  |')
	mrfprintf('logtext', '---------------------------------')
	mrfprintf('logtext', f_l1a_fgm)
	mrfprintf('logtext', f_l1b_scm)
	mrfprintf('logtext', '---------------------------------')
	mrfprintf('logtext', '')

%------------------------------------%
% Read Attitude and zMPA             %
%------------------------------------%
	[defatt, att_hdr] = mms_fdoa_read_defatt(defatt_file);
	zmpa              = att_hdr.zMPA(:,1);

%------------------------------------%
% Read Data                          %
%------------------------------------%

	% FGM & SCM
	fgm = mms_fsm_fgm_read( f_l1a_fgm, t_orbit );
	scm = mms_fsm_scm_read( f_l1b_scm, t_orbit );
	
	% Add the model directory
	fgm.model_dir = fgm_model_dir;

%------------------------------------%
% Process Data                       %
%------------------------------------%
	[t_fsm, b_omb] = mms_fsm_l3_create( fgm, scm, duration, []);
	
%------------------------------------%
% Write to File                      %
%------------------------------------%
	% Parent files
	parents      = { f_l1a_fgm f_l1b_scm };
	[~, parents] = cellfun(@fileparts, parents, 'UniformOutput', false);

	% Create a data structure
	data = struct( 'tt2000', t_fsm',  ...
	               'b',      single( sqrt( sum( b_omb.^2 ) )' ), ...
	               'b_omb',  single(b_omb)' ...
	             );

	% Write the file
	fsm_file = mms_fsm_l3_write( sc, mode, tstart, data, ...
	                             'OptDesc', outoptdesc,  ...
	                             'Parents', parents );

%------------------------------------%
% Rotate to GSE                      %
%------------------------------------%
	
	%
	% Create FGM L2 data in OMB
	%
	cmd = [ fileparts(mfilename('fullpath')) filesep() 'mms_fsm_l3_rotate.sh ' fsm_file];
	if tf_log
		cmd = [cmd ' ' logFile];
	end
	
	mrfprintf( 'logtext',  '\n\n' );
	mrfprintf( 'logtext',  '=====================================' );
	mrfprintf( 'logtext',  '| Calling IDL to Rotate Data        |' );
	mrfprintf( 'logtext',  '=====================================' );
	mrfprintf( 'logtext', '%s', cmd );
	mrfprintf( 'logtext', '\n\n' );
	
	% Call IDL
	[status, cmdout] = system(cmd);
	
	% system will capture output from both stdout and stderr
	%   - If no log file is defined, redirect cmdout to stdout
	if ~tf_log
		mrfprintf( 'stdout', '%s', cmdout );
	end
	
	mrfprintf( 'logtext',  '\n\n' );
	mrfprintf( 'logtext', '=====================================' );
	mrfprintf( 'logtext', '| Returning to MATLAB               |' );
	mrfprintf( 'logtext', '=====================================' );
	mrfprintf( 'logtext',  '\n\n' );
	
	% Check status
	if status >= 100
		mrfprintf( 'logerr', 'Error rotating L3 data.' );
		return
	end

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
	mrfprintf('logtext', fsm_file);
	mrfprintf('logtext', 'Processing time: %s', dt_text);
	mrfprintf('logtext', '');

%------------------------------------%
% Clean Up                           %
%------------------------------------%

	% Close the log file by returning to stderr
	oLog = mrstdlog('stderr');
end