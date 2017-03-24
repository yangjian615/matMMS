%
% Name
%   mms_fsm_create_l1b
%
% Purpose
%   Merge MMS fluxgate and search coil magnetometer data in the frequency
%   domain.
%
% Calling Sequence
%   [B_MERGED, T_MERGED] = mms_fsm_merge(SC. TSTART, TEND);
%     Given the MMS spacecraft number SC (e.g. 'mms1'), and a time
%     interval, [TSTART, TEND), gather all of the required search coil and
%     fluxgate magnetometer data and merge them into a single dataset
%     B_MERGED with time stamps T_MERGED.
%
%   [..., B_FG, T_FG] = mms_fsm_merge(__);
%     Also return the calibrated FGM magnetic field B_FG and its time
%     stamps T_FG.
%
%   [..., B_SC, T_SC] = mms_fsm_merge(__);
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
function status = mms_fsm_l2plus_sdc(sc, mode, tstart, tend, varargin)
% fsm_file = mms_fsm_l2plus_sdc(fgm_file, scm_file, scm_cal_file, att_file, dss_file, duration)

	% Global path variables
	global cal_path_root data_path_root dropbox_root hk_root log_path_root unh_data_root

	% Establish data paths
	mms_fsm_init();
	t0     = now();
	status = 0;

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
		error('MMS:FSM_L2Plus:SDC', 'SC must be a char array.');
	end
	if ~ischar(mode)
		error('MMS:FSM_L2Plus:SDC', 'MODE must be a char array.');
	end
	if ~ischar(tstart)
		error('MMS:FSM_L2Plus:SDC', 'TSTART must be a char array.');
	end
	
	% Examine values
	if ~ismember({'mms1', 'mms2', 'mms3', 'mms4'}, sc)
		error('MMS:FSM_L2Plus:SDC', 'Invalid spacecraft identifier: "%s".', sc)
	end
	if ~ismember({'srvy', 'brst'}, mode)
		error('MMS:FSM_L2Plus:SDC', 'Invalid telemetry mode: "%s".', mode)
	end
	if ~ismember({'srvy', 'brst'}, mode)
		error('MMS:FSM_L2Plus:SDC', 'Invalid telemetry mode: "%s".', mode)
	end
	
	% Constants for output file
	outinstr   = 'fsm';
	outlevel   = 'l2plus';
	outmode    = mode;
	
	% Convert input times to file times
	assert( MrTokens_IsMatch(tstart, '%Y%M%d%H%m%S'), 'TSTART must be formatted as YYYYMMDDhhmmss');
	if ~strcmp(mode, 'brst')
		assert( MrTokens_IsMatch(tend, '%Y%M%d%H%m%S'), 'TEND must be formatted as YYYYMMDDhhmmss');
	end

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
		% FGM L2PRE
		f_l2pre_fgm = mms_latest_file( dropbox_root, sc, fgm_instr, mode, 'l2pre', tstart, ...
		                               'RootDir', data_path_root);
		
		% FGM L1B
		f_l1a_fgm   = mms_latest_file( dropbox_root, sc, fgm_instr, mode, 'l1a', tstart, ...
		                               'RootDir', data_path_root);
		
		% SCM L1B
		f_l1b_scm   = mms_latest_file( dropbox_root, sc, 'scm', mode, 'l1b', tstart, ...
		                               'OptDesc', optdesc_scm, ...
		                               'RootDir', data_path_root);
		
		% Make sure all files are found
		if isempty(f_l2pre_fgm)
			status = 101;
			error(['No ' fgm_instr ' brst l2pre files found.']);
		end
		if isempty(f_l1a_fgm)
			status = 101;
			error(['No ' fgm_instr ' brst l1a files found.'])
		end
		if isempty(f_l1b_scm)
			status = 101;
			error(['No scm brst l1b files found.'])
		end

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
		assert(~isempty(f_l2pre_fgm),  ['No ' fgm_instr ' srvy l2pre files found.']);
		assert(~isempty(f_l1a_fgm),    ['No ' fgm_instr ' fast l1a files found.']);
		assert(~isempty(l1a_slow_fgm), ['No ' fgm_instr ' slow l1a files found.']);
		assert(~isempty(f_l1b_scm),    ['No scm l1b files found.']);
	end

%------------------------------------%
% Find NF Files                      %
%------------------------------------%
	% FGM
	fgm_nf_file = mms_find_file( sc, 'fsm', mode, 'l2plus',                  ...
	                             'Dropbox',       dropbox_root,              ...
	                             'OptDesc',       ['nf-' fgm_instr '-week'], ...
	                             'SDCroot',       data_path_root,            ...
	                             'TimeOrder',     '%Y%M%d%H%m%S',            ...
	                             'TStart',        t_orbit{1},                ...
	                             'TEnd',          t_orbit{2},                ...
	                             'RelaxedTStart', true );
	
	% SCM
	scm_nf_file = mms_find_file( sc, 'fsm', mode, 'l2plus',       ...
	                             'Dropbox',       dropbox_root,   ...
	                             'OptDesc',       'nf-scm-week',  ...
	                             'SDCroot',       data_path_root, ...
	                             'TimeOrder',     '%Y%M%d%H%m%S', ...
	                             'TStart',        t_orbit{1},     ...
	                             'TEnd',          t_orbit{2},     ...
	                             'RelaxedTStart', true );

%------------------------------------%
% SCM: Find Cal Files                %
%------------------------------------%
	
	% Determine the flight model
	switch sc
		case 'mms1'
			fm = 'scm1';
		case 'mms2'
			fm = 'scm2';
		case 'mms3'
			fm = 'scm4';
		case 'mms4'
			fm = 'scm3';
		otherwise
			error(['Invalid spacecraft ID: "' sc '".'])
	end
	
	if strcmp(mode, 'brst')
		scm_tstart = tstart;
		scm_tend   = tstart;
	else
		scm_tstart = [tstart '000000'];
		scm_tend   = [tstart '240000'];
	end

	% SCM Cal File
	scm_ftest = fullfile(scm_cal_dir, [sc '_' fm '_caltab_%Y%M%d%H%M%S_v*.txt']);
	[scm_cal_file, count] = MrFile_Search(scm_ftest, ...
	                                      'Closest',      true, ...
	                                      'TimeOrder',    '%Y%M%d%H%M%S', ...
	                                      'TStart',       scm_tstart, ...
	                                      'TEnd',         scm_tend, ...
	                                      'TPattern',     '%Y%M%d%H%m%S', ...
	                                      'VersionRegex', 'v[0-9]');
	assert(count > 0, ['No SCM calibration file found: "' scm_ftest '".']);

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
	
	% No file found
	assert(~isempty(defatt_file), 'No DEFATT file found.')

%------------------------------------%
% Parent Files                       %
%------------------------------------%
	% Write parents to log file
	mrfprintf('logtext', '')
	mrfprintf('logtext', '---------------------------------')
	mrfprintf('logtext', '| Parent Files                  |')
	mrfprintf('logtext', '---------------------------------')
	mrfprintf('logtext', f_l1a_fgm)
	mrfprintf('logtext', f_l2pre_fgm)
	mrfprintf('logtext', f_l1b_scm)
	mrfprintf('logtext', fgm_nf_file)
	mrfprintf('logtext', scm_nf_file)
	mrfprintf('logtext', scm_cal_file)
	mrfprintf('logtext', dss_file)
	mrfprintf('logtext', defatt_file)
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
	fgm = mms_fsm_fgm_read( f_l1a_fgm, f_l2pre_fgm, t_orbit );
	scm = mms_fsm_scm_read( f_l1b_scm, t_orbit );
	
	% Add the model directory
	fgm.model_dir = fgm_model_dir;

%------------------------------------%
% Read Noise Floor                   %
%------------------------------------%
	% Create weight function
%	nf_fgm = mms_fsm_nf_read(fgm_nf_file);
%	nf_scm = mms_fsm_nf_read(scm_nf_file);

%	[~, ts_tt2000] = mms_parse_time(tstart);
%	w              = mms_fsm_nf_weight(nf_fgm, nf_scm, ts_tt2000);

	w = [];
	
%------------------------------------%
% Process Data                       %
%------------------------------------%
	[t_fsm, b_omb] = mms_fsm_l2plus_create( fgm, scm, duration, w);

%------------------------------------%
% SMPA & BCS                         %
%------------------------------------%

	% Rotation Matrices
	omb2smpa = mms_fg_xomb2smpa();
	bcs2smpa = mms_fg_xbcs2smpa(zmpa);
	smpa2bcs = permute(bcs2smpa, [2,1]);
	
	% OMB --> SMPA --> BCS
	b_smpa = mrvector_rotate(omb2smpa, b_omb);
	b_bcs  = mrvector_rotate(smpa2bcs, b_smpa);

	% Delete data
	clear omb2smpa bcs2smpa b_omb

%------------------------------------%
% Despin                             %
%------------------------------------%
	% Read DSS data
	sunpulse = mms_dss_read_sunpulse(dss_file, '', '', 'UniquePulse', true);
	
	% Rotation matrix
	smpa2dmpa = mms_dss_xdespin(sunpulse, t_fsm);
	
	% SMPA --> DMPA
	b_dmpa = mrvector_rotate(smpa2dmpa, b_smpa);

	% Delete data
	clear sunpulse smpa2dmpa

%------------------------------------%
% DMPA --> GSE & GSM                 %
%------------------------------------%

	% Rotation Matrix to GEI
	[b_gsm, b_gse] = mms_rot_despun2gsm(t_fsm, b_dmpa, defatt);
	
	% Delete data
	clear defatt
	
%------------------------------------%
% Write to File                      %
%------------------------------------%
	% Parent files
	parents      = [ f_l1a_fgm f_l2pre_fgm f_l1b_scm scm_cal_file defatt_file dss_file ];
	[~, parents] = cellfun(@fileparts, parents, 'UniformOutput', false);

	% Create a data structure
	data = struct( 'tt2000', t_fsm',  ...
	               'b',      single( sqrt( sum( b_dmpa.^2 ) ) ), ...
	               'b_bcs',  single(b_bcs)',  ...
	               'b_dmpa', single(b_dmpa)', ...
	               'b_gse',  single(b_gse)',  ...
	               'b_gsm',  single(b_gsm)'   ...
	             );
	clear t_fsm b_bcs b_smpa b_gse b_gsm

	% Write the file
	fsm_file = mms_fsm_l2plus_write( sc, mode, tstart, data, ...
	                                 'OptDesc', outoptdesc,  ...
	                                 'Parents', parents );

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