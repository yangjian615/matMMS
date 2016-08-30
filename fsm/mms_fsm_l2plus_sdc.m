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
%   'Duration':     in, required, type=double
%                   The duration of each merging interval. Sets the
%                     frequency resolution of the final dataset.
%   'f_max':        in, required, type=double, default=Nyquist frequency
%                   The maximum of the frequency range to merge.
%   'f_min':        in, required, type=double, default=df
%                   The minimum ( > 0 ) of the frequency range to merge.
%   'fg_dir':       in, required, type=char, default=pwd();
%                   Directory in which to find FGM data.
%   'fg_cal_dir':   in, required, type=char, default=pwd();
%                   Directory in which to find FGM calibration data.
%   'sc_dir':       in, required, type=char, default=pwd();
%                   Directory in which to find SCM data.
%   'sc_cal_dir':   in, required, type=char, default=pwd();
%                   Directory in which to find SCM calibration data.
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-12-08      Written by Matthew Argall
%
function status = mms_fsm_l2plus_sdc(sc, mode, tstart, varargin)
% fsm_file = mms_fsm_l2plus_sdc(fgm_file, scm_file, scm_cal_file, att_file, dss_file, duration)

	% Global path variables
	global cal_path_root data_path_root dropbox_root hk_root log_path_root unh_data_root

	% Establish data paths
	mms_fsm_init();
	
	% Defaults
	duration  = [];
	fgm_instr = 'dfg';
	no_log    = true;
	
	% Optional parameters
	nOptArgs = length(varargin);
	for ii = 1 : 2 : nOptArgs
		switch varargin{ii}
			case 'Duration'
				duration = varargin{ii+1};
			case 'FGMInstr'
				fgm_instr = varargin{ii+1};
			case 'NoLog'
				no_log = varargin{ii+1};
			otherwise
				error('MMS:FSM_L2Plus:SDC', 'Optional argument not recognized: "%s"', varargin{ii})
		end
	end

%------------------------------------%
% Check Inputs                       %
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
	outoptdesc = '';

%------------------------------------%
% Create Log File                    %
%------------------------------------%
	% Dissect the input time
	ahora = datestr(now(), 'YYYYMMDD_hhmmss');
	
	% Log file name
	logFile = [sc '_' outinstr '_' mode '_' outlevel '_' tstart '_' ahora];
	
	% Parse the input time
	tvec = mms_parse_time(tstart);
	
	
	% Directory
	if strcmp(mode, 'brst')
		logRelPath = fullfile(sc, outinstr, outmode, outlevel, tvec{1}, tvec{2}, tvec{3});
	else
		logRelPath = fullfile(sc, outinstr, outmode, outlevel, tvec{1}, tvec{2});
	end
	
	% Make the directory
	if ~exist( fullfile(log_path_root, logRelPath), 'dir' )
		mkdir(log_path_root, logRelPath);
	end
	
	% Complete log file
	logFile = fullfile(log_path_root, logRelPath, logFile);
	
	% Create the log file
	if ~no_log
		oLog = mrstdlog(logFile);
	end

%------------------------------------%
% FGM: Find File                     %
%------------------------------------%
	%
	% FGM has only srvy files, so no need to check for
	% fast and slow.
	%
	fgm_file = mms_latest_file(dropbox_root, sc, fgm_instr, mode, 'l2pre', tstart, ...
	                           'RootDir', data_path_root);

%------------------------------------%
% SCM: Find Files                    %
%------------------------------------%
	%
	% As of 2015-09-01, SCM has only srvy files. Before then,
	% it had slow and fast, but no srvy.
	%
	scm_file = mms_latest_file(dropbox_root, sc, 'scm', mode, 'l1a', tstart, ...
	                           'RootDir', data_path_root);
	
	%
	% SCM Cal file
	%
	scm_cal_dir = '/home/argall/data/mms/scm_cal';
	
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
	dss_file = mms_latest_file(dropbox_root, sc, 'fields', 'hk', 'l1b', tstart, ...
	                           'OptDesc', '101', ...
	                           'RootDir', hk_root);

%------------------------------------%
% FDOA: Find File                    %
%------------------------------------%
	defatt_file = mms_anc_search(dropbox_root, sc, 'defatt', tstart, 'RootDir', data_path_root);
	
	% No file found
	assert(~isempty(defatt_file), 'No DEFATT file found.')

%------------------------------------%
% Process Data                       %
%------------------------------------%
	% Write parents to log file
	mrfprintf('logtext', '')
	mrfprintf('logtext', '---------------------------------')
	mrfprintf('logtext', '| Parent Files                  |')
	mrfprintf('logtext', '---------------------------------')
	mrfprintf('logtext', fgm_file)
	mrfprintf('logtext', scm_file)
	mrfprintf('logtext', dss_file)
	mrfprintf('logtext', defatt_file)
	mrfprintf('logtext', '---------------------------------')
	mrfprintf('logtext', '')

%------------------------------------%
% Read Attitude and zMPA             %
%------------------------------------%
	[defatt, att_hdr] = mms_fdoa_read_defatt(att_file);
	zmpa              = att_hdr.zMPA(:,1);

%------------------------------------%
% Process Data                       %
%------------------------------------%
	[t_fsm, b_omb] = mms_fsm_l2plus_create(fgm_file, scm_file, scm_cal_file, zmpa, duration);

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
try
	parents      = { fgm_file scm_file scm_cal_file att_file{:} dss_file };
	[~, parents] = cellfun(@fileparts, parents, 'UniformOutput', false);
catch ME
	keyboard
end

	% Create a data structure
	data = struct( 'tt2000', t_fsm',  ...
	               'b_bcs',  single(b_bcs)',  ...
	               'b_dmpa', single(b_dmpa)', ...
	               'b_gse',  single(b_gse)',  ...
	               'b_gsm',  single(b_gsm)'   ...
	             );
	clear t_fsm b_bcs b_smpa b_gse b_gsm

	% Write the file
	fsm_file = mms_fsm_l2plus_write( parents, data );
end