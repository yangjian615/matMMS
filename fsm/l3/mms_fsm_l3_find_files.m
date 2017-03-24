%
% Name
%   mms_fsm_l3_find_files
%
% Purpose
%   Merge MMS fluxgate and search coil magnetometer data in the frequency
%   domain.
%
% Calling Sequence
%   FILES = mms_fsm_find_files(SC, FGM_INSTR, MODE, TSTART, TEND, SCM_CAL_DIR);
%     Given the MMS spacecraft number SC (e.g. 'mms1'), FGM instrument,
%     FGM_INSTR, data rate mode, MODE, and a time interval [TSTART, TEND),
%     gather all of the files required to produce FSM L3 data and return
%     them in the structure FILES.
%
% Parameters
%   SC:             in, required, type=char
%   FGM_INSTR:      in, required, type=char
%   MODE:           in, required, type=char
%   TSTART:         in, required, type=char
%   TEND:           in, required, type=char
%
% Returns
%   FILES           out, required, type=struct
%                   A structure with the following fields:
%                       'fgm_l2pre'  -  files from afg/dfg at level l2pre
%                       'fgm_l1a'    -  files from afg/dfg at level l1a
%                       'fgm_nf'     -  files with noise floor from afg/dfg
%                       'scm_l1b'    -  files from scm at level l1b
%                       'scm_nf'     -  files with noise floor from scm
%                       'scm_cal'    -  files with calibration data from scm
%                       'dss'        -  files from the digital sun sensor
%                       'defatt'     -  files with definitive attitude data
%
% MATLAB release(s) MATLAB 9.0.0.341360 (R2016a)
% Required Products None
%
% History:
%   2017-01-30      Written by Matthew Argall
%
function files = mms_fsm_l3_find_files(sc, fgm_instr, mode, tstart, tend, scm_cal_dir)

	% Global path variables
	global cal_path_root data_path_root dropbox_root hk_root log_path_root unh_data_root

	if nargin() < 6
		scm_cal_dir = '/home/argall/data/lpp_scm/cal';
	end

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
	
	% No file found
	assert( ~isempty(dss_file), 'No DSS file found.' )

%------------------------------------%
% FDOA: Find File                    %
%------------------------------------%
	defatt_file = mms_anc_search(dropbox_root, sc, 'defatt', tstart, 'RootDir', data_path_root);
	
	% No file found
	assert(~isempty(defatt_file), 'No DEFATT file found.')

%------------------------------------%
% Return a Structure                 %
%------------------------------------%
	files = struct( 'fgm_l2pre', f_l2pre_fgm,    ...
	                'fgm_l1a',   f_l1a_fgm,      ...
	                'fgm_nf',    fgm_nf_file,    ...
	                'scm_l1b',   f_l1b_scm,      ...
	                'scm_nf',    scm_nf_file,    ...
	                'scm_cal',   scm_cal_file,   ...
	                'dss',       dss_file,       ...
	                'defatt',    { defatt_file } ...
	              )
end