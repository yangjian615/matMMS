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
%   2015-04-06      Written by Matthew Argall
%
function fsm_ql = mms_fsm_create_ql(sc, tstart, tend, varargin)

	% Inputs for testing
%	sc     = 'mms1';
%	tstart = '2015-08-15T00:00:00';
%	tend   = '2015-08-15T24:00:00';

%------------------------------------%
% Defaults                           %
%------------------------------------%

	% Default directories
	%   - SCM has only one mode starting 2015/08/02
	fgm_instr    = 'dfg';
	fgm_mode     = 'fast';
	fgm_optdesc  = '';
	save_dir     = '/nfs/fsm';
	scm_mode     = 'fast';
	scm_optdesc  = 'scf';
	scm_cal_dir  = '';
	
	% Default merging parameters
	duration = 300.0;
	f_min    = 1.0;
	f_max    = 5.0;
	
	% Get optional inputs
	nOptArgs = length(varargin);
	for ii = 1 : 2 : nOptArgs
		switch varargin{ii}
			case 'Duration'
				duration = varargin{ii+1};
			case 'f_max'
				f_max = varargin{ii+1};
			case 'f_min'
				f_min = varargin{ii+1};
			case 'fgm_instr'
				fgm_instr = varargin{ii+1};
			case 'fgm_mode'
				fgm_mode = varargin{ii+1};
			case 'fgm_optdesc'
				fgm_optdesc = varargin{ii+1};
			case 'SaveDir'
				save_dir = varargin{ii+1};
			case 'scm_mode'
				scm_mode = varargin{ii+1};
			case 'scm_optdesc'
				scm_optdesc = varargin{ii+1};
			case 'scm_cal_dir'
				scm_cal_dir = varargin{ii+1};
			otherwise
				error( ['Optional argument not recognized: "' varargin{ii} '".'] );
		end
	end
	
	% Default FG calibration directory
	if isempty(scm_cal_dir)
		scm_cal_dir = '/home/argall/data/mms/scm_cal/';
%		scm_cal_dir = fullfile('/nfs', 'scm_cal', sc);
	end
	
	% Output directory
	if isempty(save_dir)
		save_dir = '/nfs/fsm/';
	end
		
	% Constants
	fgm_level   = 'l1a';
	scm_instr   = 'scm';
	scm_level   = 'l1a';
	fgm_cal_dir = fullfile('/nfs', 'mag_cal');

%------------------------------------%
% Find & Prepare Data                %
%------------------------------------%
	% Find files
	files = mms_fsm_ql_files(sc, tstart, tend, fgm_instr, fgm_mode, fgm_optdesc, fgm_cal_dir, ...
	                                                      scm_mode, scm_optdesc, scm_cal_dir);

	% Read and calibrate FGM
	[tt2000_fgm, b_fgm_omb, sr_fgm] = mms_fsm_ql_prep_fgm(files, fgm_instr, tstart, tend);

	% Read SCM and transfer function
	[tt2000_scm, b_scm_omb, xfr_scm, sr_scm] = mms_fsm_ql_prep_scm(files, tstart, tend, duration);

%------------------------------------%
% Find Sampling Rate Differences     %
%------------------------------------%
	%
	% Fast and slow survey data overlap on the transition from
	% fast to slow. We want to make use of that data, so will be
	% processing fast and slow survey data separately, then
	% combining the results. Therefore, we do not need to check
	% sampling rate.
	%

%------------------------------------%
% Find Coninuous, Overlapping Data   %
%------------------------------------%

	% Convert data to seconds
	% t_ref     = MrCDF_Epoch_Compute([2015 03 17]);
	t_ref     = min( [tt2000_fgm(1) tt2000_scm(1)] );
	t_sec_fgm = MrCDF_epoch2sse(tt2000_fgm, t_ref);
	t_sec_scm = MrCDF_epoch2sse(tt2000_scm, t_ref);

	% Find overlapping intervals
	%   - Remove intervals of FGM that fall entirely within an SCM data gap
	%     (and vice versa).
	[fgm_int, scm_int] = MrIntervalsXY(t_sec_fgm, t_sec_scm, 'Remove', true, 'Sync', true);
	n_int              = size(fgm_int, 2);
	
	% Clear data that will not be used anymore
	clear t_ref t_sec_fgm t_sec_scm

%------------------------------------%
% Loop Over Intervals                %
%------------------------------------%

	%
	% TODO
	%   1) Noise floor
	%

	% Allocate memory to output
	t_merged = zeros(1, size(b_scm_omb,2), 'int64');
	b_merged = zeros(size(b_scm_omb), 'single');

	% Step through each interval
	for ii = 1 : n_int 
		% Absolute start and end indices.
		%   - From the beginning of the array, not the beginning of
		%     the current interval.
		is_fgm = fgm_int(1,ii);
		is_scm = scm_int(1,ii);
		ie_fgm = fgm_int(2,ii);
		ie_scm = scm_int(2,ii);
	
		% Extract the data for the current interval
		t_fgm = tt2000_fgm(   is_fgm:ie_fgm );
		t_scm = tt2000_scm(   is_scm:ie_scm );
		b_fgm = b_fgm_omb( :, is_fgm:ie_fgm );
		b_scm = b_scm_omb( :, is_scm:ie_scm );

		% Merge the data
		t_merged(is_scm:ie_scm)    = t_scm;
		b_merged(:, is_scm:ie_scm) = ...
			fsm_merge_v2(duration, b_fgm, b_scm, t_fgm, t_scm,    ...
			          'dt_fg',         1.0 / single( sr_fgm(1) ), ...
			          'dt_sc',         1.0 / single( sr_scm(1) ), ...
			          'f_max',         f_max, ...
			          'f_min',         f_min, ...
			          'ref_index_fg',  1, ...
			          'ref_index_sc',  1, ...
			          'transfr_fn_sc', xfr_scm);
	end
	clear t_fgm t_scm b_fgm b_scm sr_fgm sr_scm xfr_scm

	% Remove unwanted data
	igood    = find(t_merged ~= 0);
	t_merged = t_merged(igood);
	b_merged = b_merged(:,igood);

%------------------------------------%
% Rotate to DMPA                     %
%------------------------------------%
	% OMB -> DMPA
	[b_merged_dmpa, b_fgm_dmpa] ...
		= mms_fsm_ql_omb2smpa(files, tstart, tend, t_merged, b_merged, tt2000_fgm, b_fgm_omb);
	
	% Clear data
	clear tt2000_scm b_merged b_fgm_omb b_scm_omb

%------------------------------------%
% Rotate to GSE                      %
%------------------------------------%
	
%------------------------------------%
% Prepare Output                     %
%------------------------------------%
	% Parent files
	parents = { files.fgm, files.fgm_hical, files.fgm_local, files.fgm_stemp, ...
	            files.scm, files.scm_cal, files.defatt, files.dss };
	parents = [ parents{:} ];

	% Remove directories
	[~, names, ext] = cellfun(@fileparts, parents, 'UniformOutput', false);
	parents         = strcat( names, ext);
	
	% Create the output structure
	fsm_ql = struct( 'tt2000',     t_merged',    ...
	                 'tt2000_fgm', tt2000_fgm',     ...
	                 'b_fsm_dmpa', single( b_merged_dmpa' ), ...
	                 'b_fgm_dmpa', single( b_fgm_dmpa' ),    ...
	                 'sc',         sc,                       ...
	                 'instr',      [fgm_instr '-scm'], ...
	                 'mode',       scm_mode,       ...
	                 'tstart',     tstart,         ...
	                 'directory',  save_dir,       ...
	                 'parents',    { parents }     ...
	               );
end


%
% Name
%   mms_fsm_ql_files
%
% Purpose
%   Find files used in the merging process.
%
% Calling Sequence
%   FILES = mms_fsm_merge(SC. TSTART, TEND, FGM_INSTR, FGM_MODE, FGM_OPTDESC, FGM_CAL_DIR, ...
%                                                      SCM_MODE, SCM_OPTDESC, SCM_CAL_DIR, ...
%                                                      DEFATT_DIR, SUNPULSE_DIR);
%     Use the the MMS spacecraft number SC (e.g. 'mms1'), time interval
%     [TSTART, TEND); FGM and SCM telemetry modes FGM_MODE SCM_MODE; optional
%     descriptors FGM_OPTDESC and SCM_OPTDESC; directories in which to find calibration
%     files SCM_CAL_DIR and FGM_CAL_DIR; and directories in which to find definitive
%     attitude and sunpulse files DEFATT_DIR and SUNPULSE_DIR.
%
% Parameters
%   SC:             in, required, type=char
%   TSTART:         in, required, type=char
%   TEND:           in, required, type=char
%   FGM_INSTR:      in, required, type=char
%   FGM_MODE:       in, required, type=char
%   FGM_OPTDESC:    in, required, type=char
%   FGM_CAL_DIR:    in, required, type=char
%   SCM_MODE:       in, required, type=char
%   SCM_OPTDESC:    in, required, type=char
%   SCM_CAL_DIR:    in, required, type=char
%   DEFATT_DIR:     in, required, type=char
%   SUNPULSE_DIR:   in, required, type=char
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-10-27      Written by Matthew Argall
%
function files = mms_fsm_ql_files(sc, tstart, tend, fgm_instr, fgm_mode, fgm_optdesc, fgm_cal_dir, ...
                                                               scm_mode, scm_optdesc, scm_cal_dir)

	% Constants
	defatt_dir = fullfile('/nfs', 'ancillary', sc, 'defatt');
	hk_dir     = fullfile('/nfs', 'hk');

%------------------------------------%
% FGM Files                          %
%------------------------------------%
	% FGM
	[fgm, nfgm, fsrch] = mms_file_search(sc, fgm_instr, fgm_mode, 'l1a', ...
	                                     'OptDesc', fgm_optdesc,        ...
	                                     'TStart',  tstart,             ...
	                                     'TEnd',    tend);
	assert(nfgm > 0, ['No DFG file found: "' fsrch '".']);
	
	% FGM HiCal
	[fgm_hical, nfgm_hical, fsrch] = mms_file_search(sc, fgm_instr, 'hirangecal', 'l2pre', ...
	                                                 'RelaxedTStart', true,               ...
	                                                 'SDCroot',       fgm_cal_dir,        ...
	                                                 'SubDirs',       '',                 ...
	                                                 'TStart',        tstart,             ...
	                                                 'TEnd',          tend);
	assert(nfgm_hical > 0, ['No HiCal file found: "' fsrch '".']);
	
	% FGM LoCal
	[fgm_local, nfgm_local, fsrch] = mms_file_search(sc, fgm_instr, 'lorangecal', 'l2pre', ...
	                                                 'RelaxedTStart', true,               ...
	                                                 'SDCroot',       fgm_cal_dir,         ...
	                                                 'SubDirs',       '',                 ...
	                                                 'TStart',        tstart,             ...
	                                                 'TEnd',          tend);
	assert(nfgm_local > 0, ['No LoCal file found: "' fsrch '".']);
	
	% Sensor temperature files
	[fgm_stemp, nfgm_stemp, fsrch] = mms_file_search(sc, 'fields', 'hk', 'l1b',  ...
	                                          'OptDesc',      '10e',      ...
	                                          'SDCroot',      hk_dir,     ...
	                                          'TStart',       tstart,     ...
	                                          'TEnd',         tend);
	assert(nfgm_stemp > 0, ['No sun sensor file found: "' fsrch '".']);

%------------------------------------%
% SCM Files                          %
%------------------------------------%
	
	% SCM
	[scm, nscm, fsrch] = mms_file_search(sc, 'scm', scm_mode, 'l1a', ...
	                                     'OptDesc',   scm_optdesc,     ...
	                                     'TStart',    tstart,          ...
	                                     'TEnd',      tend);
	assert(nscm > 0, ['No SCM ' scm_mode ' file found: "' fsrch '".']);
	
	% SCM Cal File
	scm_ftest = fullfile(scm_cal_dir, [sc '_scm' sc(4) '_caltab_%Y%M%d%H%M%S_v*.txt']);
	[scm_cal, nscm_cal] = MrFile_Search(scm_ftest,                      ...
	                                    'Closest',      true,           ...
	                                    'TimeOrder',    '%Y%M%d%H%M%S', ...
	                                    'TStart',       tstart,         ...
	                                    'TEnd',         tend,           ...
	                                    'VersionRegex', 'v[0-9]');
	assert(nscm_cal > 0, ['No SCM calibration file found: "' scm_ftest '".']);

%------------------------------------%
% Attitude Files                     %
%------------------------------------%
	
	% Attitude files
	att_ftest = fullfile( defatt_dir, [upper(sc) '_DEFATT_%Y%D_%Y%D.V*'] );
	[defatt, ndefatt] = MrFile_Search(att_ftest,              ...
	                                  'Closest',      true,   ...
	                                  'TimeOrder',    '%Y%D', ...
	                                  'TStart',       tstart, ...
	                                  'TEnd',         tend,   ...
	                                  'VersionRegex', 'V[0-9]{2}');
	assert(ndefatt > 0, ['No definitive attitude file found: "' att_ftest '".']);

%------------------------------------%
% Sunpulse Files                     %
%------------------------------------%
	
	% Sunpulse files
	[dss, ndss, fsrch] = mms_file_search(sc, 'fields', 'hk', 'l1b',  ...
	                                     'OptDesc',      '101',      ...
	                                     'SDCroot',      hk_dir,     ...
	                                     'TStart',       tstart,     ...
	                                     'TEnd',         tend);
	assert(ndss > 0, ['No sun sensor file found: "' fsrch '".']);

%------------------------------------%
% Output Files                       %
%------------------------------------%

	% Prevent the formation of an array of structures by
	% create cell arrays.
	files = struct( 'nfgm',       nfgm,       ...
	                'nfgm_hical', nfgm_hical, ...
	                'nfgm_local', nfgm_local, ...
	                'nfgm_stemp', nfgm_stemp, ...
	                'nscm',       nscm,       ...
	                'nscm_cal',   nscm_cal,   ...
	                'ndefatt',    ndefatt,    ...
	                'ndss',       ndss,       ...
	                'fgm',        { fgm       }, ...
	                'fgm_hical',  { fgm_hical }, ...
	                'fgm_local',  { fgm_local }, ...
	                'fgm_stemp',  { fgm_stemp }, ...
	                'scm',        { scm       }, ...
	                'scm_cal',    { scm_cal   }, ...
	                'defatt',     { defatt    }, ...
	                'dss',        { dss       }  ...
	              );
end


%
% Name
%   mms_fsm_ql_prep_fgm
%
% Purpose
%   Read FGM L1A data and calibrate it, creating a data product in OMB coordinates.
%
% Calling Sequence
%   [TT2000_FGM, B_FGM_OMB, SR_FGM] = mms_fsm_ql_prep_fgm(FILES, FGM_INSTR, TSTART, TEND);
%     Using data files in the structure FILES from instrument FGM_INSTR between the
%     time interval of [TSTART, TEND), read and calibrate data, producing time stamps
%     TT2000_FGM, magnetic field in OMB B_FGM_OMB, and sample rate SR_FGM.
%
% Parameters
%   FILES:          in, required, type=struct
%                   Required fields:
%                       'fgm'        - FGM file name
%                       'fgm_stemp'  - Sensor temperature HK10E file
%                       'fgm_local'  - Lo-range calibration file
%                       'fgm_hical'  - Hi-range calibration file
%   FGM_INSTR:      in, required, type=char
%   TSTART:         in, required, type=char
%   TEND:           in, required, type=char
%
% Returns
%   TT2000_FGM      out, required, type=1xN INT64 (cdf_time_tt2000)
%   B_FGM_OMB       out, required, type=4xN dloat
%   SR_FGM          out, reuqired, type=1xN 
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-10-27      Written by Matthew Argall
%
function [tt2000_fgm, b_fgm_omb, sr_fgm] = mms_fsm_ql_prep_fgm(files, fgm_instr, tstart, tend)

%------------------------------------%
% Read & Calibrate FGM Data          %
%------------------------------------%
	% Read sensor temperature data
	hk_10e = mms_hk_read_0x10e(files.fgm_stemp, tstart, tend);
	if fgm_instr == 'afg'
		ttemp = hk_10e.tt2000;
		stemp = hk_10e.afg_stemp;
		etemp = hk_10e.afg_etemp;
	else
		ttemp = hk_10e.tt2000;
		stemp = hk_10e.dfg_stemp;
		etemp = hk_10e.dfg_etemp;
	end
	clear hk_10e

	% Uncalibrated FG in 123
	fgm_l1a = mms_fg_read_l1a(files.fgm, tstart, tend);
	
	% Read Calibration data
	[hiCal, hiConst] = mms_fg_read_cal(files.fgm_hical, tstart, tend);
	[loCal, loConst] = mms_fg_read_cal(files.fgm_local, tstart, tend);

	% Calibrate FG
	[b_fgm_omb, mpa] = mms_fg_calibrate(fgm_l1a.tt2000, fgm_l1a.b_123,    ...
	                                    fgm_l1a.tt2000_ts, fgm_l1a.range, ...
	                                    hiCal, loCal, hiConst, loConst, ...
	                                    ttemp, stemp, etemp);

	% Exctract other data
	tt2000_fgm = fgm_l1a.tt2000;
	sr_fgm     = fgm_l1a.sample_rate;
end


%
% Name
%   mms_fsm_ql_prep_scm
%
% Purpose
%   Read FGM L1A data and calibrate it, creating a data product in OMB coordinates.
%
% Calling Sequence
%   [TT2000_FGM, B_FGM_OMB, SR_FGM] = mms_fsm_ql_prep_scm(FILES, FGM_INSTR, TSTART, TEND);
%     Using data files in the structure FILES from instrument FGM_INSTR between the
%     time interval of [TSTART, TEND), read and calibrate data, producing time stamps
%     TT2000_FGM, magnetic field in OMB B_FGM_OMB, and sample rate SR_FGM.
%
% Parameters
%   FILES:          in, required, type=struct
%                   Required fields:
%                       'scm'      - SCM file name
%                       'scm_cal'  - Sensor temperature HK10E file
%   TSTART:         in, required, type=char
%   TEND:           in, required, type=char
%
% Returns
%   TT2000_SCM      out, required, type=1xN INT64 (cdf_time_tt2000)
%   B_SCM_OMB       out, required, type=3xN float
%   XFR_SCM         out, reuqired, type=3xN 
%   SR_SCM          out, reuqired, type=1xN 
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-10-27      Written by Matthew Argall
%
function [tt2000_scm, b_scm_omb, xfr_scm, sr_scm] = mms_fsm_ql_prep_scm(files, tstart, tend, duration)

%------------------------------------%
% Prep SCM Data                      %
%------------------------------------%
	
	% Uncalibrated SCM in 123
	scm_l1a = mms_sc_read_l1a(files.scm, tstart, tend);
	
	% Read calibration data
	[transfr_fn, freqs] = mms_sc_read_caltab(files.scm_cal);

	% Extract the time
	tt2000_scm = scm_l1a.tt2000;
	sr_scm     = scm_l1a.sample_rate;

	% Convert numbers to nano-Tesla
	%   - Call this OMB. The SCM team considers 123 to be orthogonalized already.
	%   - Technically, data in OMB is fully calibrated, but the remainder of our
	%     calibration process will be performed simultaneously as the data is merged.
	%   - SCM is inverted with respect to AFG and DFG. Negate it.
	b_scm_omb = -mms_sc_number2nT(scm_l1a.b_123);

	% Frequency resolution
	df   = 1.0 / duration;
	n_sc = duration * sr_scm(1);

	% Create the compensation function
	xfr_scm = mms_sc_tf_compensate(transfr_fn, freqs, double(n_sc), df);
end




%
% Name
%   mms_fsm_ql_omb2smpa
%
% Purpose
%   Transform the merged magnetic field from OMB to DMPA coordinates.
%
% Calling Sequence
%   [B_MERGED_DMPA, B_FGM_DMPA] = mms_fsm_ql_omb2smpa(FILES, TSTART, TEND, T_MERGED, B_MERGED, T_FGM, B_FGM);
%     Using data files in the structure FILES between the time interval of [TSTART, TEND),
%     transform the data from OMB to SMPA, then despin SMPA into DMPA. Data have
%     TT2000 time stamps of T_MERGED and T_FGM, respectively.
%
% Parameters
%   FILES:          in, required, type=struct
%                   Required fields:
%                       'defatt'  - Definitive attitude files
%                       'dss'     - Sunpulse times from HK101
%   TSTART:         in, required, type=char
%   TEND:           in, required, type=char
%   B_MERGED:       in, required, type=3xN float
%   B_FGM_OMB:      in, required, type=3xN float
%
% Returns
%   B_MERGED_DMPA   out, required, type=1xN INT64 (cdf_time_tt2000)
%   B_FGM_DMPA      out, required, type=3xN float
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-10-27      Written by Matthew Argall
%
function [b_merged_dmpa, b_fgm_dmpa] ...
	= mms_fsm_ql_omb2smpa(files, tstart, tend, t_merged, b_merged, t_fgm, b_fgm)

%------------------------------------%
% Rotate to DMPA                     %
%------------------------------------%
	% OMB --> SMPA
	omb2smpa      = mms_fg_xomb2smpa();
	b_merged_smpa = omb2smpa * b_merged;
	b_fgm_smpa    = omb2smpa * b_fgm;
	
	% Clear data
	clear b_merged b_fgm omb2smpa
	
	% Read data
	sunpulse = mms_dss_read_sunpulse( files.dss, tstart, tend, 'UniquePulse', true );
	attitude = mms_fdoa_read_defatt( files.defatt, tstart, tend );
	
	% SMPA -> DMPA
	if ~isempty(sunpulse)
		xsmpa2dmpa_fsm = mms_dss_xdespin( sunpulse, t_merged );
		xsmpa2dmpa_fgm = mms_dss_xdespin( sunpulse, t_fgm );
	elseif ~isempty(attitude)
		xsmpa2dmpa_fsm = mms_fdoa_xdespin( attitude, t_merged, 'P' );
		xsmpa2dmpa_fgm = mms_fdoa_xdespin( attitude, t_fgm, 'P' );
	else
		warning('FSM::Despin', 'No Sunpulse or Attitude data found. Cannot despin.');
	end

	% Despin
	b_merged_dmpa = mrvector_rotate( xsmpa2dmpa_fsm, b_merged_smpa );
	b_fgm_dmpa    = mrvector_rotate( xsmpa2dmpa_fgm, b_fgm_smpa );
end