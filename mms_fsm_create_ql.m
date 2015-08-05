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
function fsm_ql = mms_fsm_create_l1b(sc, tstart, tend, varargin)

	% Default directories
	current_dir  = pwd();
	fg_instr     = 'dfg';
	fg_mode      = 'f128';
	fg_loCal_dir = '';
	fg_hiCal_dir = '';
	save_dir     = '';
	scm_mode     = 'comm';
	scm_optdesc  = 'sc128';
	scm_cal_dir  = '';
	attitude_dir = '';
	sunpulse_dir = '';
	
	% Default merging parameters
	duration = 32.0;
	f_min    = [];
	f_max    = [];
	
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
			case 'fg_mode'
				fg_mode = varargin{ii+1};
			case 'fg_instr'
				fg_instr = varargin{ii+1};
			case 'fg_loCal_dir'
				fg_loCal_dir = varargin{ii+1};
			case 'fg_hiCal_dir'
				fg_hiCal_dir = varargin{ii+1};
			case 'SaveDir'
				save_dir = varargin{ii+1};
			case 'scm_mode'
				sc_mode = varargin{ii+1};
			case 'scm_optdesc'
				sc_optdesc = varargin{ii+1};
			case 'scm_cal_dir'
				scm_cal_dir = varargin{ii+1};
			otherwise
				error( ['Optional argument not recognized: "' varargin{ii} '".'] );
		end
	end
	
	% Default FG calibration directory
	if isempty(fg_loCal_dir)
		fg_loCal_dir = fullfile('/nfs', 'mag_cal', sc, fg_instr, 'lorangecal', 'l2pre');
	end
	if isempty(fg_hiCal_dir)
		fg_hiCal_dir = fullfile('/nfs', 'mag_cal', sc, fg_instr, 'hirangecal', 'l2pre');
	end
	if isempty(scm_cal_dir)
		scm_cal_dir = fullfile('/nfs', 'scm_cal', sc);
	end
	
	% Output directory
	if isempty(save_dir)
		save_dir = '/nfs/edi/fsm/';
	end
	
	% Attitude directory
	if isempty(attitude_dir)
		attitude_dir = fullfile('/nfs', 'ancillary', sc, 'defatt');
	end
	
	% Constants
	fg_level  = 'l1a';
	scm_instr = 'scm';
	scm_level = 'l1a';

%------------------------------------%
% Find Files                         %
%------------------------------------%
	% FGM
	[fg_files, count, fsrch] = mms_file_search(sc, fg_instr, fg_mode, fg_level, ...
	                                           'TStart', tstart, ...
	                                           'TEnd',   tend);
	assert(count > 0, ['No DFG file found: "' fsrch '".']);
	
	% SCM
	[scm_files, count, fsrch] = mms_file_search(sc, scm_instr, scm_mode, scm_level, ...
	                                            'OptDesc',   scm_optdesc, ...
	                                            'TStart',    tstart, ...
	                                            'TEnd',      tend);
	assert(count > 0, ['No LoCal file found: "' fsrch '".']);
	
	% HI-CAL
	[hiCal_file, count, fsrch] = mms_file_search(sc, fg_instr, 'hirangecal', 'l2pre', ...
	                                             'Directory', fg_hiCal_dir, ...
	                                             'TStart',    tstart, ...
	                                             'TEnd',      tend);
	assert(count > 0, ['No HiCal file found: "' fsrch '".']);
	
	% LO-CAL
	[loCal_file, count, fsrch] = mms_file_search(sc, fg_instr, 'lorangecal', 'l2pre', ...
	                                             'Directory', fg_loCal_dir, ...
	                                             'TStart',    tstart, ...
	                                             'TEnd',      tend);
	assert(count > 0, ['No LoCal file found: "' fsrch '".']);
	
	% SCM Cal File
	scm_ftest = fullfile(scm_cal_dir, [sc '_' scm_instr sc(4) '_caltab_%Y%M%d%H%M%S_v*.txt']);
	[scm_cal_file, count] = MrFile_Search(scm_ftest, ...
	                                      'Closest',      true, ...
	                                      'TimeOrder',    '%Y%M%d%H%M%S', ...
	                                      'TStart',       tstart, ...
	                                      'TEnd',         tend, ...
	                                      'VersionRegex', 'v[0-9]');
	assert(count > 0, ['No SCM calibration file found: "' scm_ftest '".']);
	
	% Attitude files
	att_ftest = fullfile( attitude_dir, [upper(sc) '_DEFATT_%Y%D_%Y%D.V*'] );
	[defatt_files, count] = MrFile_Search(att_ftest, ...
	                                      'Closest',      true, ...
	                                      'TimeOrder',    '%Y%D', ...
	                                      'TStart',       tstart, ...
	                                      'TEnd',         tend, ...
	                                      'VersionRegex', 'V[0-9]{2}');
	assert(count > 0, ['No definitive attitude file found: "' att_ftest '".']);
	
	% Sunpulse files
	[dss_files, count, fsrch] = mms_file_search(sc, 'fields', 'hk', 'l1b', ...
	                                            'OptDesc',      '101',  ...
	                                            'SDCroot',      '/nfs/hk/', ...
	                                            'TStart',       tstart, ...
	                                            'TEnd',         tend);
	assert(count > 0, ['No sun sensor file found: "' fsrch '".']);

%------------------------------------%
% Read & Calibrate FGM Data          %
%------------------------------------%
	% Uncalibrated FG in 123
	fg_l1a = mms_fg_read_l1a(fg_files, tstart, tend);
	
	% Read Calibration data
	hiCal = mms_fg_read_cal(hiCal_file, tstart, tend);
	loCal = mms_fg_read_cal(loCal_file, tstart, tend);
	
	% Calibrate FG
	[b_fg_omb, mpa] = mms_fg_calibrate(fg_l1a.b_123, fg_l1a.tt2000, ...
	                                   fg_l1a.range, fg_l1a.tt2000_ts, hiCal, loCal);

	% Exctract other data
	tt2000_fg = fg_l1a.tt2000;
	sr_fg     = fg_l1a.sample_rate;

	% Clear data that will not be used anymore.
	clear fg_l1a hiCal loCal

%------------------------------------%
% Prep SCM Data                      %
%------------------------------------%
	
	% Uncalibrated SCM in 123
	scm_l1a = mms_sc_read_l1a(scm_files, tstart, tend);
	
	% Read calibration data
	[transfr_fn, freqs] = mms_sc_read_caltab(scm_cal_file);

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
	tf_comp_sc = mms_sc_tf_compensate(transfr_fn, freqs, double(n_sc), df);
	
	% Clear data that will not be used anymore
	clear scm_l1a transfr_sc freqs df n_sc

%------------------------------------%
% Find Coninuous, Overlapping Data   %
%------------------------------------%

	% Convert data to seconds
	% t_ref     = MrCDF_Epoch_Compute([2015 03 17]);
	t_ref     = min( [tt2000_fg(1) tt2000_scm(1)] );
	t_sec_fgm = MrCDF_epoch2sse(tt2000_fg,  t_ref);
	t_sec_scm = MrCDF_epoch2sse(tt2000_scm, t_ref);

	% Find overlapping intervals
	%   - Remove intervals of FGM that fall entirely within an SCM data gap
	%     (and vice versa).
	[fg_int, sc_int] = MrIntervalsXY(t_sec_fgm, t_sec_scm, 'Remove', true);
	n_int            = length(fg_int(1, :));
	
	% Clear data that will not be used anymore
	clear t_ref t_sec_fgm t_sec_scm

%------------------------------------%
% Loop Over Intervals                %
%------------------------------------%

	%
	% TODO
	%   1) Check for sampling rate changes.
	%   2) Noise floor
	%

	% Allocate memory to output
	b_merged = zeros(size(b_scm_omb));

	% Step through each interval
	for ii = 1 : n_int 
		% Find the closest starting index at the beginning of each interval
		[is_fgm, is_scm] = fsm_start_index( tt2000_fg(  fg_int(1,ii):fg_int(2,ii) ), ...
		                                    tt2000_scm( sc_int(1,ii):sc_int(2,ii) ), ...
		                                    sr_fg(1), single( sr_scm(1) ) );
	
		% Absolute start and end indices.
		is_fgm = is_fgm + fg_int(1,ii) - 1;
		is_scm = is_scm + sc_int(1,ii) - 1;
		ie_fgm = fg_int(2,ii);
		ie_scm = sc_int(2,ii);
	
		% Extract the data for the current interval
		t_fgm = tt2000_fg(    is_fgm:ie_fgm );
		t_scm = tt2000_scm(   is_scm:ie_scm );
		b_fgm = b_fg_omb(  :, is_fgm:ie_fgm );
		b_scm = b_scm_omb( :, is_scm:ie_scm );

		% Merge the data
		b_merged(:, is_scm:ie_scm) = ...
			fsm_merge_v2(duration, b_fgm, b_scm, t_fgm, t_scm, ...
			          'dt_fg',         1.0 / sr_fg(1), ...
			          'dt_sc',         1.0 / single( sr_scm(1) ), ...
			          'f_max',         f_max, ...
			          'f_min',         f_min, ...
			          'ref_index_fg',  1, ...
			          'ref_index_sc',  1, ...
			          'transfr_fn_sc', tf_comp_sc);
	end

%------------------------------------%
% Rotate to DMPA                     %
%------------------------------------%
	% OMB --> SMPA
	omb2smpa      = mms_fg_xomb2smpa();
	b_merged_smpa = omb2smpa * b_merged;
	b_fg_smpa     = omb2smpa * b_fg_omb;
	
	% Clear data
	clear b_merged b_fg_omb omb2smpa
	
	% Read data
	sunpulse = mms_dss_read_sunpulse( dss_files, tstart, tend, 'UniquePulse', true );
	attitude = mms_fdoa_read_defatt( defatt_files, tstart, tend );
	
	% SMPA -> DMPA
	if ~isempty(sunpulse)
		xsmpa2dmpa_fsm = mms_dss_xdespin( sunpulse, tt2000_scm );
		xsmpa2dmpa_fgm = mms_dss_xdespin( sunpulse, tt2000_fg );
	elseif ~isempty(attitude)
		xsmpa2dmpa_fsm = mms_fdoa_xdespin( attitude, tt2000_scm, 'L' );
		xsmpa2dmpa_fgm = mms_fdoa_xdespin( attitude, tt2000_fg, 'L' );
	else
		warning('FSM::Despin', 'No Sunpulse or Attitude data found. Cannot despin.');
	end

	% Despin
	b_merged_dmpa = mrvector_rotate( xsmpa2dmpa_fsm, b_merged_smpa );
	b_fg_dmpa     = mrvector_rotate( xsmpa2dmpa_fgm, b_fg_smpa );
	
	% Clear data
	clear xsmpa2dmpa_fsm xsmpa2dmpa_fgm b_merged_smpa b_fg_smpa attitude sunpulse
	
%------------------------------------%
% Prepare Output                     %
%------------------------------------%
	% Parent files
	if ischar(fg_files)
		fg_files = { fg_files };
	end
	if ischar(scm_files)
		scm_files = { scm_files };
	end
	if ischar(dss_files)
		dss_files = { dss_files };
	end
	parents = { fg_files{:} scm_files{:} dss_files{:} defatt_files{:} hiCal_file loCal_file scm_cal_file };
	[~, names, ext] = cellfun(@fileparts, parents, 'UniformOutput', false);
	parents         = strcat( names, ext);
	
	% Create the output structure
	fsm_ql = struct( 'tt2000',     tt2000_scm',    ...
	                 'tt2000_fgm', tt2000_fg',     ...
	                 'b_fsm_dmpa', single( b_merged_dmpa' ), ...
	                 'b_fgm_dmpa', single( b_fg_dmpa' ),     ...
	                 'sc',         sc,             ...
	                 'instr',      [fg_instr '-' scm_instr], ...
	                 'mode',       sc_mode,        ...
	                 'tstart',     tstart,         ...
	                 'directory',  save_dir,       ...
	                 'parents',    { parents }     ...
	               );

	% Write to file
	fname_fsm = mms_fsm_write_ql( fsm_ql );
end