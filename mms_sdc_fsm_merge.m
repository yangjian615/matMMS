%
% Name
%   mms_fsm_merge
%
% Purpose
% 	Merge MMS fluxgate and search coil magnetometer data in the frequency
% 	domain.
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
function [b_merged, t_merged, b_fg, t_fg, b_sc, t_sc] = mms_sdc_fsm_merge(sc, tstart, tend, varargin)

	sc       = 'mms2';
	fg_instr = 'dfg';
	tstart   = '2015-03-17T15:00:00';
	tend     = '2015-03-17T24:00:00';
	sdc_root = '/nfs/mmsa/sdc/';

	duration = 32.0;
	f_min    = 0.5;
	f_max    = 3.0;
	
	fg_instr   = 'dfg';
	fg_mode    = 'srvy';
	fg_level   = 'l1b';
	fg_optdesc = '';
	
	sc_instr   = 'scm';
	sc_mode    = 'comm';
	sc_level   = 'l1a';
	sc_optdesc = 'sc256';
	sc_cal_dir = '/home/argall/data/mms/scm_cal/';

%------------------------------------%
% Find FG Files                      %
%------------------------------------%
	% Directory with cal files
	fg_cal_dir = fullfile(sdc_root, 'mag_cal', sc, fg_instr, '%((hi|lo)%)rangecal', 'l2pre');

	% FG L1B Data File
	[l1a_fname, count, str] = mms_file_search(sc, fg_instr, fg_mode, fg_level, ...
	                                          'OptDesc',   fg_optdesc, ...
	                                          'SDCroot',   sdc_root, ...
	                                          'TStart',    tstart, ...
	                                          'TEnd',      tend);
	assert(count > 0, [upper(fg_instr) ' file not found: "' str '".']);

	% FG Hi Cal Files
	[fg_hiCal, count, str] = mms_file_search(sc, fg_instr, 'hirangecal', 'l2pre', ...
	                                         'SDCroot',   sdc_root, ...
	                                         'TStart',    tstart, ...
	                                         'TEnd',      tend, ...
	                                         'Directory', fg_cal_dir);
	assert(count > 0, [upper(fg_instr) ' hirange cal file not found: "' str '".']);

	% FG Lo Cal Files
	[fg_loCal, count, str] = mms_file_search(sc, fg_instr, 'lorangecal', 'l2pre', ...
	                                         'SDCroot',   sdc_root, ...
	                                         'TStart',    tstart, ...
	                                         'TEnd',      tend, ...
	                                         'Directory', fg_cal_dir);
	assert(count > 0, [upper(fg_instr) ' lorange cal file not found: "' str '".']);

%------------------------------------%
% Find SCM Files                     %
%------------------------------------%
	
	% SCM L1B Data File
	[l1b_fname, count, str] = mms_file_search(sc, sc_instr, sc_mode, sc_level, ...
	                                          'SDCroot',   sdc_root, ...
	                                          'TStart',    tstart, ...
	                                          'TEnd',      tend, ...
	                                          'OptDesc',   sc_optdesc);
	assert(count > 0, ['SCM file not found: "' str '".']);
	
	% SCM Cal File
	str                 = fullfile(sc_cal_dir, [sc '_' sc_instr sc(4) '_caltab_%Y%M%d%H%m%S_v*.txt']);
	[cal_fname, nFiles] = MrFile_Search( str, 'VersionRegex', '([0-9])' );
	assert(nFiles > 0, ['SCM cal file not found: "' str '".']);

keyboard

%------------------------------------%
% Get FGM and SCM Data               %
%------------------------------------%

	% Prepare fluxgate data to be merged
	%   1. Read DFG time, magnetic field, sample rate
	%   2. Calibrate mag data
	%   3. Rotate into SCM_123 frame
	[b_fg, t_fg, sr_fg]  = mms_fsm_prep_fg(sc, tstart, tend, fg_dir, fg_cal_dir);
	
	% Prepare searchcoil data to be merged
	%   1. Read SCM time, magnetic field, sample rate
	%   2. Convert SCM number to nanotesla
	%   3. Read transfer function and frequencies
	[b_sc, t_sc, sr_sc, transfr_sc, f_sc] = mms_fsm_prep_sc(sc, tstart, tend, sc_dir, sc_cal_dir);

%------------------------------------%
% Find Coninuous, Overlapping Data   %
%------------------------------------%

	% Convert data to seconds
	t_ref     = min( [t_fg(1) t_sc(1)] );
	t_sec_fgm = MrCDF_epoch2sse(t_fg, t_ref);
	t_sec_scm = MrCDF_epoch2sse(t_sc, t_ref);

	% Find overlapping intervals
	%   - Remove intervals of FGM that fall entirely within an SCM data gap
	%     (and vice versa).
	[fg_int, sc_int] = fsm_intervals(t_sec_fgm, t_sec_scm, 'Remove', true);
	n_int            = length(fg_int(1, :));
	clear t_sec_fgm t_sec_scm

%------------------------------------%
% Interpolate SCM Cal Table          %
%------------------------------------%

	% Frequency resolution
	df   = 1.0 / duration;
	n_sc = duration * sr_sc(1);

	% Create the compensation function
	tf_comp_sc = mms_sc_tf_compensate(transfr_sc, f_sc, n_sc, df);
	clear transfr_sc

%------------------------------------%
% Loop Over Intervals                %
%------------------------------------%

	% Allocate memory to output
	b_merged = zeros(size(b_sc));

	% Step through each interval
	for ii = 1 : n_int 
		% Extract the data for the current interval
		t_fgm = t_fg(fg_int(1,ii):fg_int(2,ii));
		t_scm = t_sc(sc_int(1,ii):sc_int(2,ii));
		b_fgm = b_fg(:, fg_int(1,ii):fg_int(2,ii));
		b_scm = b_sc(:, sc_int(1,ii):sc_int(2,ii));

		% Merge the data
		b_merged(:, sc_int(1,ii):sc_int(2,ii)) = ...
			fsm_merge(duration, b_fgm, b_scm, t_fgm, t_scm, ...
			          'dt_fg',         1.0 / sr_fg(1), ...
			          'dt_sc',         1.0 / sr_sc(1), ...
			          'f_max',         f_max, ...
			          'f_min',         f_min, ...
			          'ref_index_fg',  1, ...
			          'ref_index_sc',  1, ...
			          'transfr_fn_sc', tf_comp_sc);
	end
	
	% The time array is the same as that of SCM
	t_merged = t_sc;
end