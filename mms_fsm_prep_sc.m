%
% Name
%   mms_fsm_prep_sc
%
% Purpose
%   Prepare MMS search coil data to be merged.
%
% Calling Sequence
%   [B, T, SR, TRANSFR_FN, F] = mms_fsm_prep_sc(SC, TSTART, TEND, DATA_DIR, CAL_DIR)
%     Read and calibrate SCM data from spacecraft SC between the start and
%     end times TSTART and TEND. Data can be found in DATA_DIR and
%     the calibrated with files found in CAL_DIR. The uncalibrated magnetic
%     field data B is returned in the 123 sensor system with TT2000 time
%     stamps T and sample rate SR. TSTART. B can be calibrated in the
%     frequency domain with the transfer function TRANSFR_FN, which is a
%     function of frequencies F. TSTART and TEND should be formatted as
%     ISO-8601 strings.
%
% Parameters
%   SC:             in, required, type = char
%   TSTART          in, required, type = char
%   TEND            in, required, type = char
%   DATA_DIR        in, required, type = char
%   CAL_DIR         in, required, type = char
%
% Returns
%   B               out, required, type=3xN double
%   T               out, required, type=1xN int64
%   SR              out, required, type=1xM single
%   TRANSFR_FN      out, required, type=3xL single
%   F               out, required, type=3xL single
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-05      Written by Matthew Argall
%
function [b, t, sr, transfr_fn, f] = mms_fsm_prep_sc(sc, tstart, tend, data_dir, cal_dir)

	% Constants
	instr    = 'scm';
	mode     = 'comm';
	level    = 'l1a_sc128';

%------------------------------------%
% Read SC Data                       %
%------------------------------------%

	% Find the file
	files = mms_file_search(sc, instr, mode, level, ...
		                      'Directory', data_dir,  ...
													'Newest',    true,      ...
												  'TStart',    tstart,    ...
													'TEnd',      tend);
	
	% Make sure exactly one file was found
	assert( length(files) == 1, 'Can only handle one file at a time. Check time range.' )

	% Dissect the file name
	[sc, instr] = mms_dissect_filename(files{1});

	% Create variable names
	b_name  = mms_construct_varname(sc, instr, 'sc128_123');
	sr_name = mms_construct_varname(sc, instr, 'samplerate_sc128');

	% Read the data
	sc_data = spdfcdfread(files{1}, ...
											  'Variables', {'Epoch' b_name sr_name}, ...
											  'CombineRecords', true, ...
											  'KeepEpochAsIs', true);

	% Extract data
	%   - Sample rate is expected to be a float
	t  = sc_data{1}';
	b  = sc_data{2}';
	sr = single( sc_data{3}' );

	% Convert numbers to nano-Tesla
	b = mms_sc_number2nT(b);

%------------------------------------%
% Transfer Function                  %
%------------------------------------%

	% Find the calibration file
	sc_cal_ftest = fullfile(cal_dir, [sc '_' instr sc(4) '_caltab_*.txt']);
	sc_cal_file  = dir(sc_cal_ftest);
	sc_cal_file  = fullfile(cal_dir, sc_cal_file.name);

	% Read the cal file
	[transfr_fn, f] = mms_sc_read_caltab(sc_cal_file);
end