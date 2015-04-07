%
% Name
%   mms_fsm_prep_fg
%
% Purpose
%   Prepare MMS fluxgate data to be merged.
%     1. Search for for spacecraft and time range provided
%     2. Read and correct sample rate (given as sec/sample)
%     3. Calibrate data -- offset, gain, orthogonalization
%     4. Rotate into SCM_123 frame.
%
% Calling Sequence
%   [B, T, SAMPLE_RATE] = mms_dss_despin(SC, TSTART, TEND, DATA_DIR, CAL_DIR)
%     Read and calibrate DFG data from spacecraft SC between the start and
%     end times TSTART and TEND. Data can be found in DATA_DIR and
%     calibrated with files found in CAL_DIR. The calibrated and
%     orthogonalized magnetic field data B is returned in the 123 sensor
%     system with TT2000 time stamps T and sample rate SAMPLE_RATE. TSTART
%     and TEND should be formatted as ISO-8601 strings.
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
%   SAMPLE_RATE     out, required, type=1xM single
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-05      Written by Matthew Argall
%
function [b, t, sample_rate] = mms_fsm_prep_fg(sc, tstart, tend, data_dir, cal_dir)

	% Constants
	instr = 'dfg';
	mode  = 'f128';
	level = 'l1a';
	
	% Find the files
	files = mms_file_search(sc, instr, mode, level, ...
		                      'Directory', data_dir,  ...
													'Newest',    true,      ...
                          'TStart',    tstart,    ...
													'TEnd',      tend);
	
	% Can only handle one file at a time for now
	assert(length(files) == 1, 'Cannot process more than one file at a time. Change time range.');	
	
	%------------------------------------%
	% Sample Rate                        %
	%------------------------------------%

	% Sample rate
	sr_name_fg  = mms_construct_varname(sc, instr, 'samplerate');
	sample_rate = spdfcdfread(files{1}, 'Variables', sr_name_fg);
	sample_rate = sample_rate';

	% FG reports the sample rate as seconds/sample. Convert to samples/second
	sample_rate = 1.0 ./ sample_rate;
	
	%------------------------------------%
	% Calibrate Data                     %
	%------------------------------------%

	% Read and calibrate FG data
	%   - Offset, gain, orthogonalization
	[b, t] = mms_fg_calibrate(files{1}, cal_dir);
	
	%------------------------------------%
	% Rotate into SCM_123                %
	%------------------------------------%
	fgm2scm = mms_instr_xxyz2instr('DFG_123', 'SCM_123');
	b       = mrvector_rotate(fgm2scm, b);
end