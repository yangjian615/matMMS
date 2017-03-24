%
% Name
%   mms_fsm_fgm_compensate
%
% Purpose
%   Apply David Fischer's filter models to AFG/DFG data. Timing delays that are removed
%   in the L1A process are added back into the data before the models are applied.
%
% Calling Sequence
%   [T_FGM, B_FGM, T_SCM, B_SCM] = mms_fsm_bkgd_compute(MODEL, T_FGM, B_FGM, T_SCM, B_SCM)
%     Determine the remaining time delay between FGM and SCM and correct the
%     magnetic field data B_FGM and B_SCM using the filter model MODEL. Time stamps
%     for FGM and SCM, T_FGM and T_SCM, should have already been up- or down-sampled
%     to the target sampling frequency FS. Time stamps are corrected for the delay. 
%
% Parameters
%   MODEL           in, required, type = struct
%   T_SCM           in, required, type = in64 (cdf_time_tt2000)
%   B_SCM           in, required, type = Nx3 single
%   T_FGM           in, required, type = in64 (cdf_time_tt2000)
%   B_FGM           in, required, type = Nx3 single
%   FS              in, required, type = int16
%
% Returns
%   T_SCM           out, required, type = in64 (cdf_time_tt2000)
%   B_SCM           out, required, type = Nx3 single
%   T_FGM           out, required, type = in64 (cdf_time_tt2000)
%   B_FGM           out, required, type = Nx3 single
%
% MATLAB release(s) MATLAB 9.0.0.341360 (R2016a)
% Required Products None
%
% History:
%   2017-01-30      Adapted from David Fischer's merge_snippet.m - MRA
%
function [t_fgm, b_fgm, t_scm, b_scm] = mms_fsm_align_time(model, t_fgm, b_fgm, t_scm, b_scm, fs)
	% Closest neighbour fitting
	[t_fgm, b_fgm, t_scm, b_scm] = align_time( t_fgm, b_fgm, t_scm, b_scm );
	
	% Calculate remaining time delta between the instruments in nanoseconds
	dt0 = double( t_fgm(1) - t_scm(1) ) / 1e9;  
	mrfprintf('logtext', ['  ','Remaining timedelta: ',num2str(dt0*1e6),' us'])
	
	% Delay the fluxgate data by the found time and the remaining (up to now uncorrected) instrument time shift
	b_fgm = delay_xfg( b_fgm, model, dt0, fs);
	
	% Correct the time stamp by the delta between the instruments
	%   - This was no time shift, but fractional sample interpolation
	t_fgm = t_fgm - dt0*1e9;

	%time align again and cut to even length
	[t_fgm, b_fgm, t_scm, b_scm] = align_time(t_fgm, b_fgm, t_scm, b_scm);
	[t_fgm, b_fgm, t_scm, b_scm] = cut_even_length(t_fgm, b_fgm, t_scm, b_scm);
end