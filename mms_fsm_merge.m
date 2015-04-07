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
function [b_merged, t_merged, b_fg, t_fg, b_sc, t_sc] = mms_fsm_merge(sc, tstart, tend, varargin)

	% Default directories
	current_dir  = pwd();
	fg_dir       = current_dir;
	sc_dir       = current_dir;
	fg_cal_dir   = current_dir;
	sc_cal_dir   = current_dir;
	
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
			case 'fg_dir'
				fg_dir = varargin{ii+1};
			case 'fg_cal_dir'
				fg_cal_dir = varargin{ii+1};
			case 'sc_dir'
				sc_dir = varargin{ii+1};
			case 'sc_cal_dir'
				sc_cal_dir = varargin{ii+1};
			otherwise
				error( ['Optional argument not recognized: "' varargin{ii} '".'] );
		end
	end

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
	% t_ref     = MrCDF_Epoch_Compute([2015 03 17]);
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