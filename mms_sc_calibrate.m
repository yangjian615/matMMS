%
% Name
%   mms_sc_calibrate
%
% Purpose
%   Calibrate MMS search coil data.
%
%   Process:
%     1. Find and read calibration file
%     2. Create compensation array for transfer function
%     3. Convert numbers to nano-Tesla
%     4. Calibrate
%       a) Fourier transform intervals
%       b) Apply compensate array
%       c) Inverse fourier transform
%
% Calling Sequence
%   [B_OUT, T_OUT] = mms_sc_calibrate(B, T, SC, DURATION)
%     Calibrate search coil magnetic field data B with time tags T,
%     observed by spacecraft SC (e.g. 'mms1'). Calibration intervals
%     will be of length 32.0 seconds. This sets the minimum frequency
%     and frequency resolution of the output magnetic field B_OUT.
%     The returned type tags T_OUT are the same as T. B_OUT and
%     T_OUT will have lengths different from B and T if an integral
%     number of 32.0s intervals do not fit within T.
%
%   [B_OUT, T_OUT] = mms_sc_calibrate(B, T, SC, DURATION)
%     Specify the duration of the calibration interval.
%
% Parameters
%   B               in, required, type = 3xN double
%   T               in, required, type = 1xN int64 (cdf_time_tt2000)
%   SC              in, required, type = char
%   DURATION        in, optional, type = double, default = 32.0
%
% Returns
%   B               out, required, type=3xN double
%   T               out, optional, type=1xN int64
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-12      Written by Matthew Argall
%
function [B_out, t_out] = mms_sc_calibrate(B, t, sc, duration)

	if nargin < 4
		duration = 32.0
	end

%------------------------------------%
% Read the Transfer Function         %
%------------------------------------%
	cal_dir = '/Users/argall/Documents/Work/Data/MMS/SCM_Cal/';
	
	% Find the calibration file
	sc_cal_ftest = fullfile(cal_dir, [sc '_' instr sc(4) '_caltab_*.txt']);
	sc_cal_file  = dir(sc_cal_ftest);
	sc_cal_file  = fullfile(cal_dir, sc_cal_file.name);

	% Read the cal file
	[transfr_fn, f] = mms_sc_read_caltab(sc_cal_file);

%------------------------------------%
% Prep for FFT                       %
%------------------------------------%

	% Frequency and time parameters
	%   - Frequency resolution
	%   - Sampling period
	%   - Nyquist frequency
	%   - Frequencies (+1 is for the DC component)
	df   = 1 / duration;
	dt   = median( diff( MrCDF_epoch2sse( t ) ) );
	NFFT = duration / dt;
	fN   = 1 / (2 * dt);
	fout = fN * linspace(0, 1, NFFT/2+1);

	% Interpolate the transfer function
	%   - Should this really be NFFT-1?
	%       * There is no 0Hz calibration point in the transfer function. It is
	%         added manually.
	compensate_array = mms_sc_tf_compensate(transfr_fn, f, NFFT, df);

%------------------------------------%
% Calibrate                          %
%------------------------------------%
	
	% Convert B to nano-Tesla
	B_out = mms_sc_number2nT(B);

	% Calibrate until end of dataset
	nPts   = length(t);
	istart = 1;
	istop  = istart + NFFT - 1;
	
	% Loop through all FFT intervals
	while istop < nPts
		% Fourier transform
		B_f = fft(B_out(:, istart:istop), [], 2);

		% Apply the transfer function
		B_f(1, :) = B_f(1, :) .* compensate_array(1, :);
		B_f(2, :) = B_f(2, :) .* compensate_array(2, :);
		B_f(3, :) = B_f(3, :) .* compensate_array(3, :);

		% Inverse transform
		B_out(:, istart:istop) = ifft(B_f, [], 2);

		% Move to next interval
		istart = istop + 1;
		istop  = istart + NFFT - 1;
	end

%------------------------------------%
% Output                             %
%------------------------------------%
	
	% Truncate output
	t_out = t(1:istart-1);
	b_out = B_out(:, 1:istart-1);
end