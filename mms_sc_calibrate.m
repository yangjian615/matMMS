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
%   NOTE: If DURATION is too short -- on the order of a
%         single spin (~20 seconds), there are not enough
%         spin revolutions in the data to accurately
%         determine
%
% Calling Sequence
%   B_OUT = mms_sc_calibrate(B, SR, TRANSFR_FN, F)
%     Calibrate search coil magnetic field data B sampled at sampling
%     rate SR using the trasfer function TRANSFR_FN as a function of
%     the frequencies F. Calibration intervals will be of length
%     32.0 seconds. This sets the minimum frequency and frequency
%     resolution of the output magnetic field B_OUT.
%
%   B_OUT = mms_sc_calibrate(B, SR, TRANSFR_FN, F, DURATION)
%     Specify the duration of the calibration interval.
%
% Parameters
%   B               in, required, type = 3xN double
%   SR              in, required, type = double
%   TRANSFR_FN      in, required, type = 3xM double
%   F               in, required, type = 3xM double
%   DURATION        in, optional, type = double, default = 600.0
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
%   2015-04-14      Do not truncate the output array. Instead, calibrate
%                     the ramaining points by extending backward from
%                     the end of the array. Read data elsewhere. - MRA
%
function [B_out, t_out] = mms_sc_calibrate(B, sr, transfr_fn, f, duration)

	if nargin < 5
		duration = 600.0
	end

%------------------------------------%
% Prep for FFT                       %
%------------------------------------%

	% Frequency and time parameters
	%   - Frequency resolution
	%   - Sampling period
	%   - Nyquist frequency
	%   - Frequencies (+1 is for the DC component)
	df   = 1 / duration;
	dt   = 1.0 / sr;
	NFFT = sr * duration;
	fN   = 1 / (2 * dt);
	fout = fN * linspace(0, 1, NFFT/2+1);

	% Make sure there are enough points in the array.
	assert( NFFT <= size(B, 2), 'FFT length greater than data interval.' )

	% Interpolate the transfer function
	compensate_array = mms_sc_tf_compensate(transfr_fn, f, NFFT, df);

%------------------------------------%
% Calibrate                          %
%------------------------------------%
	
	% Convert B to nano-Tesla
	%   - It appears as though the l1a data is already in nT
	B_out = mms_sc_number2nT(B);

	% Calibrate until end of dataset
	nPts   = size(B, 2);
	istart = 1;
	istop  = istart + NFFT - 1;
	
	% Loop through all FFT intervals
	while istop < nPts
		% Fourier transform
		B_f = fft(B_out(:, istart:istop), [], 2);

		% Apply the transfer function
		B_f(1, :) = B_f(1, :) ./ compensate_array(1, :);
		B_f(2, :) = B_f(2, :) ./ compensate_array(2, :);
		B_f(3, :) = B_f(3, :) ./ compensate_array(3, :);

		% Inverse transform
		B_out(:, istart:istop) = ifft(B_f, [], 2);

		% Move to next interval
		istart = istop + 1;
		istop  = istart + NFFT - 1;
	end

%------------------------------------%
% Calibrate End of Data              %
%------------------------------------%
	%
	% If an integral number of FFT windows do not fit
	% within the data interval, there will be 1 < N < NFFT
	% uncalibrated data points at the end of the array.
	% Calibrate them by taking an interval that extends
	% backward from the end of the array by NFFT points.
	%
%	nRemain = nPts - (istart - 1) + 1;
%	assert( nRemain > 1 && nRemain < NFFT, 'Improper calibration interval.' )
%
%	% Are there leftover data points
%	if nRemain > 0
%		% Extend backward from the end of the array
%		istop  = nPts;
%		istart = nPts - NFFT + 1;
%		
%		% Fourier transform
%		B_f = fft(B_out(:, istart:istop), [], 2);
%
%		% Apply the transfer function
%		B_f(1, :) = B_f(1, :) ./ compensate_array(1, :);
%		B_f(2, :) = B_f(2, :) ./ compensate_array(2, :);
%		B_f(3, :) = B_f(3, :) ./ compensate_array(3, :);
%
%		% Inverse transform
%		B_out(:, istart:istop) = ifft(B_f, [], 2);
%	end

%------------------------------------%
% Invert                             %
%------------------------------------%
	
	% Must multiply all components by -1.
	%   - MagCon Minutes 2015-03-25
%	B_out = -B_out;
end