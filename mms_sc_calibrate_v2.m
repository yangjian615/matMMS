%
% Name
%   mms_sc_calibrate_v2
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
%   2015-07-30      Written by Matthew Argall
%
function [B_out, t_out] = mms_sc_calibrate_v2(B, sr, transfr_fn, f, duration)

	%
	% Calibration intervals are of length DURATION. Given data sampling rate,
	% SR, this defines the number of points in the FFT, the frequency resolution,
	% and the Nyquist frequency. In following the methodology of the FGM and SCM
	% merging process, the idea was to overlap each FFT by 1/4 interval and keep
	% only the middle quarter. This strategy is based on the use of a Tukey
	% tapering window, for which the middle quarter is unity, and the outer 3/4
	% tapers to zero.
	%
	% Unfortunately, for data in the spinning reference frame, which contains a large
	% amplitude sinusoid, the Tukey window imparts in a phase shift of the spin
	% freuqency and misaligns adjacent intervals. As such, the tapering windows
	% have been replaced by 1's until a more adequate solution can be found. The
	% overall method has been left in place.
	%

	if nargin < 5
		duration = 64.0;
	end
	
	% Parameters
	tf_window = false;
	
	% Number of points
	nPts = size(B, 2);
	
	% Allocate memory to output
	B_out = zeros(size(B));

%------------------------------------%
% Convert to nano-tesla              %
%------------------------------------%

	% Convert B to nano-Tesla
	%   - Negate and scale
	B_nT = mms_sc_number2nT(B);
	B_nT = -1.0 * B_nT;

%------------------------------------%
% FFT Parameters                     %
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
	assert( NFFT <= nPts, 'FFT length greater than data interval.' )

	% Interpolate the transfer function
	compensate_array = mms_sc_tf_compensate(transfr_fn, f, NFFT, df);

%------------------------------------%
% Outline Calibration Intervals      %
%------------------------------------%
	
	%
	% Indices
	%

	% FFT interval
	istart = 1;
	istop  = NFFT;
	
	% Middle portion to be kept
%	is_middle = floor(0.375 * NFFT) + 1;
%	ie_middle = floor(0.625 * NFFT);

	%
	% Keep the middle point, only
	%
	is_middle = NFFT / 2.0;
	ie_middle = NFFT / 2.0;
	nShift    = ie_middle - is_middle + 1;
	
	% Fill section
	is_fill = istart + is_middle - 1;
	ie_fill = istart + ie_middle - 1;

	% Tapering window?
	if tf_window
		win = repmat( window(@hamming, NFFT, 'periodic')', 3, 1 );
	else
		win = ones(3, NFFT);
	end

%------------------------------------%
% Step through Each Interval         %
%------------------------------------%
	count = 0;
	while istop <= nPts
	
% Print status
if mod(count, 5000) == 0
	fprintf('data = [%d, %d], fill = [%d, %d], middle = [%d, %d] of %d\n', ...
	        istart, istop, is_fill, ie_fill, is_middle, ie_middle, nPts);
end

	%------------------------------------%
	% Calibrate                          %
	%------------------------------------%
		% Fourier transform
		B_f = fft( B_nT(:, istart:istop) .* win, [], 2);
		
		% Apply transfer function
		B_f = B_f ./ compensate_array;
		
		% Transform back into time domain
		B_f = ifft(B_f, [], 2);
		
		% Store into output array
		B_out(:, is_fill:ie_fill) = B_f(:, is_middle:ie_middle);

	%------------------------------------%
	% Next Interval                      %
	%------------------------------------%
		istart = istart + nShift;
		istop  = istop  + nShift;
	
		% Fill section
		is_fill = is_fill + nShift;
		ie_fill = ie_fill + nShift;
		
		count = count + 1;
	end
	
	% Negate and scale
%	B_out = -1.0 * B_out;
end