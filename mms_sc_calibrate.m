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
%   2015-05-08      Model FFT windowing after the FSM strategy. - MRA
%   2015-05-10      Replaced Hamming and Tukey windows with a constant window of 1.0. - MRA
%
function [B_out, t_out] = mms_sc_calibrate(B, sr, transfr_fn, f, duration)

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
		duration = 600.0;
	end
	
	% Number of points
	nPts = size(B, 2);

%------------------------------------%
% Invert                             %
%------------------------------------%
	
	% Must multiply all components by -1.
	%   - MagCon Minutes 2015-03-25
	%   - This is done only AFTER *ALL* data manipulation has been performed.
%	B_out = -B;

%------------------------------------%
% Convert to nano-tesla              %
%------------------------------------%

	% Convert B to nano-Tesla
	%   - Negate and scale
	B_nT = mms_sc_number2nT(B);
	B_nT = -1.0 * B_nT;

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
	assert( NFFT <= nPts, 'FFT length greater than data interval.' )

	% Interpolate the transfer function
	compensate_array = mms_sc_tf_compensate(transfr_fn, f, NFFT, df);

%------------------------------------%
% First Chunk of Data                %
%------------------------------------%
	B_out = zeros(size(B_nT));

	istart = 1;
	istop  = NFFT;

	% Tapering Window
	% win = repmat( window(@hamming, NFFT, 'periodic')', 3, 1 );
	win = ones(3, NFFT);

	% Fourier transform
	B_f = fft(B_nT(:, istart:istop) .* win, [], 2);

	% Apply the transfer function
	B_f(1, :) = B_f(1, :) ./ compensate_array(1, :);
	B_f(2, :) = B_f(2, :) ./ compensate_array(2, :);
	B_f(3, :) = B_f(3, :) ./ compensate_array(3, :);

	% Inverse transform
	B_f = ifft(B_f, [], 2) ./ win;
	B_out(:, istart:istop) = B_f;

%------------------------------------%
% Middle Chunks                      %
%------------------------------------%
	
	%
	% Indices
	%
	
	% Middle quarter to be kept
	is_middle = floor(0.375 * NFFT) + 1;
	ie_middle = floor(0.625 * NFFT);
	nShift    = ie_middle - is_middle + 1;
	
	% Advance by a quarter window
	istart = istart + nShift - 1;
	istop  = istop  + nShift - 1;
	
	% Fill section
	is_fill = istart + is_middle - 1;
	ie_fill = istart + ie_middle - 1;

	%
	% Calibrate
	%

	% Use a Tukey window
	% win = repmat( window(@tukeywin, NFFT, 0.75)', 3, 1 );
	win = ones(3, NFFT);

	% Loop through all FFT intervals
	while istop < nPts
		% Fourier transform
		B_f = fft(B_nT(:, istart:istop) .* win, [], 2);

		% Apply the transfer function
		B_f(1, :) = B_f(1, :) ./ compensate_array(1, :);
		B_f(2, :) = B_f(2, :) ./ compensate_array(2, :);
		B_f(3, :) = B_f(3, :) ./ compensate_array(3, :);

		% Inverse transform
		B_f = ifft(B_f, [], 2) ./ win;
		B_out(:, is_fill:ie_fill) = B_f(:, is_middle:ie_middle);

		% Save indices for plotting
		istart_old  = istart;
		istop_old   = istop;
		is_fill_old = is_fill;
		ie_fill_old = ie_fill;
		B_f_old     = B_f;

		% Move to Next Interval
		istart = istart + nShift - 1;
		istop  = istop  + nShift - 1;
	
		% Fill section
		is_fill = is_fill + nShift - 1;
		ie_fill = ie_fill + nShift - 1;
	end

%------------------------------------%
% Calibrate End of Data              %
%------------------------------------%
	% Last calibrated point
	iLast = istart;

	% Extend backward from the end of the array
	istop  = nPts;
	istart = nPts - NFFT + 1;
	
	% Take only those points that have not been calibrated
	nFill = nPts - iLast + 1;
	is_fill = iLast;
	ie_fill = nPts;
	
	is_middle = NFFT - nFill + 1;
	ie_middle = NFFT;
	
	% Tapering window
	win = repmat( window(@hamming, NFFT, 'periodic')', 3, 1 );
	
	% Fourier transform
	B_f = fft(B_nT(:, istart:istop) .* win, [], 2);

	% Apply the transfer function
	B_f(1, :) = B_f(1, :) ./ compensate_array(1, :);
	B_f(2, :) = B_f(2, :) ./ compensate_array(2, :);
	B_f(3, :) = B_f(3, :) ./ compensate_array(3, :);

	% Inverse transform
	B_f = ifft(B_f, [], 2) ./ win;
	B_out(:, is_fill:ie_fill) = B_f(:, is_middle:ie_middle);
end