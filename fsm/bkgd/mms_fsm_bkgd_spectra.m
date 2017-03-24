%
% Name
%   mms_fsm_bkgd_spectra
%
% Purpose
%   Compute the amplitude, phase, and power of magnetic field measurments
%   using a sliding square window.
%
% Calling Sequence
%   DATA = mms_fsm_bkgd_compute(T_SSM, B, T)
%     Use the time in seconds since midnight T_SSM, and the three-component
%     magnetic field, B, to compute the amplitude, phase, and power as a
%     function of frequency over a sliding window of time interval T. Data
%     is returned in structures DATA.
%
% Parameters
%   T_SSM           in, required, type = double
%   B               in, required, type = 3xN double
%   T               in, required, type = double
%
% Returns
%   DATA            out, required, type=struct
%                   Fields are:
%                       't'     - Center times of the time series
%                       'f'     - Frequencies
%                       'amp'   - Wave amplitude
%                       'phase' - Phase
%                       'psd'   - Power spectral density
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2016-05-31      Written by Matthew Argall
%
function data = mms_fsm_bkgd_spectra(t_ssm, B, T)

%------------------------------------%
% FFT Setup                          %
%------------------------------------%
	
	% FFT parameters
	%   - FGM and SCM now have the same number of points
	nPts   = size(B, 2);
	dt     = median( diff(t_ssm) );
	sr     = round( 1.0 / dt );
	dt     = 1.0 / sr;
	df     = 1 / T;
	nFFT   = int32(T / dt);
	nShift = nFFT / int32(2);
	nMax   = fft_n_max( t_ssm, double(nFFT), double(nShift) );
	freqs  = fft_freqs(df, double(nFFT) );
	i0     = int32(1);
	i1     = nFFT;

	% Allocate memory
	nFreqs = length(freqs);
	t_fft  = zeros(1, nMax);
	amp    = zeros(3, nMax, nFreqs);
	phase  = zeros(3, nMax, nFreqs);
	pwr    = zeros(3, nMax, nFreqs);

%------------------------------------%
% Loop Over All Intervals            %
%------------------------------------%
	for ii = 1 : nMax
		% Fourier transform
		b_fft = fft( B(:,i0:i1), double(nFFT), 2 );

		% Centered time stamp
		t_fft(ii) = t_ssm(i0) + (t_ssm(i1) - t_ssm(i0)) / 2.0;

		% Amplitude
		amp(:,ii,2:end) = 2.0 * dt * abs( b_fft(:, 2:nFFT/2+1) );
		amp(:,ii,1)     =       dt * abs( b_fft(:, 1) );
		
		% Phase
		phase(:,ii,1:end) = atan2( imag(b_fft(:, 1:nFFT/2+1)), real(b_fft(:, 1:nFFT/2+1)) ) * 180.0/pi;
		
		% Power
		pwr(:,ii,2:end) = (2.0 * dt / double(nFFT)) * abs( b_fft(:, 2:nFFT/2+1).^2 );
		pwr(:,ii,1)     = (      dt / double(nFFT)) * abs( b_fft(:, 1).^2 );
		
		% Advance to next interval
		i0 = i0 + nShift;
		i1 = i0 + nFFT - 1;
	end

%------------------------------------%
% Output Structures                  %
%------------------------------------%
	data = struct( 't',     t_fft,               ...
	               'f',     single( freqs ),     ...
	               'amp',   single( amp ),   ...
	               'phase', single( phase ), ...
	               'psd',   single( pwr )    ...
	             );
end