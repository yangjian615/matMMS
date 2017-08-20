%
% Name
%   mms_fsm_test_psd
%
% Purpose
%   Compare a self-written power spectral density calculation with that of PWELCH.
%
%   Results:
%       Identical when square window is used. Other windows are untested.
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2016-08-23      Written by Matthew Argall
%
function [] = mms_fsm_test_psd()
%------------------------------------%
% Read the Data                      %
%------------------------------------%

	% File to constitute test
	scm_file = '/nfs/mms2/scm/srvy/l1b/scsrvy/2015/11/mms2_scm_srvy_l1b_scsrvy_20151108_v1.1.0.cdf';
	sr       = 32.0;
	dt       = 1.0 / sr;
	nfft     = 1024;
	nShift   = nfft / 2;

	% Read data
	[B, tt2000] = MrCDF_nRead(scm_file, 'mms2_scm_acb_scm123_scsrvy_srvy_l1b');

	% Select subset
	%   - We want 8 FFTs at 50% overlap
	niter  = 8.0;
	N      = (niter-1)*nShift + nfft;
	B      = B(:, 1:N)';
	tt2000 = tt2000(1:N);

	% Convert time to SSM
	t_ssm = MrCDF_epoch2ssm(tt2000);

%------------------------------------%
% Mimick PWelch                      %
%------------------------------------%

	% FFT parameters
	df  = 1.0 / (dt * nfft);
	fN  = sr / 2.0;
	f   = 0:df:fN;
	nf  = length(f);
	psd = zeros(niter, nf, 3);

	% Compute the spectra
	i0 = 1;
	i1 = nfft;
	for ii = 1 : niter
		% Compute the fourier transform
		b_fft = fft( B(i0:i1,:), nfft, 1 );

		% Power
		psd(ii,2:end,:) = (2.0 * dt / nfft) * abs( b_fft(2:nfft/2+1, :).^2 );
		psd(ii,1,:)     = (      dt / nfft) * abs( b_fft(1, :).^2 );
	
		% Next iteration
		i0 = i0 + nShift;
		i1 = i0 + nfft - 1;
	end

	% Average the data
	psd = squeeze( mean(psd, 1) );

%------------------------------------%
% PWelch                             %
%------------------------------------%

	% Spectra via PWelch
	%   - Use a square window
	psd_welch = zeros(nfft/2+1,3);
	[psd_welch(:,1), f_welch] = pwelch(B(:,1), ones(1,nfft), nfft/2.0, nfft, sr);
	[psd_welch(:,2), f_welch] = pwelch(B(:,2), ones(1,nfft), nfft/2.0, nfft, sr);
	[psd_welch(:,3), f_welch] = pwelch(B(:,3), ones(1,nfft), nfft/2.0, nfft, sr);

%------------------------------------%
% Plot Results                       %
%------------------------------------%

	% Create figure
	f1 = figure();

	% Plot time series
	subplot(4,1,1);
	p1 = plot(t_ssm, B);
	ylim([-100, 100]);
	
	% Plot x-PSD
	subplot(4,1,2);
	p2 = loglog(f, psd(:,1), f_welch, psd_welch(:,1));
	xlim([f(1), f(end)])
	ylim([1e-5, 1e5])
	
	subplot(4,1,3);
	p3 = loglog(f, psd(:,2), f_welch, psd_welch(:,2));
	xlim([f(1), f(end)])
	ylim([1e-5, 1e5])
	
	subplot(4,1,4);
	p4 = loglog(f, psd(:,3), f_welch, psd_welch(:,3));
	xlim([f(1), f(end)])
	ylim([1e-5, 1e5])
end