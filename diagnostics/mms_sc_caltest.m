%
% Name
%   mms_sc_caltest
%
% Purpose
%   Test SCM calibration.
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
%   2015-07-23      Written by Matthew Argall
%
function [] = mms_sc_caltest(B, sr, transfr_fn, f, duration)	% Default FG calibration directory
	
	% Inputs
	sc        = 'mms2';
	instr     = 'scm';
	level     = 'l1a';
	mode      = 'comm';
	optdesc   = 'sc128';
	tstart    = '2015-06-12T19:30:00';
	tend      = '2015-06-12T24:00:00';
	tcal      = '2015-06-12T20:00:00';
	duration  = 40.0;
	
	% Plot options
	tf_linear = true;
	apply_win = false;
	component = 'Y';
	
	% Cal and attitude directories
	scm_cal_dir  = '/home/argall/data/mms/scm_cal/';
	attitude_dir = fullfile('/nfs', 'ancillary', sc, 'defatt');
	sunpulse_dir = fullfile('/nfs', 'hk');
	
%------------------------------------%
% Find Files                         %
%------------------------------------%
	% SCM L1A
	[l1a_files, count, fsrch] = mms_file_search(sc, instr, mode, 'l1a', ...
	                                            'OptDesc',   optdesc, ...
	                                            'TStart',    tstart, ...
	                                            'TEnd',      tend);
	assert(count > 0, ['No L1A file found: "' fsrch '".']);
	
	% SCM L1B
	[l1b_files, count, fsrch] = mms_file_search(sc, instr, mode, 'l1b', ...
	                                            'OptDesc',   optdesc, ...
	                                            'TStart',    tstart, ...
	                                            'TEnd',      tend);
	assert(count > 0, ['No L1B file found: "' fsrch '".']);
	
	% SCM Cal File
	scm_ftest = fullfile(scm_cal_dir, [sc '_' instr sc(4) '_caltab_%Y%M%d%H%M%S_v*.txt']);
	[scm_cal_file, count] = MrFile_Search(scm_ftest, ...
	                                      'Closest',      true, ...
	                                      'TimeOrder',    '%Y%M%d%H%M%S', ...
	                                      'TStart',       tstart, ...
	                                      'TEnd',         tend, ...
	                                      'VersionRegex', 'v[0-9]');
	assert(count > 0, ['No SCM calibration file found: "' scm_ftest '".']);
	
	% Attitude files
	att_ftest = fullfile( attitude_dir, [upper(sc) '_DEFATT_%Y%D_%Y%D.V*'] );
	[defatt_files, count] = MrFile_Search(att_ftest, ...
	                                      'Closest',      true, ...
	                                      'TimeOrder',    '%Y%D', ...
	                                      'TStart',       tstart, ...
	                                      'TEnd',         tend, ...
	                                      'VersionRegex', 'V[0-9]{2}');
	assert(count > 0, ['No definitive attitude file found: "' att_ftest '".']);
	
	% Sunpulse files
	[dss_files, count, fsrch] = mms_file_search(sc, 'fields', 'hk', 'l1b', ...
	                                            'OptDesc',      '101',  ...
	                                            'SDCroot',      sunpulse_dir, ...
	                                            'TStart',       tstart, ...
	                                            'TEnd',         tend);
	assert(count > 0, ['No sun sensor file found: "' fsrch '".']);

%------------------------------------%
% Read Data                          %
%------------------------------------%
	
	% Read L1A and L1B
	scm_l1a = mms_sc_read_l1a(l1a_files, tstart, tend);
	scm_l1b = mms_sc_read_l1b(l1b_files, tstart, tend);

	% Read calibration data
	[transfr_fn, freqs] = mms_sc_read_caltab(scm_cal_file);

	% Convert numbers to nano-Tesla
	%   - Call this OMB. The SCM team considers 123 to be orthogonalized already.
	%   - Technically, data in OMB is fully calibrated, but the remainder of our
	%     calibration process will be performed simultaneously as the data is merged.
	%   - SCM is inverted with respect to AFG and DFG. Negate it.
	b_omb_l1a = mms_sc_number2nT(scm_l1a.b_123);
	b_omb_l1a = -1.1 * b_omb_l1a;

	% Frequency resolution
	df   = 1.0 / duration;
	n_sc = duration * scm_l1a.sample_rate(1);

	% Create the compensation function
	tf_comp_sc = mms_sc_tf_compensate(transfr_fn, freqs, double(n_sc), df);

%------------------------------------%
% Prep for FFT                       %
%------------------------------------%

	% Frequency and time parameters
	%   - Frequency resolution
	%   - Sampling period
	%   - Nyquist frequency
	%   - Frequencies (+1 is for the DC component)
	sr   = scm_l1a.sample_rate(1);
	df   = 1 / duration;
	dt   = 1.0 / double(sr);
	NFFT = sr * duration;
	fN   = 1 / (2 * dt);
	fout = fN * linspace(0, 1.0, double(NFFT)/2.0);
	
	% Tapering window
	if apply_win
		win = repmat( window(@hamming, single(NFFT), 'periodic')', 3, 1 );
%		win = repmat( window(@tukeywin, single(NFFT), 0.75)', 3, 1 );
	else
		win = ones(3, NFFT);
	end

%------------------------------------%
% Take FFT                           %
%------------------------------------%

	% Get indicated data interval
	tvec   = datevec(tcal, 'yyyy-mm-ddTHH:MM:SS');
	t0     = spdfcomputett2000([tvec 0 0 0]);
	istart = find( scm_l1a.tt2000 >= t0, 1, 'first');
	istop  = istart + double(NFFT) - 1;

	% Fourier transform
	Bf_l1a = fft( b_omb_l1a(:, istart:istop)     .* win, [], 2 );
	Bf_l1b = fft( scm_l1b.b_omb(:, istart:istop) .* win, [], 2 );

	% Apply the transfer function
	Bf_l1a(1, :) = Bf_l1a(1, :) ./ tf_comp_sc(1, :);
	Bf_l1a(2, :) = Bf_l1a(2, :) ./ tf_comp_sc(2, :);
	Bf_l1a(3, :) = Bf_l1a(3, :) ./ tf_comp_sc(3, :);

	% Amplitude
	amp_l1a = 2.0 * dt * abs( Bf_l1a(:, 1:NFFT/2) );
	amp_l1b = 2.0 * dt * abs( Bf_l1b(:, 1:NFFT/2) );
	
	% Phase
	phase_l1a = atan( imag(Bf_l1a) ./ real(Bf_l1a) ) * 180.0/pi;
	phase_l1b = atan( imag(Bf_l1b) ./ real(Bf_l1b) ) * 180.0/pi;
	
	% Power (not including DC)
	pwr_l1a = 2 * dt * abs( Bf_l1a(:, 2:NFFT/2).^2 );
	pwr_l1b = 2 * dt * abs( Bf_l1b(:, 2:NFFT/2).^2 );

%------------------------------------%
% Create Time Series                 %
%------------------------------------%

	% L1A -> L1B
	t_unh = scm_l1a.tt2000(istart:istop);
	b_unh = ifft( Bf_l1a, [], 2 );
	
	% L1B
	t_scm   = scm_l1b.tt2000(istart:istop);
	b_scm_a = ifft( Bf_l1b, [], 2 );
	b_scm_b = scm_l1b.b_omb(:, istart:istop);
	
	clear scm_l1a scm_l1b

	% Convert time to datenumber for plotting purposes
	t_ref = t_unh(1);
	t_unh = double( t_unh - t_ref ) * 1e-9;
	t_scm = double( t_scm - t_ref ) * 1e-9;
%------------------------------------%
% Plot Results                       %
%------------------------------------%
	fig = figure();
	
	if tf_linear
		pwr_plot = @semilogy;
		phs_plot = @plot;
	else
		pwr_plot = @loglog;
		phs_plot = @semilogx;
	end
	
	% X-Component
	if component == 'X'
		subplot(4,1,1);
		plot( t_unh, b_unh(1, :), t_scm, b_scm(1, :) );
		legend('UNH', 'L1B');
		title('SCM Calibration Comparison');
		ylabel('B (nT)');
		xlabel('Time (TT2000)');
		
		subplot(4,1,2);
		pwr_plot( fout, amp_l1a(1,1:NFFT/2), fout, amp_l1b(1,1:NFFT/2) );
		legend('UNH', 'SCM');
		title('SCM Calibration Comparison');
		ylabel('X Amplitude');
		xlabel('Freuqency (Hz)');
	
		subplot(4,1,3);
		phs_plot( fout, phase_l1a(1,1:NFFT/2), fout, phase_l1b(1,1:NFFT/2) );
		ylabel('X Phase');
		ylim([-90, 90]);
		xlabel('Freuqency (Hz)');
	
		subplot(4,1,4);
		pwr_plot( fout(2:end), pwr_l1a(1,:), fout(2:end), pwr_l1b(1,:) );
		ylabel('X Power');
		xlabel('Freuqency (Hz)');
	end
	
	% Y-Component
	if component == 'Y'
		subplot(4,1,1);
		plot( t_unh, b_unh(2, :), t_scm, b_scm_a(2, :), t_scm, b_scm_b(2, :) );
		legend('UNH', 'SCM (ifft)', 'SCM' );
		title('SCM Calibration Comparison');
		ylabel('B (nT)');
		xlabel('Time (s)');
		
		subplot(4,1,2);
		pwr_plot( fout, amp_l1a(2,:), fout, amp_l1b(2,:) );
		legend('UNH', 'SCM');
		title('SCM Calibration Comparison');
		ylabel('Y Amplitude');
		xlabel('Freuqency (Hz)');
	
		subplot(4,1,3);
		phs_plot( fout, phase_l1a(2,1:NFFT/2), fout, phase_l1b(2,1:NFFT/2) );
		ylabel('Y Phase');
		ylim([-90, 90]);
		xlabel('Freuqency (Hz)');
	
		subplot(4,1,4);
		pwr_plot( fout(2:end), pwr_l1a(2,:), fout(2:end), pwr_l1b(2,:) );
		ylabel('Y Power');
		xlabel('Freuqency (Hz)');
	end
	
	% Z-Component
	if component == 'Z'
		subplot(4,1,1);
		plot( t_unh, b_unh(3, :), t_scm, b_scm(3, :) );
		legend('UNH', 'L1B');
		title('SCM Calibration Comparison');
		ylabel('B (nT)');
		xlabel('Time (TT2000)');
		
		subplot(4,1,2);
		pwr_plot( fout, amp_l1a(3,1:NFFT/2), fout, amp_l1b(3,1:NFFT/2) );
		legend('UNH', 'SCM');
		title('SCM Calibration Comparison');
		ylabel('Z Amplitude');
		xlabel('Freuqency (Hz)');
	
		subplot(4,1,3);
		phs_plot( fout, phase_l1a(3,1:NFFT/2), fout, phase_l1b(3,1:NFFT/2) );
		ylabel('Z Phase');
		ylim([-90, 90]);
		xlabel('Freuqency (Hz)');
	
		subplot(4,1,4);
		pwr_plot( fout(2:end), pwr_l1a(3,:), fout(2:end), pwr_l1b(3,:) );
		ylabel('Z Power');
		xlabel('Freuqency (Hz)');
	end
end