%
% Name
%   mms_sc_caltest_v2
%
% Purpose
%   Test SCM calibration.
%
%   Process:
%     1. Go through the UNH calibration process to produce l1b and ql data.
%        i)  Calibrate    point-by-point    with mms_sc_calibrate_v2.
%        ii) Calibrate interval-by-interval with mms_sc_calibrate
%               (see mms_sc_create_l1b to switch between the two)
%     2. Read SCM file for L1B and QL data
%     3. Compare
%       a) amplitude
%       b) phase
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
function [] = mms_sc_caltest_v2(sc, tstart, duration)
	
	if nargin < 3
		duration = 64.0;
	end
	if nargin < 2
		tstart = '2015-06-12T18:45:00';
	end
	if nargin < 1
		sc = 'mms2';
	end
	
	% Plot options
	tf_restore = false;
	tf_linear  = true;
	apply_win  = false;
	component  = 'Y';
	
	% Cal and attitude directories
	scm_cal_dir  = '/home/argall/data/mms/scm_cal/';
	attitude_dir = fullfile('/nfs', 'ancillary', sc, 'defatt');
	sunpulse_dir = fullfile('/nfs', 'hk');
	
%------------------------------------%
% Find Files                         %
%------------------------------------%
	% Reformat the time
	fdate  = MrTimeParser(tstart, '%Y-%M-%dT%H:%m:%S', '%Y%M%d');
	fstart = [tstart(1:10) 'T18:40:00'];
	fend   = [tstart(1:10) 'T19:00:00'];

	instr   = 'scm';
	mode    = 'comm';
	optdesc = 'sc128';

	% SCM L1A
	[l1a_files, count, fsrch] = mms_file_search(sc, instr, mode, 'l1a', ...
	                                            'OptDesc',   optdesc, ...
	                                            'TStart',    fstart, ...
	                                            'TEnd',      fend);
	assert(count > 0, ['No L1A file found: "' fsrch '".']);
	
	% SCM L1B
	[l1b_files, count, fsrch] = mms_file_search(sc, instr, mode, 'l1b', ...
	                                            'OptDesc',   optdesc, ...
	                                            'TStart',    fstart, ...
	                                            'TEnd',      fend);
	assert(count > 0, ['No L1B file found: "' fsrch '".']);
	
	% SCM Cal File
	scm_ftest = fullfile(scm_cal_dir, [sc '_' instr sc(4) '_caltab_%Y%M%d%H%M%S_v*.txt']);
	[scm_cal_file, count] = MrFile_Search(scm_ftest, ...
	                                      'Closest',      true, ...
	                                      'TimeOrder',    '%Y%M%d%H%M%S', ...
	                                      'TStart',       fstart, ...
	                                      'TEnd',         fend, ...
	                                      'VersionRegex', 'v[0-9]');
	assert(count > 0, ['No SCM calibration file found: "' scm_ftest '".']);
	
	% Attitude files
	att_ftest = fullfile( attitude_dir, [upper(sc) '_DEFATT_%Y%D_%Y%D.V*'] );
	[defatt_files, count] = MrFile_Search(att_ftest, ...
	                                      'Closest',      true, ...
	                                      'TimeOrder',    '%Y%D', ...
	                                      'TStart',       fstart, ...
	                                      'TEnd',         fend, ...
	                                      'VersionRegex', 'V[0-9]{2}');
	assert(count > 0, ['No definitive attitude file found: "' att_ftest '".']);
	
	% Sunpulse files
	[dss_files, count, fsrch] = mms_file_search(sc, 'fields', 'hk', 'l1b', ...
	                                            'OptDesc',      '101',  ...
	                                            'SDCroot',      sunpulse_dir, ...
	                                            'TStart',       fstart, ...
	                                            'TEnd',         fend);
	assert(count > 0, ['No sun sensor file found: "' fsrch '".']);

%------------------------------------%
% Read Data                          %
%------------------------------------%

	% Save data to mat file
	fmat = mms_construct_filename( sc, instr, mode, 'all',           ...
	                               'OptDesc', [optdesc '-calv1'], ...
	                               'TStart',  fdate,                 ...
	                               'Version', 'v0.0.0' );
	fmat(end-2:end) = 'mat';

	% Takes a while. Prefer to restore from matfile.
	if tf_restore
		load(fmat)
	else
		% Attitude and Sunpulse
		[attitude, att_hdr] = mms_fdoa_read_defatt( defatt_files, fstart, fend );
		sunpulse            = mms_dss_read_sunpulse( dss_files, fstart, fend, 'UniquePulse', true );
		
		% Read SCM L1B
		scm_l1b = mms_sc_read_l1b(l1b_files, fstart, fend);

		% Create L2 data 
		[t, ~, ~, b_dmpa, ~, b_omb] = mms_sc_create_l2( l1a_files, scm_cal_file, fstart, fend, ...
		                                                'Duration', duration, ...
		                                                'Attitude', attitude, ...
		                                                'SunPulse', sunpulse, ...
		                                                'zMPA',     att_hdr.zMPA(:,end) );

	
		% Rotate SCM data into DMPA
		omb2smpa   = mms_fg_xomb2smpa();
		b_smpa_scm = mrvector_rotate( omb2smpa, scm_l1b.b_omb );
		
		% Despin SCM data
		smpa2dmpa  = mms_dss_xdespin( sunpulse, scm_l1b.tt2000 );
		b_dmpa_scm = mrvector_rotate( smpa2dmpa, b_smpa_scm );
		
		% Create L2 structure
		scm_ql     = struct( 'tt2000', scm_l1b.tt2000, 'b_dmpa', b_dmpa_scm );
		clear omb2smpa b_smpa_scm smpa2dmpa b_dmpa_scm
		
		% Save data to matfile
		save(fmat, 't', 'b_dmpa', 'b_omb', 'scm_l1b', 'scm_ql', 'sunpulse', 'attitude', 'att_hdr');
	end

%------------------------------------%
% Extract Time Interval              %
%------------------------------------%
	% Find data in time range
	t0    = spdfparsett2000( [tstart '.000000000'] );
	t1    = t0 + int64( duration * 1e9 );
	ikeep = find( t >= t0 & t <= t1 );

	% Select the data
	tt2000     = t(ikeep);
	b_omb_unh  = b_omb(:, ikeep);
	b_dmpa_unh = b_dmpa(:, ikeep);
	
	b_omb_scm  = scm_l1b.b_omb(:, ikeep);
	b_dmpa_scm = scm_ql.b_dmpa(:, ikeep);

%------------------------------------%
% Prep for FFT                       %
%------------------------------------%

	% Frequency and time parameters
	%   - Frequency resolution
	%   - Sampling period
	%   - Nyquist frequency
	%   - Frequencies (+1 is for the DC component)
	dt   = double( median( diff(tt2000) ) ) * 1e-9;
	sr   = round(1.0 / dt);
	dt   = 1.0 / sr;

	df   = 1 / duration;
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

	% Fourier transform
	Bf_omb_unh  = fft( b_omb_unh  .* win, [], 2 );
	Bf_omb_scm  = fft( b_omb_scm  .* win, [], 2 );
	Bf_dmpa_unh = fft( b_dmpa_unh .* win, [], 2 );
	Bf_dmpa_scm = fft( b_dmpa_scm .* win, [], 2 );

	% Amplitude
	amp_omb_unh  = 2.0 * dt * abs( Bf_omb_unh(:, 1:NFFT/2) );
	amp_omb_scm  = 2.0 * dt * abs( Bf_omb_scm(:, 1:NFFT/2) );
	amp_dmpa_unh = 2.0 * dt * abs( Bf_dmpa_unh(:, 1:NFFT/2) );
	amp_dmpa_scm = 2.0 * dt * abs( Bf_dmpa_scm(:, 1:NFFT/2) );
	
	% Phase
	phase_omb_unh  = atan( imag(Bf_omb_unh)  ./ real(Bf_omb_unh)  ) * 180.0/pi;
	phase_omb_scm  = atan( imag(Bf_omb_scm)  ./ real(Bf_omb_scm)  ) * 180.0/pi;
	phase_dmpa_unh = atan( imag(Bf_dmpa_unh) ./ real(Bf_dmpa_unh) ) * 180.0/pi;
	phase_dmpa_scm = atan( imag(Bf_dmpa_scm) ./ real(Bf_dmpa_scm) ) * 180.0/pi;
	
	% Power (not including DC)
	pwr_omb_unh  = 2 * dt * abs( Bf_omb_unh( :, 2:NFFT/2).^2 );
	pwr_omb_scm  = 2 * dt * abs( Bf_omb_scm( :, 2:NFFT/2).^2 );
	pwr_dmpa_unh = 2 * dt * abs( Bf_dmpa_unh(:, 2:NFFT/2).^2 );
	pwr_dmpa_scm = 2 * dt * abs( Bf_dmpa_scm(:, 2:NFFT/2).^2 );

%------------------------------------%
% Create Time Series                 %
%------------------------------------%

	% L1A -> L1B
	bif_omb_unh  = ifft( Bf_omb_unh,  [], 2 );
	bif_omb_scm  = ifft( Bf_omb_scm,  [], 2 );
	bif_dmpa_unh = ifft( Bf_dmpa_unh, [], 2 );
	bif_dmpa_scm = ifft( Bf_dmpa_scm, [], 2 );

	% Convert time to datenumber for plotting purposes
	t_sec = double( tt2000 - tt2000(1) ) * 1e-9;
	
%------------------------------------%
% Plot OMB Results                   %
%------------------------------------%
	fig_omb = figure();
	
	% Take the correct plot handle
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
		plot( t_sec, b_omb_unh(1, :), t_sec, bif_omb_unh(1, :), t_sec, b_omb_scm(1, :), t_sec, bif_omb_scm(1, :) );
		legend('UNH', 'UNH (ifft)', 'SCM', 'SCM (ifft)' );
		title('SCM Calibration Comparison');
		ylabel('B (nT)');
		xlabel('Time (s)');
		
		subplot(4,1,2);
		pwr_plot( fout, amp_omb_unh(1,:), fout, amp_omb_scm(1,:) );
		legend('UNH', 'SCM');
		title('SCM Calibration Comparison');
		ylabel('X Amplitude');
		xlabel('Freuqency (Hz)');
	
		subplot(4,1,3);
		phs_plot( fout, phase_omb_unh(1,1:NFFT/2), fout, phase_omb_scm(1,1:NFFT/2) );
		ylabel('X Phase');
		ylim([-90, 90]);
		xlabel('Freuqency (Hz)');
	
		subplot(4,1,4);
		pwr_plot( fout(2:end), pwr_omb_unh(1,:), fout(2:end), pwr_omb_scm(1,:) );
		ylabel('X Power');
		xlabel('Freuqency (Hz)');
	end
	
	% Y-Component
	if component == 'Y'
		subplot(4,1,1);
		plot( t_sec, b_omb_unh(2, :), t_sec, bif_omb_unh(2, :), t_sec, b_omb_scm(2, :), t_sec, bif_omb_scm(2, :) );
		legend('UNH', 'UNH (ifft)', 'SCM', 'SCM (ifft)' );
		title('SCM Calibration Comparison');
		ylabel('B (nT)');
		xlabel('Time (s)');
		
		subplot(4,1,2);
		pwr_plot( fout, amp_omb_unh(2,:), fout, amp_omb_scm(2,:) );
		legend('UNH', 'SCM');
		title('SCM Calibration Comparison');
		ylabel('Y Amplitude');
		xlabel('Freuqency (Hz)');
	
		subplot(4,1,3);
		phs_plot( fout, phase_omb_unh(2,1:NFFT/2), fout, phase_omb_scm(2,1:NFFT/2) );
		ylabel('Y Phase');
		ylim([-90, 90]);
		xlabel('Freuqency (Hz)');
	
		subplot(4,1,4);
		pwr_plot( fout(2:end), pwr_omb_unh(2,:), fout(2:end), pwr_omb_scm(2,:) );
		ylabel('Y Power');
		xlabel('Freuqency (Hz)');
	end
	
	% Z-Component
	if component == 'Z'
		subplot(4,1,1);
		plot( t_sec, b_omb_unh(3, :), t_sec, bif_omb_unh(3, :), t_sec, b_omb_scm(3, :), t_sec, bif_omb_scm(3, :) );
		legend('UNH', 'UNH (ifft)', 'SCM', 'SCM (ifft)' );
		title('SCM Calibration Comparison');
		ylabel('B (nT)');
		xlabel('Time (s)');
		
		subplot(4,1,2);
		pwr_plot( fout, amp_omb_unh(3,:), fout, amp_omb_scm(3,:) );
		legend('UNH', 'SCM');
		title('SCM Calibration Comparison');
		ylabel('Z Amplitude');
		xlabel('Freuqency (Hz)');
	
		subplot(4,1,3);
		phs_plot( fout, phase_omb_unh(3,1:NFFT/2), fout, phase_omb_scm(3,1:NFFT/2) );
		ylabel('Z Phase');
		ylim([-90, 90]);
		xlabel('Freuqency (Hz)');
	
		subplot(4,1,4);
		pwr_plot( fout(2:end), pwr_omb_unh(3,:), fout(2:end), pwr_omb_scm(3,:) );
		ylabel('Z Power');
		xlabel('Freuqency (Hz)');
	end
	
%------------------------------------%
% Plot DMPA Results                  %
%------------------------------------%
	fig_dmpa = figure();
	
	% X-Component
	if component == 'X'
		subplot(4,1,1);
		plot( t_sec, b_dmpa_unh(1, :), t_sec, bif_dmpa_unh(1, :), t_sec, b_dmpa_scm(1, :), t_sec, bif_dmpa_scm(1, :) );
		legend('UNH', 'UNH (ifft)', 'SCM', 'SCM (ifft)' );
		title('SCM Calibration Comparison');
		ylabel('B (nT)');
		xlabel('Time (s)');
		
		subplot(4,1,2);
		pwr_plot( fout, amp_dmpa_unh(1,:), fout, amp_dmpa_scm(1,:) );
		legend('UNH', 'SCM');
		title('SCM Calibration Comparison');
		ylabel('X Amplitude');
		xlabel('Freuqency (Hz)');
	
		subplot(4,1,3);
		phs_plot( fout, phase_dmpa_unh(1,1:NFFT/2), fout, phase_dmpa_scm(1,1:NFFT/2) );
		ylabel('X Phase');
		ylim([-90, 90]);
		xlabel('Freuqency (Hz)');
	
		subplot(4,1,4);
		pwr_plot( fout(2:end), pwr_dmpa_unh(1,:), fout(2:end), pwr_dmpa_scm(1,:) );
		ylabel('X Power');
		xlabel('Freuqency (Hz)');
	end
	
	% Y-Component
	if component == 'Y'
		subplot(4,1,1);
		plot( t_sec, b_dmpa_unh(2, :), t_sec, bif_dmpa_unh(2, :), t_sec, b_dmpa_scm(2, :), t_sec, bif_dmpa_scm(2, :) );
		legend('UNH', 'UNH (ifft)', 'SCM', 'SCM (ifft)' );
		title('SCM Calibration Comparison');
		ylabel('B (nT)');
		xlabel('Time (s)');
		
		subplot(4,1,2);
		pwr_plot( fout, amp_dmpa_unh(2,:), fout, amp_dmpa_scm(2,:) );
		legend('UNH', 'SCM');
		title('SCM Calibration Comparison');
		ylabel('Y Amplitude');
		xlabel('Freuqency (Hz)');
	
		subplot(4,1,3);
		phs_plot( fout, phase_dmpa_unh(2,1:NFFT/2), fout, phase_dmpa_scm(2,1:NFFT/2) );
		ylabel('Y Phase');
		ylim([-90, 90]);
		xlabel('Freuqency (Hz)');
	
		subplot(4,1,4);
		pwr_plot( fout(2:end), pwr_dmpa_unh(2,:), fout(2:end), pwr_dmpa_scm(2,:) );
		ylabel('Y Power');
		xlabel('Freuqency (Hz)');
	end

	% Z-Component
	if component == 'Z'
		subplot(4,1,1);
		plot( t_sec, b_dmpa_unh(3, :), t_sec, bif_dmpa_unh(3, :), t_sec, b_dmpa_scm(3, :), t_sec, bif_dmpa_scm(3, :) );
		legend('UNH', 'UNH (ifft)', 'SCM', 'SCM (ifft)' );
		title('SCM Calibration Comparison');
		ylabel('B (nT)');
		xlabel('Time (s)');
		
		subplot(4,1,2);
		pwr_plot( fout, amp_dmpa_unh(3,:), fout, amp_dmpa_scm(3,:) );
		legend('UNH', 'SCM');
		title('SCM Calibration Comparison');
		ylabel('Z Amplitude');
		xlabel('Freuqency (Hz)');
	
		subplot(4,1,3);
		phs_plot( fout, phase_dmpa_unh(3,1:NFFT/2), fout, phase_dmpa_scm(3,1:NFFT/2) );
		ylabel('Z Phase');
		ylim([-90, 90]);
		xlabel('Freuqency (Hz)');
	
		subplot(4,1,4);
		pwr_plot( fout(2:end), pwr_dmpa_unh(3,:), fout(2:end), pwr_dmpa_scm(3,:) );
		ylabel('Z Power');
		xlabel('Freuqency (Hz)');
	end
end