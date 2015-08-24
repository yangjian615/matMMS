%
% Name
%   mms_fg_caltest_tamphase
%
% Purpose
%   Compare UNH and FGM calibration results in the time
%   and frequency domains.
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-08-23      Written by Matthew Argall
%

get_data  = false;
component = 'Z';
tf_linear = false;
apply_win = false;

if get_data
%------------------------------------%
% Inputs                             %
%------------------------------------%
	sc         = 'mms2';
	tstart     = '2015-06-30T00:00:00';
	tend       = '2015-06-30T24:00:00';
	
	% Directories
	fg_cal_dir = '/nfs/mag_cal';
	hk_dir     = '/nfs/hk/';
	att_dir    = fullfile('/nfs', 'ancillary', sc, 'defatt');

	% Instrument parameters
	fg_instr   = 'dfg';
	fg_mode    = 'f128';
	fg_optdesc = '';
	
	% FFT Params
	T         = 64.0;
	fft_start = '2015-06-30T17:22:00';

%------------------------------------%
% Find Files                         %
%------------------------------------%
	% FG QL Data File
	[fg_ql_fname, count, str] = mms_file_search(sc, fg_instr, 'srvy', 'ql', ...
	                                            'TStart',    tstart,        ...
	                                            'TEnd',      tend);
	assert(count > 0, ['FG file not found: "' str '".']);


	% FG L1A Data File
	[fg_l1a_fname, count, str] = mms_file_search(sc, fg_instr, fg_mode, 'l1a', ...
	                                             'TStart',    tstart,          ...
	                                             'TEnd',      tend);
	assert(count > 0, ['FG file not found: "' str '".']);
	
	% FG Hi-Cal File
	[fg_hical_fname, count, str] = mms_file_search(sc, fg_instr, 'hirangecal', 'l2pre', ...
	                                               'RelaxedTStart', true,               ...
	                                               'SDCroot',       fg_cal_dir,         ...
	                                               'SubDirs',       '',                 ...
	                                               'TStart',        tstart,             ...
	                                               'TEnd',    tend);
	assert(count > 0, ['FG Hi-Cal file not found: "' str '".']);
	
	% FG Lo-Cal File
	[fg_local_fname, count, str] = mms_file_search(sc, fg_instr, 'lorangecal', 'l2pre', ...
	                                               'RelaxedTStart', true,               ...
	                                               'SDCroot',       fg_cal_dir,         ...
	                                               'SubDirs',       '',                 ...
	                                               'TStart',        tstart,             ...
	                                               'TEnd',          tend);
	assert(count > 0, ['FG Lo-Cal file not found: "' str '".']);
	
	% HK 0x10e files
	[hk_0x10e_fname, count, str] = mms_file_search(sc, 'fields', 'hk', 'l1b', ...
	                                                'OptDesc', '10e',         ...
	                                                'SDCroot', '/nfs/hk/',    ...
	                                                'TStart',  tstart,        ...
	                                                'TEnd',    tend);
	assert(count > 0, ['HK 0x10e file not found: "' str '".']);
	
	
	% Attitude files
	att_ftest = fullfile( att_dir, [upper(sc) '_DEFATT_%Y%D_%Y%D.V*'] );
	[defatt_files, count] = MrFile_Search(att_ftest,              ...
	                                      'Closest',      true,   ...
	                                      'TimeOrder',    '%Y%D', ...
	                                      'TStart',       tstart, ...
	                                      'TEnd',         tend,   ...
	                                      'VersionRegex', 'V[0-9]{2}');
	assert(count > 0, ['No definitive attitude file found: "' att_ftest '".']);

%------------------------------------%
% Read Data                          %
%------------------------------------%
	
	% Attitude
	[attitude, att_hdr] = mms_fdoa_read_defatt(defatt_files, tstart, tend);
	
	% FG L1A
	[t_fg, ~, ~, b_fg_dmpa] ...
		= mms_fg_create_l2(fg_l1a_fname, fg_hical_fname, fg_local_fname, tstart, tend, ...
		                   'Attitude', attitude, ...
		                   'hk_file',  hk_0x10e_fname);
	
	% FG QL
	fg_ql = mms_fg_read_ql(fg_ql_fname, tstart, tend);
end

% FFT Pareters
dt_fgm = double( median( diff(t_fgm) ) ) * 1e-9;
sr_fgm = round( 1.0 / dt_fgm );
dt_fgm = 1.0 / sr_fgm;

dt_unh = double( median( diff(t_unh) ) ) * 1e-9;
sr_unh = round( 1.0 / dt_unh );
dt_unh = 1.0 / sr_unh;

df   = 1 / T;
nfft_unh = sr_unh * T;
fN_unh   = sr_unh / 2.0;
fout_unh = fN_unh * linspace(0, 1.0, double(nfft_unh)/2.0);

nfft_fgm = sr_fgm * T;
fN_fgm   = sr_fgm / 2.0;
fout_fgm = fN_fgm * linspace(0, 1.0, double(nfft_fgm)/2.0);
	
% Tapering window
if apply_win
	win = repmat( window(@hamming, single(nfft), 'periodic')', 3, 1 );
%	win = repmat( window(@tukeywin, single(nfft), 0.75)', 3, 1 );
else
	win_fgm = ones(3, nfft_fgm);
	win_unh = ones(3, nfft_unh);
end

%------------------------------------%
% Select Subset of Data              %
%------------------------------------%

% TT2000 values of FFT interval
t0 = MrCDF_Epoch_Parse( [fft_start '.000000000'], 'cdf_time_tt2000');

% Index range of FFT interval
ifgm = find(fg_ql.tt2000 >= t0, 1, 'first');
iunh = find(t_fg         >= t0, 1, 'first');

% Subset of time-series data
t_fgm      = fg_ql.tt2000(ifgm:ifgm+nfft_fgm-1);
t_unh      = t_fg(iunh:iunh+nfft_unh-1);
b_dmpa_fgm = fg_ql.b_dmpa(1:3,ifgm:ifgm+nfft_fgm-1);
b_dmpa_unh = b_fg_dmpa(:,iunh:iunh+nfft_unh-1);

%------------------------------------%
% Take FFT                           %
%------------------------------------%
t0        = min( [t_fgm(1) t_unh(1)] );
t_sse_fgm = MrCDF_epoch2sse(t_fgm, t0);
t_sse_unh = MrCDF_epoch2sse(t_unh, t0);

% Fourier transform
Bf_dmpa_unh = fft( b_dmpa_unh .* win_unh, [], 2 );
Bf_dmpa_fgm = fft( b_dmpa_fgm .* win_fgm, [], 2 );

% Amplitude
amp_dmpa_unh = 2.0 * dt * abs( Bf_dmpa_unh(:, 1:nfft_unh/2) );
amp_dmpa_fgm = 2.0 * dt * abs( Bf_dmpa_fgm(:, 1:nfft_fgm/2) );

% Phase
phase_dmpa_unh = atan2( imag(Bf_dmpa_unh), real(Bf_dmpa_unh) ) * 180.0/pi;
phase_dmpa_fgm = atan2( imag(Bf_dmpa_fgm), real(Bf_dmpa_fgm) ) * 180.0/pi;

% Power (not including DC)
pwr_dmpa_unh = 2 * dt * abs( Bf_dmpa_unh(:, 2:nfft_unh/2).^2 );
pwr_dmpa_fgm = 2 * dt * abs( Bf_dmpa_fgm(:, 2:nfft_fgm/2).^2 );


%------------------------------------%
% Plot DMPA in Time Domain           %
%------------------------------------%
f_t = figure();

% X-component
subplot(3,1,1)
plot( t_sse_fgm, b_dmpa_fgm(1,:), ...
      t_sse_unh, b_dmpa_unh(1,:) );
title([ 'Magnetometers in DMPA '] );
xlabel( 'Time UTC' );
ylabel( {'B_{X}', '(nT)'} );
legend('B_{FGM}', 'B_{UNH}');

% Y-component
subplot(3,1,2)
plot( t_sse_fgm, b_dmpa_fgm(2,:), ...
      t_sse_unh, b_dmpa_unh(2,:) );
title([ 'Magnetometers in DMPA '] );
xlabel( 'Time UTC' );
ylabel( {'B_{Y}', '(nT)'} );

% Z-component
subplot(3,1,3)
plot( t_sse_fgm, b_dmpa_fgm(3,:), ...
      t_sse_unh, b_dmpa_unh(3,:) );
title([ 'Magnetometers in DMPA '] );
xlabel( 'Time UTC' );
ylabel( {'B_{Z}', '(nT)'} );

%------------------------------------%
% Plot DMPA in Frequency Domain      %
%------------------------------------%
fig_f = figure();
	
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
	plot( t_sse_unh, b_dmpa_unh(1, :), t_sse_fgm, b_dmpa_fgm(1, :) );
	legend('UNH', 'FGM' );
	title('FGM Calibration Comparison');
	ylabel('B (nT)');
	xlabel('Time (s)');
	
	subplot(4,1,2);
	pwr_plot( fout_unh, amp_dmpa_unh(1,:), fout_fgm, amp_dmpa_fgm(1,:) );
	legend('UNH', 'FGM');
	title('FGM Calibration Comparison');
	ylabel('X Amplitude');
	xlabel('Freuqency (Hz)');

	subplot(4,1,3);
	phs_plot( fout_unh, phase_dmpa_unh(1,1:nfft_unh/2), fout_fgm, phase_dmpa_fgm(1,1:nfft_fgm/2) );
	ylabel('X Phase');
	ylim([-200,200]);
	xlabel('Freuqency (Hz)');

	subplot(4,1,4);
	pwr_plot( fout_unh(2:end), pwr_dmpa_unh(1,:), fout_fgm(2:end), pwr_dmpa_fgm(1,:) );
	ylabel('X Power');
	xlabel('Freuqency (Hz)');
end

% Y-Component
if component == 'Y'
	subplot(4,1,1);
	plot( t_sse_unh, b_dmpa_unh(2, :), t_sse_fgm, b_dmpa_fgm(2, :) );
	legend('UNH', 'FGM' );
	title('FGM Calibration Comparison');
	ylabel('B (nT)');
	xlabel('Time (s)');
	
	subplot(4,1,2);
	pwr_plot( fout_unh, amp_dmpa_unh(2,:), fout_fgm, amp_dmpa_fgm(2,:) );
	legend('UNH', 'FGM');
	title('FGM Calibration Comparison');
	ylabel('Y Amplitude');
	xlabel('Freuqency (Hz)');

	subplot(4,1,3);
	phs_plot( fout_unh, phase_dmpa_unh(2,1:nfft_unh/2), fout_fgm, phase_dmpa_fgm(2,1:nfft_fgm/2) );
	ylabel('Y Phase');
	ylim([-200,200]);
	xlabel('Freuqency (Hz)');

	subplot(4,1,4);
	pwr_plot( fout_unh(2:end), pwr_dmpa_unh(2,:), fout_fgm(2:end), pwr_dmpa_fgm(2,:) );
	ylabel('Y Power');
	xlabel('Freuqency (Hz)');
end

% Z-Component
if component == 'Z'
	subplot(4,1,1);
	plot( t_sse_unh, b_dmpa_unh(3, :), t_sse_fgm, b_dmpa_fgm(3, :) );
	legend('UNH', 'FGM' );
	title('FGM Calibration Comparison');
	ylabel('B (nT)');
	xlabel('Time (s)');
	
	subplot(4,1,2);
	pwr_plot( fout_unh, amp_dmpa_unh(3,:), fout_fgm, amp_dmpa_fgm(3,:) );
	legend('UNH', 'FGM');
	title('FGM Calibration Comparison');
	ylabel('Z Amplitude');
	xlabel('Freuqency (Hz)');

	subplot(4,1,3);
	phs_plot( fout_unh, phase_dmpa_unh(3,1:nfft_unh/2), fout_fgm, phase_dmpa_fgm(3,1:nfft_fgm/2) );
	ylabel('Z Phase');
	ylim([-200,200]);
	xlabel('Freuqency (Hz)');

	subplot(4,1,4);
	pwr_plot( fout_unh(2:end), pwr_dmpa_unh(3,:), fout_fgm(2:end), pwr_dmpa_fgm(3,:) );
	ylabel('Z Power');
	xlabel('Freuqency (Hz)');
end

%------------------------------------%
% Amplitude, Phase, & Power Diff     %
%------------------------------------%
fig_f = figure();
	
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
	plot( t_sse_unh, b_dmpa_unh(1, :), t_sse_fgm, b_dmpa_fgm(1, :) );
	legend('UNH', 'FGM' );
	title('FGM Calibration Comparison');
	ylabel('B (nT)');
	xlabel('Time (s)');
	
	subplot(4,1,2);
	pwr_plot( fout_fgm, abs( amp_dmpa_unh(1,1:nfft_fgm/2) - amp_dmpa_fgm(1,:) ) );
	legend('UNH - FGM');
	title('FGM Calibration Comparison');
	ylabel('X Amplitude');
	xlabel('Freuqency (Hz)');

	subplot(4,1,3);
	phs_plot( fout_fgm, phase_dmpa_unh(1,1:nfft_fgm/2) - phase_dmpa_fgm(1,1:nfft_fgm/2) );
	ylabel('X Phase');
	ylim([-200,200]);
	xlabel('Freuqency (Hz)');

	subplot(4,1,4);
	pwr_plot( fout_fgm(2:end), abs( pwr_dmpa_unh(1,2:nfft_fgm/2) - pwr_dmpa_fgm(1,:) ) );
	ylabel('X Power');
	xlabel('Freuqency (Hz)');
end

% Y-Component
if component == 'Y'
	subplot(4,1,1);
	plot( t_sse_unh, b_dmpa_unh(2, :), t_sse_fgm, b_dmpa_fgm(2, :) );
	legend('UNH', 'FGM' );
	title('FGM Calibration Comparison');
	ylabel('B (nT)');
	xlabel('Time (s)');
	
	subplot(4,1,2);
	pwr_plot( fout_fgm, abs( amp_dmpa_unh(2,1:nfft_fgm/2) - amp_dmpa_fgm(2,:) ) );
	legend('UNH - FGM' );
	title('FGM Calibration Comparison');
	ylabel('Y Amplitude');
	xlabel('Freuqency (Hz)');

	subplot(4,1,3);
	phs_plot( fout_fgm, phase_dmpa_unh(2,1:nfft_fgm/2) - phase_dmpa_fgm(2,1:nfft_fgm/2) );
	ylabel('Y Phase');
	ylim([-200,200]);
	xlabel('Freuqency (Hz)');

	subplot(4,1,4);
	pwr_plot( fout_fgm(2:end), abs( pwr_dmpa_unh(2,nfft_fgm/2) - pwr_dmpa_fgm(2,:) ) );
	ylabel('Y Power');
	xlabel('Freuqency (Hz)');
end

% Z-Component
if component == 'Z'
	subplot(4,1,1);
	plot( t_sse_unh, b_dmpa_unh(3, :), t_sse_fgm, b_dmpa_fgm(3, :) );
	legend('UNH', 'FGM' );
	title('FGM Calibration Comparison');
	ylabel('B (nT)');
	xlabel('Time (s)');

	subplot(4,1,2);
	pwr_plot( fout_fgm, abs( amp_dmpa_unh(3,1:nfft_fgm/2) - amp_dmpa_fgm(3,:) ) );
	legend('UNH - FGM');
	title('FGM Calibration Comparison');
	ylabel('Z Amplitude');
	xlabel('Freuqency (Hz)');

	subplot(4,1,3);
	phs_plot( fout_fgm, phase_dmpa_unh(3,1:nfft_fgm/2) - phase_dmpa_fgm(3,1:nfft_fgm/2) );
	ylabel('Z Phase');
	ylim([-200,200]);
	xlabel('Freuqency (Hz)');

	subplot(4,1,4);
	pwr_plot( fout_fgm(2:end), abs( pwr_dmpa_unh(3,2:nfft_fgm/2) - pwr_dmpa_fgm(3,:) ) );
	ylabel('Z Power');
	xlabel('Freuqency (Hz)');
end
