%
% Name
%   mms_fsm_test_tanh_filter
%
% Purpose
%   Display the impulse, step, and frequency response of the filter used to
%   merge FSM data.
%
% Calling Sequence
%   [] = MrFilter_WinSinc()
%     Display the effects of windowing on a sinc function.
%
% Parameters
%
% Returns
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2017-02-08      Written by Matthew Argall
%

f_low  = 1.0;
f_high = 64.0;

%------------------------------------%
% Load Data                          %
%------------------------------------%

% Load data from mat file
file    = '/home/argall/mms4_fsm_brst_l3_cal_20151101000000.mat';
file    = 'mms4_fsm_brst_l2plus_cal-dfg-month_20151101000000_v0.0.0.mat';
file    = 'mms4_fsm_brst_l2plus_cal-floor-month_20151101000000_v0.0.0.mat';
png_dir = '/home/argall/figures/merging/'; % '/home/argall/figures/merging/'
load( file );

[~, ~, ~, ~, ~, ~, optdesc] = mms_dissect_filename(file);
interval = regexp(optdesc, '-', 'split');
interval = interval{3};

% Pick relevant data and clear the rest.
f_dfg    = dfg_f;
f_scm    = scm_f;
bins_dfg = dfg_bins;
bins_scm = scm_bins;
fit_dfg  = dfg_floor;
fit_scm  = scm_floor;
clear -regex 'scm_*' 'dfg_*'

%------------------------------------%
% Filter Range                       %
%------------------------------------%

% Size of dimensions
nFreq = length(f_scm);

% Index range
if0_dfg = find(f_dfg >= f_low,  1, 'first');
if0_scm = find(f_scm >= f_low,  1, 'first');
if1_dfg = find(f_dfg <= f_high, 1, 'last');
if1_scm = find(f_scm <= f_high, 1, 'last');

% DFG Nyquist
if_end = length(f_dfg);

%------------------------------------%
% Create Weight Filter               %
%------------------------------------%

% Start by smoothing the fit
fit_dfg_smooth = MrSmooth( fit_dfg, 11, 'truncate' );
fit_scm_smooth = MrSmooth( fit_scm, 11, 'truncate' );

% Weighting factor
%   - Invert noise floor
temp_dfg = squeeze( 1.0 ./ fit_dfg_smooth(if0_dfg:end)    );
temp_scm = squeeze( 1.0 ./ fit_scm_smooth(if0_scm:if_end) );

% DFG Weight function
w_dfg                  = zeros(1,nFreq);
w_dfg(1:if0_scm-1)     = 1.0;
w_dfg(if0_scm:if1_scm) = temp_scm ./ (temp_dfg + temp_scm);
w_dfg(if1_scm+1:end)   = 0.0;


% SCM Weight Function
w_scm = 1.0 - w_dfg;

%------------------------------------%
% Impulse Response                   %
%------------------------------------%

% Reflect to negative frequencies
w_ir_dfg = ifft( [w_dfg, fliplr( w_dfg(2:end) )] );
w_ir_scm = ifft( [w_scm, fliplr( w_scm(2:end) )] );

% Trim
% w_ir_dfg = w_ir_dfg(1:nFreq);
% w_ir_scm = w_ir_scm(1:nFreq);

% Time
time     = 1:2*nFreq-1;
w_ir_dfg = [ w_ir_dfg(nFreq+1:end) w_ir_dfg(1:nFreq) ];

%------------------------------------%
% Windowed Impulse Response          %
%------------------------------------%
N = length(time); % Number of points in the filter
M = floor(N/6);
if mod(M, 2) ~= 0
	M = M-1;      % Width of the square impulse
end

% Blackman & Hamming windows
blackwin = 0.42 - 0.50 * cos(2*pi*(0:M)/M) + 0.08 * cos(4*pi*(0:M)/M);
hamwin   = 0.54 - 0.46 * cos(2*pi*(0:M)/M);

% Pad windows with 0's
blackwin = [ zeros(1,nFreq-M/2-1) blackwin zeros(1,nFreq-M/2-1) ];
hamwin   = [ zeros(1,nFreq-M/2-1)  hamwin  zeros(1,nFreq-M/2-1) ];

% Apply Blackman window
w_ir_win_dfg = w_ir_dfg .* blackwin;

%------------------------------------%
% Windowed Frequency Response        %
%------------------------------------%

w_fr_win_dfg = fft(w_ir_win_dfg);
w_fr_win_dfg = w_fr_win_dfg(1:nFreq);

amp_fr_win_dfg          = 0.5 * abs(w_fr_win_dfg);
amp_fr_win_dfg(2:end-1) = 2.0 * amp_fr_win_dfg(2:end-1);

%------------------------------------%
% Frequency Response (Original)      %
%------------------------------------%
fig = figure( 'OuterPosition', [10, 10, 800, 650] );

subplot(3,1,1);
h = loglog(f_dfg, 10.0.^fit_dfg, f_dfg, 10.0.^fit_dfg_smooth, f_scm, 10.0.^fit_scm);
legend('dfg', 'dfg-smooth', 'scm');
xlim([0.9, 64]);
ylabel({'PSD', 'nT^2/Hz'});
title('Noise Floor')


subplot(3,1,2);
h = semilogx(f_scm, w_dfg, f_scm, amp_fr_win_dfg);
h(1).LineWidth = 2.0;
h(2).LineWidth = 2.0;
line( get(gca, 'XLim'), [0.5, 0.5], 'LineStyle', '--' );
legend('original', 'windowed');
xlabel('f (Hz)');
xlim([0.9, 64]);
ylim([0.4, 0.6]);
ylabel('Filter');
title( ['Frequency Response: N=' num2str(N) ' M=' num2str(M)] )


subplot(3,1,3);
h = plot(time-nFreq, w_ir_dfg, time-nFreq, blackwin*max(w_ir_dfg), time-nFreq, w_ir_win_dfg);
legend('original', 'windowed');
xlabel('Sample Number');
xlim([-M/2,+M/2]);
ylabel('Amplitude');
title( ['Impulse Response: N=' num2str(N) ' M=' num2str(M)] )

%------------------------------------%
% Save Figure                        %
%------------------------------------%
if ~isempty(png_dir)
	set(gcf, 'PaperPositionMode', 'auto');
	filename = fullfile( png_dir, ['mms_fsm_tanh_filter_' interval '_n' num2str(N) '_m' num2str(M) '.png'] );
	print( filename, '-dpng', '-r300' );
	disp( ['Saving file to: "' filename '".'] );
end