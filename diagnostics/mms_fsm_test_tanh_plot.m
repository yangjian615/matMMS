%
% Name
%   mms_fsm_test_savgol
%
% Purpose
%   
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2016-09-09      Written by Matthew Argall
%
%***************************************************************************

% Load data from mat file
instr   = 'dfg';
iComp   = 3;
png_dir = '/home/argall/figures/merging/'; % '/home/argall/figures/merging/';
f0      = 4.0;
file    = '/home/argall/mms4_fsm_brst_l3_cal_20151101000000.mat';
file    = '/home/argall/mms4_fsm_brst_l2plus_cal-floor-month_20151101000000_v0.0.0.mat';
load( file );

% Amount of binned data
[sc, ~, ~, ~, ~, ~, optdesc] = mms_dissect_filename(file);
interval = regexp(optdesc, '-', 'split');
interval = interval{3};

% LPP Noise Floor
lpp_noise = mms_lpp_read_floor();
f_lpp     = lpp_noise.([sc '_f_nemi']);
nf_lpp    = lpp_noise.([sc '_b_nemi'])(iComp,:).^2.0;


% Data
if strcmp(instr, 'dfg')
	f           = dfg_f;
	bins        = dfg_bins;
	hfit        = dfg_fit;
	dist        = dfg_dist;
	dist_smooth = dfg_dist_smooth;
	dist_fit    = dfg_dist_fit;
	step        = dfg_step;
	step_smooth = dfg_step_smooth;
	step_fit    = dfg_step_fit;
	nf          = dfg_floor;
else
	f           = scm_f;
	bins        = scm_bins;
	hfit        = scm_fit;
	dist        = scm_dist;
	dist_smooth = scm_dist_smooth;
	dist_fit    = scm_dist_fit;
	step        = scm_step;
	step_smooth = scm_step_smooth;
	step_fit    = scm_step_fit;
	nf          = scm_floor;
end

% Clear the data
clear -regex 'scm_*' 'dfg_*'


%------------------------------------%
% Scale Data                         %
%------------------------------------%

% Scale the PSD
%sclPSD = MrRescale( log10( psd ), 1, 64, ...
%                    'MinValue', -6, ...
%                    'MaxValue', -1, ...
%                    'Class',    'uint8' );

%
% Histogram
%
sclDist = MrRescale( dist, 1, 64, ...
                     'Class', 'uint8' );

sclDist_smooth = MrRescale( dist_smooth, 1, 64, ...
                            'Class', 'uint8' );


dist_fit( dist_fit > 1.5*max(dist_smooth(:))  ) = 0;
dist_fit( dist_fit < 0 ) = 0;
sclDist_fit = MrRescale( dist_fit, 1, 64, ...
                         'Class', 'uint8' );

%
% Step response
%
sclStep_hist = MrRescale( step, 1, 64, ...
                          'Class', 'uint8' );

sclStep_smooth = MrRescale( step_smooth, 1, 64, ...
                            'Class', 'uint8' );


step_fit( step_fit > 1 ) = 0;
step_fit( step_fit < 0 ) = 0;
sclStep_fit = MrRescale( step_fit, 1, 64, ...
                         'Class', 'uint8' );

%------------------------------------%
% Plot Distributions                 %
%------------------------------------%

set(0,'defaultLineLineWidth', 1.5)

% Create a figure
fig1 = figure( 'OuterPosition', [10, 10, 800, 650] );

layout = [3,2];
[inPos, outPos] = MrLayout( layout,              ...
                            'Figure',   fig1,    ...
                            'OXMargin', [4,1],   ...
                            'OYMargin', [2,0],   ...
                            'IYMargin', [1.5,2], ...
                            'XGap',     6,       ...
                            'YGap',     0 );

% Convert back to log-scale
bins = 10.0.^bins;
nf   = 10.0.^nf;

% Original Histogram
idx = 1;
subplot(layout(1), layout(2), idx, 'OuterPosition', outPos(idx,:), 'Position', inPos(idx,:));
h = pcolor( f, bins, sclDist );
h.EdgeColor = 'none';
ax = gca();
ax.XScale = 'log';
ax.YScale = 'log';
colorbar
ax.Position = inPos(idx,:);
hold on
h = plot( f, nf, '--m' );
h.LineWidth = 2.0;
h = plot( f_lpp, nf_lpp, '--c' );
h.LineWidth = 2.0;
line([f0, f0], [bins(1) bins(end)], 'LineStyle', '--', 'Color', [0,0,0])
hold off
ylabel({'PSD', 'nT^2/Hz'})
title( ['Noise Distribution: ' upper(instr)] )


% Smoothed Histogram
idx = 3;
subplot(layout(1), layout(2), idx, 'OuterPosition', outPos(idx,:), 'Position', inPos(idx,:));
h = pcolor( f, bins, sclDist_smooth );
h.EdgeColor = 'none';
ax = gca();
ax.XScale = 'log';
ax.YScale = 'log';
colorbar
ax.Position = inPos(idx,:);
hold on
h = plot(f, nf, '--m');
h.LineWidth = 2.0;
h = plot( f_lpp, nf_lpp, '--c' );
h.LineWidth = 2.0;
hold off
ylabel({'PSD', 'nT^2/Hz'})
title('Smoothed')


% Fitted Histogram
idx = 5;
subplot(layout(1), layout(2), idx, 'OuterPosition', outPos(idx,:), 'Position', inPos(idx,:));
h = pcolor( f, bins, sclDist_fit );
h.EdgeColor = 'none';
ax = gca();
ax.XScale = 'log';
ax.YScale = 'log';
colorbar
ax.Position = inPos(idx,:);
hold on
h = plot(f, nf, '--m');
h = plot( f_lpp, nf_lpp, '--c' );
hold off
xlabel('f (Hz)')
ylabel({'PSD', 'nT^2/Hz'})
title('Fitted')


% Original Step Response
idx = 2;
subplot(layout(1), layout(2), idx, 'OuterPosition', outPos(idx,:), 'Position', inPos(idx,:));
h = pcolor( f, bins, sclStep_hist );
h.EdgeColor = 'none';
ax = gca();
ax.XScale = 'log';
ax.YScale = 'log';
hold on
h = plot(f, nf, '--m');
hold off
title('Probability Distribution')


% Smoothed Step Response
idx = 4;
subplot(layout(1), layout(2), idx, 'OuterPosition', outPos(idx,:), 'Position', inPos(idx,:));
h = pcolor( f, bins, sclStep_smooth );
h.EdgeColor = 'none';
ax = gca();
ax.XScale = 'log';
ax.YScale = 'log';
hold on
h = plot(f, nf, '--m');
hold off
title('Smoothed')


% Fitted Step Response
idx = 6;
subplot(layout(1), layout(2), idx, 'OuterPosition', outPos(idx,:), 'Position', inPos(idx,:));
h = pcolor( f, bins, sclStep_fit );
h.EdgeColor = 'none';
ax = gca();
ax.XScale = 'log';
ax.YScale = 'log';
hold on
h = plot(f, nf, '--m');
hold off
xlabel('f (Hz)')
title('TanH Fit')


%------------------------------------%
% Plot A Cut                         %
%------------------------------------%
if0 = find(f >= f0, 1, 'first');

fig2 = figure( 'OuterPosition', [50, 10, 800, 650] );

% Distribution
subplot(2,1,1);
h = plot( bins, dist(:,if0), bins, dist_smooth(:,if0), bins, dist_fit(:,if0) );
ax = gca();
ax.XScale = 'log';
legend('Signal', 'Smooth', 'TanH');
xlabel('PSD (nT^2/Hz)')
ylabel('Occurrence')
title(['Distribution at f=' num2str(f0) 'Hz'])



% Probability distribution
subplot(2,1,2);
h = plot( bins, step(:,if0), bins, step_smooth(:,if0), bins, step_fit(:,if0) );
ax = gca();
ax.XScale = 'log';
legend('Signal', 'Smooth', 'TanH');
xlabel('PSD (nT^2/Hz)')
ylabel('Probability')
ylim([-0.5,1.5])
title(['Probability Distribution at f=' num2str(f0) 'Hz'])



%------------------------------------%
% Save Figure                        %
%------------------------------------%
if ~isempty(png_dir)
	% Distributions
	set( fig1, 'PaperPositionMode', 'auto');
	filename = fullfile( png_dir, ['mms_fsm_tanh_dist-' instr '-' interval '.png'] );
	print( fig1, filename, '-dpng', '-r300' );
	disp( ['Saving file to: "' filename '".'] );
	
	% Cuts
	set( fig2, 'PaperPositionMode', 'auto');
	filename = fullfile( png_dir, ['mms_fsm_tanh_dist-' instr '-' interval '_f' strrep( num2str(f0), '.', 'p' ) '.png'] );
	print( fig2, filename, '-dpng', '-r300' );
	disp( ['Saving file to: "' filename '".'] );
end
