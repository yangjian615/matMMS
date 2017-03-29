%
% Name
%   mms_fsm_plot_cal_psd_hist
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
fgm_file    = '/nfs/fsm/temp/mms1_fsm_brst_l2plus_cal-dfg_20150901121114_v2.0.0.cdf';
scm_file    = '/nfs/fsm/temp/mms1_fsm_brst_l2plus_cal-scm_20150901121114_v2.0.0.cdf';
theFlag_fgm = 5;
theFlag_scm = 1;
png_dir     = '/home/argall/figures/merging/'; % '/home/argall/figures/merging/';


% Amount of binned data
[sc, instr, mode, level, tstart, ~, optdesc] = mms_dissect_filename(fgm_file);
fgm_instr = regexp(optdesc, '-', 'split');
fgm_instr = fgm_instr{2};

% LPP Noise Floor
lpp_noise = mms_lpp_read_floor();
f_lpp     = lpp_noise.([sc '_f_nemi']);
nf_lpp    = lpp_noise.([sc '_b_nemi']).^2.0;

% IWF Noise Floor
iwf_noise = mms_iwf_read_floor();
f_iwf     = iwf_noise.f;
nf_iwf    = [ iwf_noise.Bx; iwf_noise.By; iwf_noise.Bz ].^2.0;

%------------------------------------%
% Read Data                          %
%------------------------------------%

% Variable names
psd_vname  = [sc '_fsm_psd_omb_'  mode '_' level];
hist_vname = [sc '_fsm_psd_hist_' mode '_' level];

% FGM
[psd_fgm, t_fgm, f_fgm]           = MrCDF_Read( fgm_file, psd_vname  );
[hist_fgm, ~, bins_fgm, flag_fgm] = MrCDF_Read( fgm_file, hist_vname );

% SCM
[psd_scm, t_scm, f_scm]           = MrCDF_Read( scm_file, psd_vname  );
[hist_scm, ~, bins_scm, flag_scm] = MrCDF_Read( scm_file, hist_vname );


% DEC 32/64 or ADC A/B
fgm_title = '';
if strcmp(fgm_instr, 'dfg')
	if bitget(theFlag_fgm, 3)
		fgm_title = [fgm_title ' DEC32'];
	else
		fgm_title = [fgm_title ' DEC64'];
	end
else
	if bitget(theFlag_fgm, 3)
		fgm_title = [fgm_title ' ADCB'];
	else
		fgm_title = [fgm_title ' ADCA'];
	end
end

%------------------------------------%
% Select Flag                        %
%------------------------------------%
	% Flag index
	iflag_fgm = find(flag_fgm == theFlag_fgm);
	iflag_scm = find(flag_scm == theFlag_scm);
	assert( ~isempty(iflag_fgm), [ 'No such flag for FGM (' num2str(theFlag_fgm) ').' ] );
	assert( ~isempty(iflag_scm), [ 'No such flag for SCM (' num2str(theFlag_scm) ').' ] );
	
	% Trim the data
	psd_fgm  = log10( permute( psd_fgm, [2,1,3] ) );
	psd_scm  = log10( permute( psd_scm, [2,1,3] ) );
	hist_fgm = permute( squeeze( hist_fgm(:, :, iflag_fgm, :) ), [2,1,3] );
	hist_scm = permute( squeeze( hist_scm(:, :, iflag_scm, :) ), [2,1,3] );
	
%------------------------------------%
% Plot PSD                           %
%------------------------------------%

% Create a figure
fig1 = figure( 'OuterPosition', [10, 10, 1010, 660] );

layout = [3,2];
[inPos, outPos] = MrLayout( layout,              ...
                            'Figure',   fig1,    ...
                            'OXMargin', [4,16],   ...
                            'OYMargin', [2,0],   ...
                            'IYMargin', [1.5,2], ...
                            'XGap',     18,       ...
                            'YGap',     0 );

% Convert back to log-scale
bins_fgm = 10.0.^bins_fgm;
bins_scm = 10.0.^bins_scm;

% Convert to seconds
t_dn_fgm = MrCDF_epoch2datenum(t_fgm);
t_dn_scm = MrCDF_epoch2datenum(t_scm);

% FGM BX
idx = 1;
subplot(layout(1), layout(2), idx, 'OuterPosition', outPos(idx,:), 'Position', inPos(idx,:));
h = pcolor( t_dn_fgm, f_fgm, psd_fgm(:,:,1) );
h.EdgeColor = 'none';
ax = gca();
ax.XLim          = [t_dn_fgm(1) t_dn_fgm(end)];
ax.YLabel.String = 'f (Hz)';
ax.YLim          = [0.5, f_fgm(end)];
ax.YScale        = 'log';
h = colorbar();
h.YLabel.String  = 'PSD (nT^{2}/Hz)';
ax.Position      = inPos(idx,:);
title( ['PSD Bx ' upper(fgm_instr) fgm_title] );

% FGM BY
idx = 3;
subplot(layout(1), layout(2), idx, 'OuterPosition', outPos(idx,:), 'Position', inPos(idx,:));
h = pcolor( t_dn_fgm, f_fgm, psd_fgm(:,:,3) );
h.EdgeColor = 'none';
ax = gca();
ax.XLim          = [t_dn_fgm(1) t_dn_fgm(end)];
ax.YLabel.String = 'f (Hz)';
ax.YLim          = [0.5, f_fgm(end)];
ax.YScale        = 'log';
h = colorbar();
h.YLabel.String  = 'PSD (nT^{2}/Hz)';
ax.Position      = inPos(idx,:);
title( ['PSD By ' upper(fgm_instr) fgm_title] )

% FGM BZ
idx = 5;
subplot(layout(1), layout(2), idx, 'OuterPosition', outPos(idx,:), 'Position', inPos(idx,:));
h = pcolor( t_dn_fgm, f_fgm, psd_fgm(:,:,2) );
h.EdgeColor = 'none';
ax = gca();
ax.XLabel.String = 'Time';
ax.XLim          = [t_dn_fgm(1) t_dn_fgm(end)];
ax.YLabel.String = 'f (Hz)';
ax.YLim          = [0.5, f_fgm(end)];
ax.YScale        = 'log';
h = colorbar();
h.YLabel.String  = 'PSD (nT^{2}/Hz)';
ax.Position      = inPos(idx,:);
title( ['PSD Bz ' upper(fgm_instr) fgm_title] );


% SCM BX
idx = 2;
subplot(layout(1), layout(2), idx, 'OuterPosition', outPos(idx,:), 'Position', inPos(idx,:));
h = pcolor( t_dn_scm, f_scm, psd_scm(:,:,1) );
h.EdgeColor = 'none';
ax = gca();
ax.XLim          = [t_dn_fgm(1) t_dn_fgm(end)];
ax.YLabel.String = 'f (Hz)';
ax.YLim          = [0.5, f_fgm(end)];
ax.YScale        = 'log';
h = colorbar();
h.YLabel.String  = 'PSD (nT^{2}/Hz)';
ax.Position      = inPos(idx,:);
title( 'PSD Bx SCM' );

% SCM BY
idx = 4;
subplot(layout(1), layout(2), idx, 'OuterPosition', outPos(idx,:), 'Position', inPos(idx,:));
h = pcolor( t_dn_scm, f_scm, psd_scm(:,:,3) );
h.EdgeColor = 'none';
ax = gca();
ax.XLim          = [t_dn_fgm(1) t_dn_fgm(end)];
ax.YLabel.String = 'f (Hz)';
ax.YLim          = [0.5, f_fgm(end)];
ax.YScale        = 'log';
h = colorbar();
h.YLabel.String  = 'PSD (nT^{2}/Hz)';
ax.Position      = inPos(idx,:);
title( 'PSD By SCM' );

% SCM BZ
idx = 6;
subplot(layout(1), layout(2), idx, 'OuterPosition', outPos(idx,:), 'Position', inPos(idx,:));
h = pcolor( t_dn_scm, f_scm, psd_scm(:,:,2) );
h.EdgeColor = 'none';
ax = gca();
ax.XLabel.String = 'Time';
ax.XLim          = [t_dn_fgm(1) t_dn_fgm(end)];
ax.YLabel.String = 'f (Hz)';
ax.YLim          = [0.5, f_fgm(end)];
ax.YScale        = 'log';
h = colorbar();
h.YLabel.String  = 'PSD (nT^{2}/Hz)';
ax.Position      = inPos(idx,:);
title( 'PSD Bz SCM' );

%------------------------------------%
% Plot Histogram                     %
%------------------------------------%

% Create a figure
fig2 = figure( 'OuterPosition', [10, 10, 1010, 660] + 50 );

layout = [3,2];
[inPos, outPos] = MrLayout( layout,              ...
                            'Figure',   fig2,    ...
                            'OXMargin', [4,12],   ...
                            'OYMargin', [2,0],   ...
                            'IYMargin', [1.5,2], ...
                            'XGap',     18,       ...
                            'YGap',     0 );

% FGM BX
idx = 1;
subplot(layout(1), layout(2), idx, 'OuterPosition', outPos(idx,:), 'Position', inPos(idx,:));
h = pcolor( f_fgm, bins_fgm, hist_fgm(:,:,1) );
h.EdgeColor = 'none';
ax = gca();
ax.XScale        = 'log';
ax.XLim          = [0.5, f_fgm(end)];
ax.YLabel.String = {'PSD', 'nT^2/Hz'};
ax.YScale        = 'log';
h = colorbar();
h.YLabel.String  = 'Occurrence';;
ax.Position      = inPos(idx,:);
title( ['PSD Bx ' upper(fgm_instr) fgm_title] );

hold on
plot( f_iwf, nf_iwf(1,:), '--r', 'LineWidth', 1.5 );
hold off


% FGM BY
idx = 3;
subplot(layout(1), layout(2), idx, 'OuterPosition', outPos(idx,:), 'Position', inPos(idx,:));
h = pcolor( f_fgm, bins_fgm, hist_fgm(:,:,2) );
h.EdgeColor = 'none';
ax = gca();
ax.XScale        = 'log';
ax.XLim          = [0.5, f_fgm(end)];
ax.YLabel.String = {'PSD', 'nT^2/Hz'};
ax.YScale        = 'log';
h = colorbar();
h.YLabel.String  = 'Occurrence';
ax.Position      = inPos(idx,:);
title( ['PSD By ' upper(fgm_instr) fgm_title] );

hold on
plot( f_iwf, nf_iwf(2,:), '--r', 'LineWidth', 1.5 );
hold off


% FGM BZ
idx = 5;
subplot(layout(1), layout(2), idx, 'OuterPosition', outPos(idx,:), 'Position', inPos(idx,:));
h = pcolor( f_fgm, bins_fgm, hist_fgm(:,:,3) );
h.EdgeColor = 'none';
ax = gca();
ax.XScale        = 'log';
ax.XLabel.String = 'f (Hz)';
ax.XLim          = [0.5, f_fgm(end)];
ax.YLabel.String = {'PSD', 'nT^2/Hz'};
ax.YScale        = 'log';
h = colorbar();
h.YLabel.String  = 'Occurrence';
ax.Position      = inPos(idx,:);
title( ['PSD Bz ' upper(fgm_instr) fgm_title] );

hold on
plot( f_iwf, nf_iwf(3,:), '--r', 'LineWidth', 1.5 );
hold off


% SCM BX
idx = 2;
subplot(layout(1), layout(2), idx, 'OuterPosition', outPos(idx,:), 'Position', inPos(idx,:));
h = pcolor( f_scm, bins_scm, hist_scm(:,:,1) );
h.EdgeColor = 'none';
ax = gca();
ax.XScale        = 'log';
ax.XLim          = [1.0, f_scm(end)];
ax.YLabel.String = {'PSD', 'nT^2/Hz'};
ax.YScale        = 'log';
h = colorbar();
h.YLabel.String  = 'Occurrence';
ax.Position      = inPos(idx,:);
title( 'PSD Bx SCM' )

hold on
plot( f_lpp, nf_lpp(1,:), '--r', 'LineWidth', 1.5 );
hold off


% SCM BY
idx = 4;
subplot(layout(1), layout(2), idx, 'OuterPosition', outPos(idx,:), 'Position', inPos(idx,:));
h = pcolor( f_scm, bins_scm, hist_scm(:,:,2) );
h.EdgeColor = 'none';
ax = gca();
ax.XScale        = 'log';
ax.XLim          = [1.0, f_scm(end)];
ax.YLabel.String = {'PSD', 'nT^2/Hz'};
ax.YScale        = 'log';
h = colorbar();
h.YLabel.String  = 'Occurrence';
ax.Position      = inPos(idx,:);
title( 'PSD By SCM' )

hold on
plot( f_lpp, nf_lpp(2,:), '--r', 'LineWidth', 1.5 );
hold off


% SCM BZ
idx = 6;
subplot(layout(1), layout(2), idx, 'OuterPosition', outPos(idx,:), 'Position', inPos(idx,:));
h = pcolor( f_scm, bins_scm, hist_scm(:,:,3) );
h.EdgeColor = 'none';
ax = gca();
ax.XScale        = 'log';
ax.XLabel.String = 'f (Hz)';
ax.XLim          = [1.0, f_scm(end)];
ax.YLabel.String = {'PSD', 'nT^2/Hz'};
ax.YScale        = 'log';
h = colorbar();
h.YLabel.String  = 'Occurrence';
ax.Position      = inPos(idx,:);
title( 'PSD Bz SCM' )

hold on
plot( f_lpp, nf_lpp(3,:), '--r', 'LineWidth', 1.5 );
hold off


%------------------------------------%
% Save Figure                        %
%------------------------------------%
if ~isempty(png_dir)
	% Distributions
	set( fig1, 'PaperPositionMode', 'auto');
	filename = fullfile( png_dir, [sc '_fsm_' mode '_' level '_' optdesc '-psd_' tstart '.png'] );
	print( fig1, filename, '-dpng', '-r300' );
	disp( ['Saving file to: "' filename '".'] );
	
	% Cuts
	set( fig2, 'PaperPositionMode', 'auto');
	filename = fullfile( png_dir, [sc '_fsm_' mode '_' level '_' optdesc '-hist_' tstart '.png'] );
	print( fig2, filename, '-dpng', '-r300' );
	disp( ['Saving file to: "' filename '".'] );
end
