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
file    = '/nfs/fsm/temp/mms1_fsm_brst_l2plus_cal-dfg-scm_20150901121114_v0.0.0.cdf';
theFlag = 5;
png_dir = '/home/argall/figures/merging/'; % '/home/argall/figures/merging/';


% Amount of binned data
[sc, instr, mode, level, tstart, ~, optdesc] = mms_dissect_filename(file);
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
gain_vname       = [sc '_fsm_psdrat_omb_'      mode '_' level];
gain_hist_vname  = [sc '_fsm_psdrat_hist_'     mode '_' level];
phase_vname      = [sc '_fsm_phaseshift_omb_'  mode '_' level];
phase_hist_vname = [sc '_fsm_phaseshift_hist_' mode '_' level];

% GAIN
[gain, t, f]                    = MrCDF_Read( file, gain_vname  );
[gain_hist, ~, gain_bins, flag] = MrCDF_Read( file, gain_hist_vname );

% PHASE
phase                       = MrCDF_Read( file, phase_vname  );
[phase_hist, ~, phase_bins] = MrCDF_Read( file, phase_hist_vname );


% DEC 32/64 or ADC A/B
fgm_title = '';
if strcmp(fgm_instr, 'dfg')
	if bitget(theFlag, 3)
		fgm_title = [fgm_title ' DEC32'];
	else
		fgm_title = [fgm_title ' DEC64'];
	end
else
	if bitget(theFlag, 3)
		fgm_title = [fgm_title ' ADCB'];
	else
		fgm_title = [fgm_title ' ADCA'];
	end
end

%------------------------------------%
% Select Flag                        %
%------------------------------------%
switch mode
	case 'brst'
		frange = [1.0, 64.0];
	case {'fast', 'srvy'}
		frange = [0.5, 16.0];
	case 'slow'
		frange = [0.5, 4.0];
	otherwise
		error( ['Invalid mode: "' mode '".'] );
end

% Trim frequency range for quicker plotting
ifrange    = find( f >= frange(1) & f <= frange(2) );
f          = f(ifrange);
gain       = gain(:,ifrange,:);
phase      = phase(:,ifrange,:);
gain_hist  = gain_hist(ifrange,:,:,:);
phase_hist = phase_hist(ifrange,:,:,:);


% Flag index
iflag = find(flag == theFlag);
assert( ~isempty(iflag), [ 'No such flag (' num2str(theFlag) ').' ] );

% Trim the data
gain       = log10( permute( gain, [2,1,3] ) );
phase      = permute( phase, [2,1,3] );
gain_hist  = permute( squeeze( gain_hist(:, :, iflag, :) ), [2,1,3] );
phase_hist = permute( squeeze( phase_hist(:, :, iflag, :) ), [2,1,3] );
	
%------------------------------------%
% Plot PSD                           %
%------------------------------------%

% Create a figure
fig1 = figure( 'OuterPosition', [10, 10, 1010, 660] );

layout = [3,2];
[inPos, outPos] = MrLayout( layout,              ...
                            'Figure',   fig1,    ...
                            'OXMargin', [4,12],   ...
                            'OYMargin', [2,0],   ...
                            'IYMargin', [1.5,2], ...
                            'XGap',     18,       ...
                            'YGap',     0 );

% Convert back to log-scale
gain_bins = 10.0.^gain_bins;

% Convert to seconds
t_dn = MrCDF_epoch2datenum(t);

% GAIN BX
idx = 1;
subplot(layout(1), layout(2), idx, 'OuterPosition', outPos(idx,:), 'Position', inPos(idx,:));
h = pcolor( t_dn, f, gain(:,:,1) );
h.EdgeColor = 'none';
ax = gca();
ax.XLim          = [t_dn(1) t_dn(end)];
ax.YLabel.String = 'f (Hz)';
ax.YLim          = frange;
ax.YScale        = 'log';
h = colorbar();
h.YLabel.String  = 'Log_{10}(Gain)';
ax.Position      = inPos(idx,:);
title( ['Gain Bx ' upper(fgm_instr) '/SCM' fgm_title] );

% GAIN BY
idx = 3;
subplot(layout(1), layout(2), idx, 'OuterPosition', outPos(idx,:), 'Position', inPos(idx,:));
h = pcolor( t_dn, f, gain(:,:,2) );
h.EdgeColor = 'none';
ax = gca();
ax.XLim          = [t_dn(1) t_dn(end)];
ax.YLabel.String = 'f (Hz)';
ax.YLim          = frange;
ax.YScale        = 'log';
h = colorbar();
h.YLabel.String  = 'Log_{10}(Gain)';
ax.Position      = inPos(idx,:);
title( ['Gain By ' upper(fgm_instr) '/SCM' fgm_title] )

% GAIN BZ
idx = 5;
subplot(layout(1), layout(2), idx, 'OuterPosition', outPos(idx,:), 'Position', inPos(idx,:));
h = pcolor( t_dn, f, gain(:,:,3) );
h.EdgeColor = 'none';
ax = gca();
ax.XLabel.String = 'Time';
ax.XLim          = [t_dn(1) t_dn(end)];
ax.YLabel.String = 'f (Hz)';
ax.YLim          = frange;
ax.YScale        = 'log';
h = colorbar();
h.YLabel.String  = 'Log_{10}(Gain)';
ax.Position      = inPos(idx,:);
title( ['Gain Bz ' upper(fgm_instr) '/SCM' fgm_title] );

% Gain Hist BX
idx = 2;
subplot(layout(1), layout(2), idx, 'OuterPosition', outPos(idx,:), 'Position', inPos(idx,:));
h = pcolor( f, gain_bins, gain_hist(:,:,1) );
h.EdgeColor = 'none';
ax = gca();
ax.XScale       = 'log';
ax.XLim         = frange;
ax.YLabel.String = 'Gain';
ax.YScale       = 'log';
h = colorbar();
h.YLabel.String = 'Occurrence';
ax.Position     = inPos(idx,:);
title( ['Gain Bx ' upper(fgm_instr) '/SCM' fgm_title] );

% Gain Hist BY
idx = 4;
subplot(layout(1), layout(2), idx, 'OuterPosition', outPos(idx,:), 'Position', inPos(idx,:));
h = pcolor( f, gain_bins, gain_hist(:,:,2) );
h.EdgeColor = 'none';
ax = gca();
ax.XScale        = 'log';
ax.XLim          = frange;
ax.YLabel.String = 'Gain';
ax.YScale        = 'log';
h = colorbar();
h.YLabel.String  = 'Occurrence';
ax.Position      = inPos(idx,:);
title( ['Gain By ' upper(fgm_instr) '/SCM' fgm_title] );

% Gain Hist BZ
idx = 6;
subplot(layout(1), layout(2), idx, 'OuterPosition', outPos(idx,:), 'Position', inPos(idx,:));
h = pcolor( f, gain_bins, gain_hist(:,:,3) );
h.EdgeColor = 'none';
ax = gca();
ax.XScale        = 'log';
ax.XLim          = frange;
ax.YLabel.String = 'Gain';
ax.YScale        = 'log';
h = colorbar();
h.YLabel.String  = 'Occurrence';
ax.Position      = inPos(idx,:);
title( ['Gain Bz ' upper(fgm_instr) '/SCM' fgm_title] );


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

% Phase BX
idx = 1;
subplot(layout(1), layout(2), idx, 'OuterPosition', outPos(idx,:), 'Position', inPos(idx,:));
h = pcolor( t_dn, f, phase(:,:,1) );
h.EdgeColor = 'none';
ax = gca();
ax.XLim          = [t_dn(1) t_dn(end)];
ax.YLabel.String = 'f Hz';
ax.YLim          = frange;
ax.YScale        = 'log';
h = colorbar();
h.YLabel.String  = '\Delta\theta (deg)';
ax.Position      = inPos(idx,:);
title( ['\DeltaPhase Bx ' upper(fgm_instr) fgm_title] );

% Phase BY
idx = 3;
subplot(layout(1), layout(2), idx, 'OuterPosition', outPos(idx,:), 'Position', inPos(idx,:));
h = pcolor( t_dn, f, phase(:,:,1) );
h.EdgeColor = 'none';
ax = gca();
ax.XLim          = [t_dn(1) t_dn(end)];
ax.YLabel.String = 'f Hz';
ax.YLim          = frange;
ax.YScale        = 'log';
h = colorbar();
h.YLabel.String  = '\Delta\theta (deg)';
ax.Position      = inPos(idx,:);
title( ['\DeltaPhase By ' upper(fgm_instr) '-SCM' fgm_title] );

% Phase BZ
idx = 5;
subplot(layout(1), layout(2), idx, 'OuterPosition', outPos(idx,:), 'Position', inPos(idx,:));
h = pcolor( t_dn, f, phase(:,:,1) );
h.EdgeColor = 'none';
ax = gca();
ax.XLabel.String = 'Time';
ax.XLim          = [t_dn(1) t_dn(end)];
ax.YLabel.String = 'f Hz';
ax.YLim          = frange;
ax.YScale        = 'log';
h = colorbar();
h.YLabel.String  = '\Delta\theta (deg)';
ax.Position      = inPos(idx,:);
title( ['\DeltaPhase Bz ' upper(fgm_instr) '-SCM' fgm_title] );

hold on
plot( f_iwf, nf_iwf(3,:), '--r', 'LineWidth', 1.5 );
hold off


% Phase Hist BX
idx = 2;
subplot(layout(1), layout(2), idx, 'OuterPosition', outPos(idx,:), 'Position', inPos(idx,:));
h = pcolor( f, phase_bins, phase_hist(:,:,1) );
h.EdgeColor = 'none';
ax = gca();
ax.XScale        = 'log';
ax.XLim          = frange;
ax.YLabel.String = '\Delta\theta';
h = colorbar();
h.YLabel.String  = 'Occurrence';
ax.Position      = inPos(idx,:);
title( ['\DeltaPhase Bx ' upper(fgm_instr) '-SCM' fgm_title] )

% Phase Hist BY
idx = 4;
subplot(layout(1), layout(2), idx, 'OuterPosition', outPos(idx,:), 'Position', inPos(idx,:));
h = pcolor( f, phase_bins, phase_hist(:,:,2) );
h.EdgeColor = 'none';
ax = gca();
ax.XScale        = 'log';
ax.XLim          = frange;
ax.YLabel.String = '\Delta\theta';
h = colorbar();
h.YLabel.String  = 'Occurrence';
ax.Position      = inPos(idx,:);
title( ['\DeltaPhase By ' upper(fgm_instr) '-SCM' fgm_title] )

% Phase Hist BZ
idx = 6;
subplot(layout(1), layout(2), idx, 'OuterPosition', outPos(idx,:), 'Position', inPos(idx,:));
h = pcolor( f, phase_bins, phase_hist(:,:,3) );
h.EdgeColor = 'none';
ax = gca();
ax.XScale        = 'log';
ax.XLabel.String = 'f (Hz)';
ax.XLim          = frange;
ax.YLabel.String = '\Delta\theta';
h = colorbar();
h.YLabel.String  = 'Occurrence';
ax.Position      = inPos(idx,:);
title( ['\DeltaPhase Bz ' upper(fgm_instr) '-SCM' fgm_title] )


%------------------------------------%
% Save Figure                        %
%------------------------------------%
if ~isempty(png_dir)
	% Distributions
	set( fig1, 'PaperPositionMode', 'auto');
	filename = fullfile( png_dir, [sc '_fsm_' mode '_' level '_' optdesc '-gain_' tstart '.png'] );
	print( fig1, filename, '-dpng', '-r300' );
	disp( ['Saving file to: "' filename '".'] );
	
	% Cuts
	set( fig2, 'PaperPositionMode', 'auto');
	filename = fullfile( png_dir, [sc '_fsm_' mode '_' level '_' optdesc '-phase_' tstart '.png'] );
	print( fig2, filename, '-dpng', '-r300' );
	disp( ['Saving file to: "' filename '".'] );
end
