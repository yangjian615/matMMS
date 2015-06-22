function [] = edi_drift_step_plot ( ...
	obsID, ...
	B_tt2000, ...
	gd_virtual_bpp, ...
	gd_fv_bpp, ...
	DMPA2BPP, ...
	GrubbsBeamIntercepts, GrubbsBeamInterceptMean, GrubbsBeamInterceptMean_stdDev, ...
	P0);

	myLibAppConstants % custom colors; set default axis colors
	cPathSep = pathsep;
	cFileSep = filesep;
	global dotVersion

	hBPP_figure = figure ('Position', [ 400   150   850   800 ]);
	set (hBPP_figure, 'WindowStyle', 'normal')
	set (hBPP_figure, 'DockControls', 'off')
	set (gcf, 'name', 'B_E_field_GDU_locs_beams', 'NumberTitle', 'off', 'visible', 'on');
	set (0, 'CurrentFigure', hBPP_figure) % hBPP_plotElements
	clf

	EDI1gunLoc = [ -1.45598,  1.11837, 0.0 ]; % EDI2 gun atan2(-1.11837, 1.45598)*180/pi ~> -37.52865°
	EDI1detLoc = [ -1.35885,  1.03395, 0.0 ]; % EDI2 det atan2(-1.03395, 1.35885)*180/pi ~> -37.26753°
	EDI2gunLoc = [  1.45598, -1.11837, 0.0 ]; % EDI1detLoc:EDI1gunLoc angle = atan2(1.11837-1.03395, -1.45598+1.35885)*180/pi ~> -40.995°
	EDI2detLoc = [  1.35885, -1.03395, 0.0 ]; % norm(EDI1gunLoc-EDI1detLoc,2) = 0.128689
	mmsEDI_VirtualRadius = norm (EDI1gunLoc-EDI2detLoc, 2);

	theta = 0.0: 1.0: 360.0;
	virtual_instr_plane = mmsEDI_VirtualRadius * ...
		[ cosd(theta); sind(theta); zeros(1,length(theta),'double') ];

	BPP_plane = virtual_instr_plane; % re-use vip, because the BPP_plane looks normal in BPP
	hBPP_plotElements (1) = plot3 ( ...
		BPP_plane (1,:), ...
		BPP_plane (2,:), ...
		BPP_plane (3,:), ...
		'LineStyle', '-', 'LineWidth', 1.0, 'Color', myDarkBlue);

	AxisMax = 4;
	axis ([ -AxisMax AxisMax  -AxisMax AxisMax  -AxisMax AxisMax ]); % expanded axes for viewing larger drift steps
	axis square
	axis vis3d
	axis on
	grid on
	set (gca, 'XColor', myLightGrey4, 'YColor', myLightGrey4, 'ZColor', myLightGrey4)
	hold on
	xlabel ('X');
	ylabel ('Y');
	zlabel ('Z');
	view ([ 0 90 ])
	set (gcf, 'Units', 'normal')
	set (gca, 'Position', [0.01 0.01 0.99 0.99])
% keyboard
	GDU_planeInBPP = DMPA2BPP * virtual_instr_plane;
	hBPP_plotElements (2) = plot3 ( ...
		GDU_planeInBPP (1,:), ...
		GDU_planeInBPP (2,:), ...
		GDU_planeInBPP (3,:), ...
		'LineStyle', '-', 'LineWidth', 1.0, 'Color', myDarkRed);
	disp 'step 2 ~> DMPA in red'

	% 	edi_gun1_virtual_bpp = DMPA2BPP * edi_gun1_virtual_dmpa;
	% 	edi_gun2_virtual_bpp = DMPA2BPP * edi_gun2_virtual_dmpa;
% 		all_guns_virtual_bpp = [ edi_gun1_virtual_bpp, edi_gun2_virtual_bpp ];
	% 	edi_gd12_fv_bpp  = DMPA2BPP * edi_gd12_fv_dmpa;
	% 	edi_gd21_fv_bpp  = DMPA2BPP * edi_gd21_fv_dmpa;
% 		allgdxx_fv_bpp = [ edi_gd12_fv_bpp, edi_gd21_fv_bpp ];

	% GDU_BPP_coords = DMPA2BPP * GDU_OCS_coords;
	hBPP_plotElements (3) = plot3 ( ...
		gd_virtual_bpp (1,:), ...
		gd_virtual_bpp (2,:), ...
		gd_virtual_bpp (3,:), ...
		'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', myDarkRed, 'MarkerEdgeColor', myDarkRed, 'MarkerSize', 5.0);
	disp 'step 3 ~> GDUs as seen in BPP'
% keyboard
	for i=1: size (gd_virtual_bpp, 2)
		% The beams are // already [supposed to be] // parallel to BPP, so all we need to do is move them to the GDUs.
		% This should place them all in the same plane.
		GDU_Loc   = gd_virtual_bpp (:,i);
		BeamStart = GDU_Loc - 6.0 * gd_fv_bpp (:,i);
		BeamEnd   = GDU_Loc + 6.0 * gd_fv_bpp (:,i);
		BeamStartBPP = BeamStart; % OCS2BPP * BeamStart;
		BeamEndBPP   = BeamEnd; % OCS2BPP * BeamEnd;

		hBPP_bplotElements (4) = line ( ...
			[ BeamStartBPP(1) BeamEndBPP(1) ], ...
			[ BeamStartBPP(2) BeamEndBPP(2) ], ...
			[ BeamStartBPP(3) BeamEndBPP(3) ], ...
			'LineStyle', ':', 'LineWidth', 1.0, 'Color', myDarkBlue); % Beam
	end
	disp 'step 4 ~> firing vectors as seen in BPP'

	hBPP_plotElements (5) = plot3 (GrubbsBeamInterceptMean (1), GrubbsBeamInterceptMean (2), 0.0, ...
		'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', myOrange, 'MarkerEdgeColor', myOrange, 'MarkerSize', 5.0);

% 	disp ( sprintf ('Beam convergence: %+8.3f %+8.3f %+8.3f %+8.3f %+8.3f %+8.3f', ...
% 		GrubbsBeamInterceptMean, GrubbsBeamInterceptStdDev, GrubbsBeamInterceptMean_stdDev) )

	ConfidenceBounds = norminv (1.0 - (1.0 - P0) / 2.0);
	ConfidenceIntervals = [ ...
		GrubbsBeamInterceptMean-(ConfidenceBounds*GrubbsBeamInterceptMean_stdDev),... % x, y lower limits
		GrubbsBeamInterceptMean+(ConfidenceBounds*GrubbsBeamInterceptMean_stdDev) ];  % x, y upper limits
	disp (['Grubbs 84% confidence intervals   = ', sprintf('( %g, %g )', ConfidenceIntervals) ])
	if ( ((ConfidenceIntervals(1,2) - ConfidenceIntervals(1,1)) > 0.0) & ...
	     ((ConfidenceIntervals(2,2) - ConfidenceIntervals(2,1)) > 0.0) )
		r = rectangle ('Position', [ ...
			ConfidenceIntervals(1,1), ...
			ConfidenceIntervals(2,1), ...
			ConfidenceIntervals(1,2)-ConfidenceIntervals(1,1), ...
			ConfidenceIntervals(2,2)-ConfidenceIntervals(2,1) ], ...
			'LineStyle', '--', 'LineWidth', 1);
		set (r, 'edgecolor', myOrange)
	end

	for i = 1: size (GrubbsBeamIntercepts, 2)
		hBPP_plotElements (6) = plot3 ( ...
			GrubbsBeamIntercepts (1,i), ...
			GrubbsBeamIntercepts (2,i), ...
			0.0, ...
			'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', myGold, 'MarkerEdgeColor', myGold, 'MarkerSize', 2.0);
% 			zeros (1, length (GrubbsBeamIntercepts), 'double'), ...
	end

	nBeamIntercepts = length (GrubbsBeamIntercepts);
	bTimeStr = datestr (spdftt2000todatenum (B_tt2000), 'yyyy-mm-dd HH:MM:ss');
	title ( { ...
		[ 'EDI drift step using beam convergence, ', dotVersion, ', BPP'];
		[ 'MMS', obsID, ': B (UTC): ', bTimeStr, ', P_{0} = ', num2str(P0, 3), ' Points = ', num2str(nBeamIntercepts) ];
		[ '~> click plot to advance...'];
		[] }, 'Fontname', 'Times');

	% Example: 'mms4_edi_slow_l1a_efield_20150506_SDP__EDI_2D_driftstep_E_field__B_avg_3Dd'
	SavePlotFilename = [ ...
		'.' cFileSep 'mms', obsID, '_edi_', ...
		bTimeStr(1:4) bTimeStr(6:7) bTimeStr(9:10) '_drift_E__SDP__Bavg_', ...
		bTimeStr(12:13) bTimeStr(15:16) bTimeStr(18:19), ...
		'-', dotVersion, 'a.png' ];
	hgexport (gcf, SavePlotFilename, EDI_presentation_beam_plot_style);

% 	saveas (hBPP_figure, [ SavePlotFilename, '-v0101a.png' ], 'png');
% keyboard

% 	dummy = waitforbuttonpress;
	close (hBPP_figure);
end

% 00 - default
% 01 - Expand axes to fill figure - clips badly - NO
% 02 - Custom color = w for white frame
% 03a- Auto
% 03b- 600 dpi
% 04a- auto, painters
% 04b- a600 dpi, painters, not much change from 03a,b
% 05a- screen ?= auto
% 06 - 300 looks as good as 600
% 07 - lines >= 1 point - good
% 08 - lines, scale, min 1 point, - good
% assume your style sheet's name is <foo.txt>
%{
% create a fig
     plot(1:10);
     fnam='your_fig.png'; % your file name
% the engine
% ...get style sheet info
     snam='foo'; % note: NO extension...
     s=hgexport('readstyle',snam);
% ...apply style sheet info
     hgexport(gcf,fnam,s);


% create an example fig that we want to format with style file 'foo'
 plot(rand(14,10));

 % get style sheet info
 snam='foo'; % The name of your style file (NO extension)
 s = hgexport ('readstyle',snam);

 %apply style sheet info
 fnam='myfig.jpeg'; % your file name
 s.Format = 'jpeg'; %I needed this to make it work but maybe you wont.
 hgexport (gcf, fnam, s);

if you just want to apply the style to the figure itself,
the matlab-command hgexport(gcf,'temp_dummy','mystyle','applystyle', true);

s = hgexport ('readstyle','EDI beam plots')
s =
            Version: '1'
             Format: 'eps'
            Preview: 'none'
              Width: 'auto'
             Height: 'auto'
              Units: 'inches'
              Color: 'rgb'
         Background: 'w'
      FixedFontSize: '10'
     ScaledFontSize: 'auto'
           FontMode: 'scaled'
        FontSizeMin: '8'
     FixedLineWidth: '1'
    ScaledLineWidth: 'auto'
           LineMode: 'scaled'
       LineWidthMin: '1'
           FontName: 'auto'
         FontWeight: 'auto'
          FontAngle: 'auto'
       FontEncoding: 'latin1'
            PSLevel: '2'
           Renderer: 'painters'
         Resolution: '300'
       LineStyleMap: 'none'
         ApplyStyle: '0'
             Bounds: 'loose'
           LockAxes: 'on'
      LockAxesTicks: 'off'
             ShowUI: 'on'
       SeparateText: 'off'

s = hgexport ('readstyle','EDI presentation beam plots')
EDI_presentation_beam_plot_style =
            Version: '1'
             Format: 'eps'
            Preview: 'none'
              Width: 'auto'
             Height: 'auto'
              Units: 'inches'
              Color: 'rgb'
         Background: 'w'
      FixedFontSize: '10'
     ScaledFontSize: 'auto'
           FontMode: 'scaled'
        FontSizeMin: '8'
     FixedLineWidth: '1'
    ScaledLineWidth: 'auto'
           LineMode: 'fixed'
       LineWidthMin: '1'
           FontName: 'auto'
         FontWeight: 'auto'
          FontAngle: 'auto'
       FontEncoding: 'latin1'
            PSLevel: '2'
           Renderer: 'painters'
         Resolution: '300'
       LineStyleMap: 'none'
         ApplyStyle: '0'
             Bounds: 'loose'
           LockAxes: 'on'
      LockAxesTicks: 'on'
             ShowUI: 'on'
       SeparateText: 'off'
%}