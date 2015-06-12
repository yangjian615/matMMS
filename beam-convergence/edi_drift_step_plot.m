function [] = edi_drift_step_plot ( ...
	b_tt2000, ...
	edi_gun1_virtual_bpp, edi_gun2_virtual_bpp, ...
	edi_gd12_fv_bpp, edi_gd21_fv_bpp, ...
	DMPA2BPP, ...
	GrubbsBeamInterceptMean, GrubbsBeamInterceptMean_stdDev);

	myColors        % MMS custom colors; set default axis colors
	cPathSep = pathsep;
	cFileSep = filesep;

	% if it doesn't exist, create it
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
	title ( { ...
		['1. The BPP coordinate system'];
		['~> click plot to advance...'];
		[] }, 'Fontname', 'Times');

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

	GDU_planeInBPP = DMPA2BPP * virtual_instr_plane;
	hBPP_plotElements (2) = plot3 ( ...
		GDU_planeInBPP (1,:), ...
		GDU_planeInBPP (2,:), ...
		GDU_planeInBPP (3,:), ...
		'LineStyle', '-', 'LineWidth', 1.0, 'Color', myDarkRed);
	disp 'step 2 ~> DMPA in red'

	% 	edi_gun1_virtual_bpp = DMPA2BPP * edi_gun1_virtual_dmpa;
	% 	edi_gun2_virtual_bpp = DMPA2BPP * edi_gun2_virtual_dmpa;
		all_guns_virtual_bpp = [ edi_gun1_virtual_bpp, edi_gun2_virtual_bpp ];
	% 	edi_gd12_fv_bpp  = DMPA2BPP * edi_gd12_fv_dmpa;
	% 	edi_gd21_fv_bpp  = DMPA2BPP * edi_gd21_fv_dmpa;
		allgdxx_fv_bpp = [ edi_gd12_fv_bpp, edi_gd21_fv_bpp ];

	% GDU_BPP_coords = DMPA2BPP * GDU_OCS_coords;
	hBPP_plotElements (3) = plot3 ( ...
		all_guns_virtual_bpp (1,:), ...
		all_guns_virtual_bpp (2,:), ...
		all_guns_virtual_bpp (3,:), ...
		'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', myDarkRed, 'MarkerEdgeColor', myDarkRed, 'MarkerSize', 5.0);
	disp 'step 3 ~> GDUs as seen in BPP'

	for i=1: length (all_guns_virtual_bpp)
		% The beams are // already [supposed to be] // parallel to BPP, so all we need to do is move them to the GDUs.
		% This should place them all in the same plane.
		GDU_Loc   = all_guns_virtual_bpp (:,i);
		BeamStart = GDU_Loc - 6.0 * allgdxx_fv_bpp (:,i);
		BeamEnd   = GDU_Loc + 6.0 * allgdxx_fv_bpp (:,i);
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

	ConfidenceBounds_z84 = 1.405; % norminv (0.92); % norminv (1.0 - (1.0 - 0.84) / 2.0);
	ConfidenceIntervals = [ ...
		GrubbsBeamInterceptMean-(ConfidenceBounds_z84*GrubbsBeamInterceptMean_stdDev),... % x, y lower limits
		GrubbsBeamInterceptMean+(ConfidenceBounds_z84*GrubbsBeamInterceptMean_stdDev) ];  % x, y upper limits
	disp (['Grubbs 84% confidence intervals   = ', sprintf('( %g, %g )', ConfidenceIntervals) ])
	r = rectangle ('Position', [ ...
		ConfidenceIntervals(1,1), ...
		ConfidenceIntervals(2,1), ...
		ConfidenceIntervals(1,2)-ConfidenceIntervals(1,1), ...
		ConfidenceIntervals(2,2)-ConfidenceIntervals(2,1) ], ...
		'LineStyle', '--', 'LineWidth', 3);
	set (r, 'edgecolor', myOrange)

	% mms4_edi_slow_l1a_efield_20150506_SDP__EDI_2D_driftstep_E_field__B_avg_3Dd
	bTimeStr = datestr (spdftt2000todatenum (b_tt2000), 'yyyy-mm-dd HH:MM:ss.fff');
	SavePlotFilename = [ ...
		'.' cFileSep 'mms4_edi_', ...
		bTimeStr(1:4) bTimeStr(6:7) bTimeStr(9:10) '_drift_E__SDP__Bavg_', ...
		bTimeStr(12:13) bTimeStr(15:16) bTimeStr(18:23), ...
		];
	saveas (hBPP_figure, [ SavePlotFilename, '.png' ], 'png');

	dummy = waitforbuttonpress;
	close (hBPP_figure);
end
