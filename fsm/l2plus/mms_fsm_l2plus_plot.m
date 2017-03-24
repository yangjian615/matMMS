%
% Name
%   mms_fsm_l2plus_plot
%
% Purpose
%   Create a comparison plot between the FGM and FSM data products. Comparison
%   is performed in the time domian by plotting FGM on top of FSM.
%
% Calling Sequence
%   mms_fsm_l2plus_plot( SC, MODE, TSTART )
%     Create a time series plot of Bx, By, Bz from the FSM and FGM data products.
%     Data is from spacecraft SC, in telemetry mode MODE, with file start time of
%     TSTART, formatted as yyyy-mm-dd for survey mode and yyyy-mm-ddTHH:MM:SS for
%     burst mode.
%
%   FNAMES = mms_fsm_l2plus_plot( ..., PLOT_DIR )
%     Save the plot to the directory specified by PLOT_DIR and return the file names
%     to FNAMES.
%
% Inputs
%   SC          in, required, type=char
%   MODE        in, required, type=char
%   TSTART      in, required, type=char
%   PLOT_DIR    in, optional, type=char, default=''
%
% Outputs
%   FNAMES      out, optional, type=char
%
% MATLAB release(s) 9.0.0.341360 (R2016a)
% Required Products None
%
% History:
%   2016-10-02      Written by Matthew Argall
%
%***************************************************************************
function fname = mms_fsm_l2plus_plot(sc, mode, tstart, tend, plot_dir)

	% Declare global attributes
	global cal_path_root data_path_root dropbox_root hk_root log_path_root unh_data_root

	% Defaults
	if nargin() < 6
		plot_dir = '';
	end
	
	% Initialize
	mms_fsm_init();
	
	% FGM
%	sc        = 'mms1';
%	mode      = 'brst';
%	tstart    = '2015-11-01T08:28:04';
%	tend      = '2015-11-01T08:28:04';
	instr     = 'fsm';
	level     = 'l2plus';
	fgm_instr = 'dfg';
	coord_sys = 'dmpa';
	
%------------------------------------%
% Find Files                         %
%------------------------------------%
	switch mode
		case 'srvy'
			optdesc_scm = 'scsrvy';
		case 'brst'
			optdesc_scm = 'scb';
		otherwise
			error( ['Invalid value for MODE: "' mode '".'] );
	end
	
	
	% FSM
	fsm_files = mms_find_file( sc, instr, mode, level,   ...
	                           'Dropbox',   dropbox_root,   ...
	                           'SDCroot',   data_path_root, ...
	                           'TimeOrder', '%Y%M%d%H%m%S', ...
	                           'TStart' ,   tstart,         ...
	                           'TEnd',      tend );
	
	% FGM
	fgm_files = mms_find_file( sc, fgm_instr, mode, 'l2pre',  ...
	                           'Dropbox',   dropbox_root,     ...
	                           'SDCroot',   data_path_root,   ...
	                           'TStart' ,   tstart,           ...
	                           'TEnd',      tend );
	
	% SCM
%	scm_files = mms_find_file( sc, 'scm', mode, 'l2',         ...
%	                           'Dropbox',   dropbox_root,     ...
%	                           'OptDesc',   optdesc_scm,      ...
%	                           'SDCroot',   data_path_root,   ...
%	                           'TStart' ,   tstart,           ...
%	                           'TEnd',      tend );

%------------------------------------%
% Variable Names                     %
%------------------------------------%
	fsm_prefix = [sc '_fsm_'];
	fgm_prefix = [sc '_' fgm_instr '_'];
	scm_prefix = [sc '_scm_'];
	
	fsm_suffix = ['_' mode '_l2plus'];
	fgm_suffix = ['_' mode '_l2pre_' coord_sys];
	scm_suffix = ['_' mode '_l2'];

	% Parse file name to get variable name prefix and suffix
	b_fsm_vname = [ fsm_prefix 'b' '_' coord_sys fsm_suffix        ];
	b_fgm_vname = [ fgm_prefix                   fgm_suffix(2:end) ];
%	b_scm_vname = [ scm_prefix 'b' '_' coord_sys scm_suffix        ];

%------------------------------------%
% Read Data                          %
%------------------------------------%
	[b_fsm, t_fsm] = MrCDF_Read( fsm_files, b_fsm_vname );
	[b_fgm, t_fgm] = MrCDF_Read( fgm_files, b_fgm_vname );
%	[b_scm, t_scm] = MrCDF_Read( scm_files, b_scm_vname );

	% Convert time to seconds since midnight
	[~, ~, ~, ~, fstart] = mms_dissect_filename( fsm_files );
	[~, t_ref]           = mms_parse_time( fstart );
%	t_ref     = min( [ t_fsm(1) t_fgm(1) ] );
	t_fsm_ssm = MrCDF_epoch2sse( t_fsm, t_ref );
	t_fgm_ssm = MrCDF_epoch2sse( t_fgm, t_ref );
	
	% Trim data
	trange = [ max( t_fsm_ssm(1), t_fgm_ssm(1) ) min( t_fsm_ssm(end), t_fgm_ssm(end) ) ];
	ifsm   = MrValue_Locate( t_fsm_ssm, trange );
	ifgm   = MrValue_Locate( t_fgm_ssm, trange );
	t_fsm_ssm = t_fsm_ssm( ifsm(1):ifsm(2) );
	t_fgm_ssm = t_fgm_ssm( ifgm(1):ifgm(2) );
	b_fsm  = b_fsm( :, ifsm(1):ifsm(2) );
	b_fgm  = b_fgm( :, ifgm(1):ifgm(2) );

%------------------------------------%
% FSM - FGM                          %
%------------------------------------%
	% Interpolate FGM to FSM
	db      = zeros( size(b_fsm) );
	db(1,:) = spline( t_fgm_ssm, b_fgm(1,:), t_fsm_ssm );
	db(2,:) = spline( t_fgm_ssm, b_fgm(2,:), t_fsm_ssm );
	db(3,:) = spline( t_fgm_ssm, b_fgm(3,:), t_fsm_ssm );

	% Take the difference with FSM
	db = abs( db - b_fsm );

%------------------------------------%
% Setup Figure                       %
%------------------------------------%
	% Create the figures
	f1 = figure( 'OuterPosition', [ 25, 25, 1000, 650] );
	
	% Rows and columns
	nRows = 3;
	nCols = 2;
	
	% Create plot positions
	[inPos, outPos] = MrLayout( [nRows, nCols],      ...
	                            'Figure',   f1,      ...
	                            'IXMargin', [0,0],   ...
	                            'IYMargin', [0,0],   ...
	                            'OXMargin', [12,15], ...
	                            'OYMargin', [4,2],   ...
	                            'XGap',     8,       ...
	                            'YGap',     2 );

%------------------------------------%
% Plot Time Series                   %
%------------------------------------%
	yrange = [-10,50];
	dyrange = [0,1];
%	xrange  = hms_to_ssm({'15:15:00', '18:15:00'});
	xrange = hms_to_ssm({'15:56:36', '15:56:38'});

	
	% Bx Position
	subplot( nRows, nCols, 1,              ...
	         'OuterPosition', outPos(1,:), ...
	         'Position',      inPos(1,:) );

	% Bx Plot
	h = plot( t_fsm_ssm, b_fsm(1,:), t_fgm_ssm, b_fgm(1,:) );
	title( [upper(sc) ' ' fstart] );
	legend( {'FSM', 'FGM'} );
	ax               = gca();
	ax.YLabel.String = {'Bx' '(nT)'};
	ax.XTickLabel    = [];
	ax.YLim          = yrange;
	ax.XLim          = xrange;

	% By Position
	subplot( nRows, nCols, 3,              ...
	         'OuterPosition', outPos(3,:), ...
	         'Position',      inPos(3,:) );
	
	% By Plot
	h = plot( t_fsm_ssm, b_fsm(2,:), t_fgm_ssm, b_fgm(2,:) );
	legend( {'FSM', 'FGM'} );
	ax               = gca();
	ax.YLabel.String = {'By' '(nT)'};
	ax.XTickLabel    = [];
	ax.YLim          = yrange;
	ax.XLim          = xrange;
	
	% Bz Position
	subplot( nRows, nCols, 5,              ...
	         'OuterPosition', outPos(5,:), ...
	         'Position',      inPos(5,:) );

	% Bz Plot
	plot( t_fsm_ssm, b_fsm(3,:), t_fgm_ssm, b_fgm(3,:) );
	legend( {'FSM', 'FGM'} );
	ax               = gca();
	ax.YLabel.String = {'Bz' '(nT)'};
	ax.XLabel.String = 'Time (ssm)';
	ax.YLim          = yrange;
	ax.XLim          = xrange;

%------------------------------------%
% Plot Differences                   %
%------------------------------------%
	
	% dBx Position
	subplot( nRows, nCols, 2,              ...
	         'OuterPosition', outPos(2,:), ...
	         'Position',      inPos(2,:) );
	

	% dBx Plot
	h = plot(t_fsm_ssm, db(1,:));
	ax = gca();
	ax.Title.String  = {'|B_{FGM} - B_{FSM}|'};
	ax.YLabel.String = {'dB_{x} (nT)'};
	ax.YLim = dyrange;
	ax.XLim = xrange;
	
	% dBy Position
	subplot( nRows, nCols, 4,              ...
	         'OuterPosition', outPos(4,:), ...
	         'Position',      inPos(4,:) );
	
	% dBz Plot
	h = plot(t_fsm_ssm, db(2,:));
	ax = gca();
	ax.YLabel.String = {'dB_{y} (nT)'};
	ax.YLim = dyrange;
	ax.XLim = xrange;
	
	% dBx Position
	subplot( nRows, nCols, 6,              ...
	         'OuterPosition', outPos(6,:), ...
	         'Position',      inPos(6,:) );
	
	% dBz Plot
	h = plot(t_fsm_ssm, db(3,:));
	ax = gca();
	ax.XLabel.String = {'Time (ssm)'};
	ax.YLabel.String = {'dB_{x} (nT)'};
	ax.YLim = [2,7];
	ax.XLim = xrange;
	
%------------------------------------%
% Save                               %
%------------------------------------%
	if ~isempty(plot_dir)
		% Start and end times without delimiters
		fstart = regexp(tstart, '[-T:]', 'split');
		fstart = strcat(fstart{:});
		fend   = regexp(tend, '[-T:]', 'split');
		fend   = strcat(fend{:});
		
		% File name
		fname = sprintf('%s_%s_%s_%s_comp-fgm_%s_%s.png', ...
		                sc, 'fsm', mode, 'l2plus', fstart, fend);
		
		% Save
		print(f1, fullfile(plot_dir, fname), '-dpng');
	end
end
