%
% Name
%   mms_fsm_l2plus_plot_psd
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
function fname = mms_fsm_l2plus_plot_psd(sc, mode, tstart, tend, plot_dir)

	% Declare global attributes
	global cal_path_root data_path_root dropbox_root hk_root log_path_root unh_data_root

	% Defaults
	if nargin() < 4
		plot_dir = '';
	end
	
	% Initialize
	mms_fsm_init();
	
	% FGM
%	sc          = 'mms4';
%	mode        = 'brst';
%	tstart      = '2015-11-01T08:28:04';
%	tend        = '2015-11-01T08:28:04';
	fgm_instr   = 'dfg';
	coord_sys   = 'gse';
	f_HeavySide = 4.0;
	plot_dir    = '/nfs/fsm/figures/';
	
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
	fsm_files = mms_find_file( sc, 'fsm', mode, 'l2plus',   ...
	                           'Dropbox',   dropbox_root,   ...
	                           'SDCroot',   data_path_root, ...
	                           'TStart' ,   tstart,         ...
	                           'TEnd',      tend );
	
	% FGM
	fgm_files = mms_find_file( sc, fgm_instr, mode, 'l2pre',  ...
	                           'Dropbox',   dropbox_root,     ...
	                           'SDCroot',   data_path_root,   ...
	                           'TStart' ,   tstart,           ...
	                           'TEnd',      tend );
	
	% SCM
	scm_files = mms_find_file( sc, 'scm', mode, 'l2',         ...
	                           'Dropbox',   dropbox_root,     ...
	                           'OptDesc',   optdesc_scm,      ...
	                           'SDCroot',   data_path_root,   ...
	                           'TStart' ,   tstart,           ...
	                           'TEnd',      tend );

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
	b_fsm_vname = [ fsm_prefix 'b'   '_' coord_sys                 fsm_suffix        ];
	b_fgm_vname = [ fgm_prefix                                     fgm_suffix(2:end) ];
	b_scm_vname = [ scm_prefix 'acb' '_' coord_sys '_' optdesc_scm scm_suffix        ];

%------------------------------------%
% Read Data                          %
%------------------------------------%
	[b_fsm, t_fsm] = MrCDF_Read( fsm_files, b_fsm_vname );
	[b_fgm, t_fgm] = MrCDF_Read( fgm_files, b_fgm_vname );
	[b_scm, t_scm] = MrCDF_Read( scm_files, b_scm_vname );

	% Convert time to seconds since midnight
	t_ref     = min( [ t_fsm(1) t_fgm(1) ] );
	t_fsm_ssm = MrCDF_epoch2ssm( t_fsm, t_ref );
	t_fgm_ssm = MrCDF_epoch2ssm( t_fgm, t_ref );
	t_scm_ssm = MrCDF_epoch2ssm( t_scm, t_ref );

%------------------------------------%
% Select Sub-Interval                %
%------------------------------------%
	%
	% The first and last segment of FSM data are not correct. Trim
	% the first and last 10% of the data, then find the corresponding
	% segments in FGM and SCM.
	%
	N    = size( b_fsm, 2 );
	Ncut = floor( 0.1 * N );
	
	% FSM
	t_fsm_ssm = t_fsm_ssm( Ncut+1:end-Ncut );
	b_fsm     = b_fsm(:, Ncut+1:end-Ncut );
	
	% FGM
	iRange    = MrValue_Locate( t_fgm_ssm, [ t_fsm_ssm(1) t_fsm_ssm(end) ] );
	t_fgm_ssm = t_fgm_ssm( iRange(1)+1:iRange(2) );
	b_fgm     = b_fgm(1:3, iRange(1)+1:iRange(2) );
	
	% SCM
	iRange    = MrValue_Locate( t_scm_ssm, [ t_fsm_ssm(1) t_fsm_ssm(end) ] );
	t_scm_ssm = t_scm_ssm( iRange(1)+1:iRange(2) );
	b_scm     = b_scm( :, iRange(1)+1:iRange(2) );
	
	% Free data
	clear t_ref N Ncut iRange t_fsm t_fgm t_scm

%------------------------------------%
% High-Pass Filter FSM               %
%------------------------------------%
	% Cutoff frequency
	fc = 0.1;

	% Create filter
	%   - butter() + filtfilt() was producing NaNs. This was the solution
	%   - https://www.mathworks.com/matlabcentral/answers/62553-filtfilt-function-returning-nan-at-certain-frequencies
	sr = round( 1.0 / median( diff(t_fsm_ssm) ) );
	fN = sr / 2.0;
	d  = fdesign.highpass('N,F3dB', 9, fc/fN);
	h  = design(d, 'butter');

	% Apply the filter
	b_fsm(1,:) = filtfilt( h.sosMatrix, h.ScaleValues, double( b_fsm(1,:) ) );
	b_fsm(2,:) = filtfilt( h.sosMatrix, h.ScaleValues, double( b_fsm(2,:) ) );
	b_fsm(3,:) = filtfilt( h.sosMatrix, h.ScaleValues, double( b_fsm(3,:) ) );

%------------------------------------%
% High-Pass Filter FGM               %
%------------------------------------%
	% Create filter
	%   - butter() + filtfilt() was producing NaNs. This was the solution
	%   - https://www.mathworks.com/matlabcentral/answers/62553-filtfilt-function-returning-nan-at-certain-frequencies
	sr = round( 1.0 / median( diff(t_fgm_ssm) ) );
	fN = sr / 2.0;
	d  = fdesign.highpass('N,F3dB', 9, fc/fN);
	h  = design(d, 'butter');

	% Apply the filter
	b_fgm(1,:) = filtfilt( h.sosMatrix, h.ScaleValues, double( b_fgm(1,:) ) );
	b_fgm(2,:) = filtfilt( h.sosMatrix, h.ScaleValues, double( b_fgm(2,:) ) );
	b_fgm(3,:) = filtfilt( h.sosMatrix, h.ScaleValues, double( b_fgm(3,:) ) );

%------------------------------------%
% Compute PSD                        %
%------------------------------------%
	% Sampling rates
	sr_fsm = round( 1.0 / median( diff( t_fsm_ssm ) ) );
	sr_fgm = round( 1.0 / median( diff( t_fgm_ssm ) ) );
	sr_scm = round( 1.0 / median( diff( t_scm_ssm ) ) );
	
	% Sampling interval
	dt_fsm = 1.0 / sr_fsm;
	dt_fgm = 1.0 / sr_fgm;
	dt_scm = 1.0 / sr_scm;
	
	% Number of points
	N_fsm = length( t_fsm_ssm );
	N_fgm = length( t_fgm_ssm );
	N_scm = length( t_scm_ssm );
	nfft_fsm = N_fsm - mod(N_fsm, 2);
	nfft_fgm = N_fgm - mod(N_fgm, 2);
	nfft_scm = N_scm - mod(N_scm, 2);
	
	% Frequency interval
	df_fsm = 1.0 / (nfft_fsm * dt_fsm);
	df_fgm = 1.0 / (nfft_fgm * dt_fgm);
	df_scm = 1.0 / (nfft_scm * dt_scm);
	
	% Frequencies
	f_fsm = [0:nfft_fsm/2] * df_fsm;
	f_fgm = [0:nfft_fgm/2] * df_fgm;
	f_scm = [0:nfft_scm/2] * df_scm;
	
	% FFT
	b_fft_fsm = fft( b_fsm, [], 2 );
	b_fft_fgm = fft( b_fgm, [], 2 );
	b_fft_scm = fft( b_scm, [], 2 );

	% Power
	psd_fsm = ( 2.0 * dt_fsm / nfft_fsm ) * abs( b_fft_fsm(:, 1:nfft_fsm/2+1).^2 );
	psd_fgm = ( 2.0 * dt_fgm / nfft_fgm ) * abs( b_fft_fgm(:, 1:nfft_fgm/2+1).^2 );
	psd_scm = ( 2.0 * dt_scm / nfft_scm ) * abs( b_fft_scm(:, 1:nfft_scm/2+1).^2 );
	
	% DC Component
	psd_fsm(:, [1,end]) = psd_fsm(:, [1,end]) / 2.0;
	psd_fgm(:, [1,end]) = psd_fgm(:, [1,end]) / 2.0;
	psd_scm(:, [1,end]) = psd_scm(:, [1,end]) / 2.0;

	% Clear data
	clear sr_fsm sr_fgm sr_scm dt_fsm dt_fgm dt_scm N_fsm N_fgm N_scm ...
	      nfft_fsm nfft_fgm nfft_scm df_fsm df_fgm df_scm b_fft_fsm b_fft_fgm b_fft_scm

%------------------------------------%
% Setup Figure                       %
%------------------------------------%
	% Create the figures
	f1 = figure( 'OuterPosition', [ 25, 25, 1000, 650] );
	
	% Rows and columns
	nRows = 3;
	nCols = 1;
	
	% Create plot positions
	[inPos, outPos] = MrLayout( [nRows, nCols],      ...
	                            'Figure',   f1,      ...
	                            'IXMargin', [0,0],   ...
	                            'IYMargin', [0,0],   ...
	                            'OXMargin', [12,15], ...
	                            'OYMargin', [4,2],   ...
	                            'XGap',     1,       ...
	                            'YGap',     2 );

%------------------------------------%
% Plot Time Series                   %
%------------------------------------%
	
	% Bx Position
	subplot( nRows, nCols, 1,              ...
	         'OuterPosition', outPos(1,:), ...
	         'Position',      inPos(1,:) );

	% Bx Plot
	h = loglog( f_fsm, psd_fsm(1,:), f_scm, psd_scm(1,:), f_fgm, psd_fgm(1,:) );
	title( [upper(sc) ' ' tstart] );
	legend( {'FSM', 'SCM', 'FGM'} );
	ax               = gca();
	ax.YLabel.String = {'Bx' '(nT^2/Hz)'};
	ax.XTickLabel    = [];
	
	% Line at HeavySide freuqency
	line( [f_HeavySide, f_HeavySide], ax.YLim, 'Color', 'r', 'LineStyle', '--')
	
	% By Position
	subplot( nRows, nCols, 2,              ...
	         'OuterPosition', outPos(2,:), ...
	         'Position',      inPos(2,:) );

	% By Plot
	h = loglog( f_fsm, psd_fsm(2,:), f_scm, psd_scm(2,:), f_fgm, psd_fgm(2,:) );
	legend( {'FSM', 'SCM', 'FGM'} );
	ax               = gca();
	ax.YLabel.String = {'By' '(nT^2/Hz)'};
	ax.XTickLabel    = [];
	
	% Line at HeavySide freuqency
	line( [f_HeavySide, f_HeavySide], ax.YLim, 'Color', 'r', 'LineStyle', '--')
	
	% Bz Position
	subplot( nRows, nCols, 3,              ...
	         'OuterPosition', outPos(3,:), ...
	         'Position',      inPos(3,:) );

	% Bz Plot
	loglog( f_fsm, psd_fsm(3,:), f_scm, psd_scm(3,:), f_fgm, psd_fgm(3,:) );
	legend( {'FSM', 'SCM', 'FGM'} );
	ax               = gca();
	ax.YLabel.String = {'Bz' '(nT^2/Hz)'};
	ax.XLabel.String = 'f (Hz)';
	
	% Line at HeavySide freuqency
	line( [f_HeavySide, f_HeavySide], ax.YLim, 'Color', 'r', 'LineStyle', '--')
	
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
		fname = sprintf('%s_%s_%s_%s_comp-psd_%s_%s.png', ...
		                sc, 'fsm', mode, 'l2plus', fstart, fend);
		
		% Save
		print(f1, fullfile(plot_dir, fname), '-dpng');
	end
end
