%
% Name
%   mms_fsm_bkgd_plot_psdhist
%
% Purpose
%   Test different methods for fitting SCM data. Two methods are currently chosen:
%     1) fit guassian via the curve fitting toolbox
%     2) minimize least squared error to gaussian
%   The first method is much slower than the second but provides more statistical
%   information about the results. Initial conditions for the fit have been hard-
%   coded below. No quantitative comparison has yet been done.
%
%   Other methods can be found at this website (2016-09-09)
%     http://www.walkingrandomly.com/?p=5196
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2016-09-13      Written by Matthew Argall
%
%***************************************************************************
function fnames = mms_fsm_bkgd_plot_psdfloor(sc, mode, optdesc, tstart, tend, plot_dir)

	% Declare global attributes
	global cal_path_root data_path_root dropbox_root hk_root log_path_root unh_data_root

	% Defaults
	if nargin() < 6
		plot_dir = '';
	end
	
	% Initialize
	mms_fsm_init();
	
	% Constants
	instr     = 'fsm';
	level     = 'l2plus';
	mag_instr = regexp(optdesc, '-', 'split');
	mag_instr = mag_instr{2};
	
	% Find files
	files = mms_find_file( sc, instr, mode, level,      ...
	                       'Dropbox',   dropbox_root,   ...
	                       'OptDesc',   optdesc,        ...
	                       'SDCroot',   data_path_root, ...
	                       'TimeOrder', '%Y%M%d%H%m%S', ...
	                       'TStart' ,   tstart,         ...
	                       'TEnd',      tend );
	
	if ischar(files)
		files = { files };
	end
	
%------------------------------------%
% Variable Names                     %
%------------------------------------%

	% Parse file name to get variable name prefix and suffix
	prefix    = [sc '_' instr '_'];
	suffix    = ['_' mode '_' level];
	mag_instr = optdesc(5:end);

	floor_vname = [prefix 'psd' '_floor' suffix];

%------------------------------------%
% Read Data                          %
%------------------------------------%
	nFiles     = length(files);
	noiseFloor = cell(1, nFiles);
	hFlag      = cell(1, nFiles);

	% Loop through each file.
	for ii = 1 : nFiles
		[noiseFloor{ii}, f, hFlag{ii}, comp] = MrCDF_Read(files{ii}, floor_vname, 'RowMajor', true);
	end

	% Concatenate data along flag dimension
	hFlag      = cat( 2, hFlag{:} );
	noiseFloor = cat( 1, noiseFloor{:} );
	
	% LPP Noise Floor
	lpp_noise = mms_lpp_read_floor();

%------------------------------------%
% Setup Figure                       %
%------------------------------------%
	% Unique flags
	uniqFlags = unique(hFlag);
	nFlags    = length(uniqFlags);
	
	% Number of columns
	if all( ismember([0,1], bitget(uniqFlags, 6)) )
		nCols = 2;
	else
		nCols = 1;
	end

	% Number of rows
	%   - Number of uniq flags when the Deck32/64 flag is off
	%   - Plot one component per row for SCM
	%   - Plot one flag per row for FGM
	if strcmp(mag_instr, 'scm')
		nRows = 3;
	else
		fnames  = cell(1,3);
		rowFlag = unique( bitset(uniqFlags, 6, 0) );
		nRows   = length( rowFlag );
	end
	
	% Create figure
	xyz       = {'X', 'Y', 'Z'};
	xrange    = [ f(1) f(end) ];
	yrange    = [ min(noiseFloor(:)) max(noiseFloor(:)) ];
	yrange(1) = max( yrange(1), -10 );

%------------------------------------%
% Components                         %
%------------------------------------%
	% Loop over each component
	for iComp = 1 : length(comp)
	
		% Create a new figure for each component
		%   - One figure per component for FGM
		%   - One figure total for SCM
		if iComp == 1 || ~strcmp(mag_instr, 'scm')
			fig = figure( 'OuterPosition', [ 5, 5, 1000, 650] );
			
			% Create plot positions
			[inPos, outPos] = MrLayout( [nRows, nCols],      ...
			                            'Figure',   fig,     ...
			                            'IXMargin', [0,0],   ...
			                            'IYMargin', [0,0],   ...
			                            'OXMargin', [12,15], ...
			                            'OYMargin', [4,2],   ...
			                            'XGap',     1.5,     ...
			                            'YGap',     2 );
		end

	%------------------------------------%
	% Flags                              %
	%------------------------------------%
		for iFlag = 1 : nFlags
			% Elements with matching flags
			theFlag = uniqFlags(iFlag);
			idx     = find( hFlag == theFlag );
		
		%------------------------------------%
		% Current Plot Location              %
		%------------------------------------%
			% Column
			%   1 = Deck32
			%   2 = Deck64
			if nCols == 2 && bitget(theFlag, 6)
				col = 2;
			else
				col = 1;
			end
			
			% Row
			%   - Use component for SCM
			%   - Use flag for FGM
			if strcmp(mag_instr, 'scm')
				row = iComp;
			else
				row = find( rowFlag == bitset(theFlag, 6, 0) );
			end
			
			% Titles
			if bitget( theFlag, 1)
				theTitle = 'fast';
			else
				theTitle = 'slow';
				
				% Slow-lo-pre
				if all( bitget( theFlag, [2,4] ) )
					theTitle = [theTitle ' lo pre'];
					
				% Slow-hi
				elseif bitget( theFlag, 3 )
					theTitle = [theTitle ' hi'];
					
				% Slow-lo-post
				elseif all( bitget( theFlag, [2,5] ) )
					theTitle = [theTitle ' lo post'];
				else
					theTitle = [theTitle ' ???'];
%					error( ['Unknown bit combination (' num2str(uniqFlags(jj)) ').'] )
				end
			end
			if bitget(theFlag, 6)
				theTitle = ['Deck64 ' theTitle];
			else
				theTitle = ['Deck32 ' theTitle];
			end
		
		%------------------------------------%
		% Plot Individual Lines              %
		%------------------------------------%
			
			% Plot index
			iPlot = (row-1)*nCols + col;
			subplot( nRows, nCols, iPlot,        ...
			         'OuterPosition', outPos(iPlot,:), ...
			         'Position',      inPos(iPlot,:) );

			% Plot noise floor
			plot(f, squeeze(noiseFloor(idx,iComp,:)) );
			ax = gca();

		%------------------------------------%
		% Plot Mean & StdDev                 %
		%------------------------------------%
			% Mean and standard deviation
			meanFloor = squeeze( mean( noiseFloor(idx,iComp,:), 1 ) );
			stdFloor  = squeeze( std( noiseFloor(idx,iComp,:), 1, 1 ) );
			
			% Overplot w/error bars
			hold on
			h = plot(f, meanFloor, '--k');
			h.LineWidth = 2.0;
			hold off

		%------------------------------------%
		% Plot LPP Noise Floor               %
		%------------------------------------%
			hold on
			switch sc
				case 'mms1'
					h = plot( lpp_noise.mms1_f_nemi, log10(lpp_noise.mms1_b_nemi(iComp,:).^2), '--c' );
				case 'mms2'
					h = plot( lpp_noise.mms2_f_nemi, log10(lpp_noise.mms2_b_nemi(iComp,:).^2), '--c' );
				case 'mms3'
					h = plot( lpp_noise.mms3_f_nemi, log10(lpp_noise.mms3_b_nemi(iComp,:).^2), '--c' );
				case 'mms4'
					% Ground
					h = plot( lpp_noise.mms4_f_nemi, log10(lpp_noise.mms4_b_nemi(iComp,:).^2), '--c' );
					
					% In-situ
					if iComp == 3
						if strcmp(mode, 'srvy')
							hz = plot( lpp_noise.f_srvy_insitu, log10(lpp_noise.bz_srvy_insitu.^2), '--m' );
						else
							hz = plot( lpp_noise.f_brst_insitu, log10(lpp_noise.bz_brst_insitu.^2), '--m' );
						end
						hz.LineWidth = 2.0;
					end
				otherwise
					error( 'Invalid value for SC.' );
			end
			h.LineWidth = 2.0;
			hold off
		
		%------------------------------------%
		% Make Pretty                        %
		%------------------------------------%
			title( [upper(mag_instr) ' ' xyz{iComp} ' ' theTitle] );
			ax.YLim = [-7,-2]; %yrange;
			ax.XLim = xrange;
			
			% BOTTOM
			if row == nRows
				ax.XLabel.String = 'f (Hz)';
			else
				ax.XTickLabel = [];
			end
			
			% RIGHT
			if col == 1
				ax.YLabel.String = {'PSD', 'log_{10}nT^{2}/Hz'};
			else
				ax.YTickLabel = [];
			end
		end % iFlag

	%------------------------------------%
	% Save Component (FGM)               %
	%------------------------------------%
		if ~strcmp(mag_instr, 'scm') && ~isempty(plot_dir)
			% Start and end times without delimiters
			fstart = regexp(tstart, '[-T:]', 'split');
			fstart = strcat(fstart{:});
			fend   = regexp(tend, '[-T:]', 'split');
			fend   = strcat(fend{:});

			% File name
			fnames{iComp} = sprintf('%s_%s_%s_%s_%s-b%s-psd-floor-zoom_%s_%s.png', ...
			                        sc, instr, mode, level, optdesc, lower(xyz{iComp}), fstart, fend);

			% Save
			print(fig, fullfile(plot_dir, fnames{iComp}), '-dpng');
		end
	end % iComp
	

%------------------------------------%
% Save Figure (FGM)                  %
%------------------------------------%
	if strcmp(mag_instr, 'scm') && ~isempty(plot_dir)
		% Start and end times without delimiters
		fstart = regexp(tstart, '[-T:]', 'split');
		fstart = strcat(fstart{:});
		fend   = regexp(tend, '[-T:]', 'split');
		fend   = strcat(fend{:});
		
		% File name
		fnames = sprintf('%s_%s_%s_%s_%s-b-psd-floor_%s_%s.png', ...
		                 sc, instr, mode, level, optdesc, fstart, fend);
		
		% Save
		print(fig, fullfile(plot_dir, fnames), '-dpng');
	end
end
