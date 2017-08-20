%
% Name
%   mms_fsm_nf_read
%
% Purpose
%   Read FSM noise floor data.
%
% Calling Sequence
%   DATA = mms_fsm_nf_read(FILES)
%     Read FSM noise floor data from files named FILES and return at
%     a data structure DATA.
%
% Parameters
%   FILES:          in, required, type=string/cell
%
% Returns
%   DATA            out, required, type=struct
%                   Structure with the following fields:
%                     't'     -  Time (cdf_time_tt2000)
%                     'dt'    -  Time interval over which the noise floor is valid.
%                     'nf'    -  Noise floor
%                     'std'   -  Standard deviation of the noise floor
%                     'flag'  -  Operational flag
%                     'comp'  -  Component index (0=x, 1=y, 2=z)
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2016-10-23      Written by Matthew Argall
%
function [] = mms_fsm_nf_plot(sc, mode, tstart, tend)

	% Declare global attributes
	global cal_path_root data_path_root dropbox_root hk_root log_path_root unh_data_root
	
	% Initialize
	mms_fsm_init();

	sc        = 'mms4';
	instr     = 'fsm';
	mode      = 'srvy';
	level     = 'l2plus';
	tstart    = '2015-11-01T00:00:00';
	tend      = '2015-11-01T24:00:00';
	fgm_instr = 'dfg';
	plot_dir  = '/nfs/fsm/figures/';
	
	% Initial time
	fstart  = MrTimeParser( tstart, '%Y-%M-%dT%H:%m:%S', '%Y%M%d%H%m%S' );
	[~, t0] = mms_parse_time(fstart);
	
%------------------------------------%
% Find & Read Files                  %
%------------------------------------%
	% FGM
	f_nf_fgm = mms_find_file( sc, instr, mode, level,                     ...
	                          'Dropbox',       dropbox_root,              ...
	                          'OptDesc',       ['nf-' fgm_instr '-week'], ...
	                          'SDCroot',       data_path_root,            ...
	                          'TimeOrder',     '%Y%M%d%H%m%S',            ...
	                          'TStart',        tstart,                    ...
	                          'TEnd',          tend,                      ...
	                          'RelaxedTStart', true );
	
	% SCM
	f_nf_scm = mms_find_file( sc, instr, mode, level,         ...
	                          'Dropbox',       dropbox_root,   ...
	                          'OptDesc',       'nf-scm-week',  ...
	                          'SDCroot',       data_path_root, ...
	                          'TimeOrder',     '%Y%M%d%H%m%S', ...
	                          'TStart',        tstart,         ...
	                          'TEnd',          tend,           ...
	                          'RelaxedTStart', true );

	% Read
	nf_fgm = mms_fsm_nf_read( f_nf_fgm );
	nf_scm = mms_fsm_nf_read( f_nf_scm );
	
%------------------------------------%
% Create the Weight Function         %
%------------------------------------%
	weight = mms_fsm_nf_weight( nf_fgm, nf_scm, t0 );
	comp   = [0 1 2];


%------------------------------------%
% Setup Figure                       %
%------------------------------------%
	% Unique flags
	nFlags = length( weight.flag );
	
	% Number of columns
	if all( ismember([0,1], bitget(weight.flag, 6)) )
		nCols = 2;
	else
		nCols = 1;
	end

	% Number of rows
	%   - Number of uniq flags when the Deck32/64 flag is off
	%   - Plot one flag per row
	fnames  = cell(1,3);
	rowFlag = unique( bitset(weight.flag, 6, 0) );
	nRows   = length( rowFlag );
	
	% Create figure
	xyz       = {'X', 'Y', 'Z'};
	xrange    = [ weight.f(1) weight.f(end) ];
	yrange    = [ -0.1, 1.1 ];

%------------------------------------%
% Components                         %
%------------------------------------%
	% Loop over each component
	for iComp = 1 : length(comp)
	
		% Create a new figure for each component
		%   - One figure per component
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

	%------------------------------------%
	% Flags                              %
	%------------------------------------%
		for iFlag = 1 : nFlags
			% Elements with matching flags
			theFlag = weight.flag(iFlag);
		
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
			
			% Row (== flag)
			row = find( rowFlag == bitset(theFlag, 6, 0) );
			
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

			% Weight function
			plot( weight.f, squeeze( weight.w(:,iComp,iFlag)) );
			ax = gca();
		
		%------------------------------------%
		% Make Pretty                        %
		%------------------------------------%
			title( [ 'W(f) ' xyz{iComp} ' ' theTitle] );
			ax.YLim = yrange;
			ax.XLim = xrange;
			
			% BOTTOM
			if row == nRows
				ax.XLabel.String = 'f (Hz)';
			else
				ax.XTickLabel = [];
			end
			
			% RIGHT
			if col == 1
				ax.YLabel.String = 'Weight';
			else
				ax.YTickLabel = [];
			end
		end % iFlag

	%------------------------------------%
	% Save Component (FGM)               %
	%------------------------------------%
		if ~isempty(plot_dir)
			% Start and end times without delimiters
			fstart = regexp(tstart, '[-T:]', 'split');
			fstart = strcat(fstart{:});
			fend   = regexp(tend, '[-T:]', 'split');
			fend   = strcat(fend{:});

			% File name
			fnames{iComp} = sprintf('%s_%s_%s_%s_%s-b%s-weight_%s_%s.png', ...
			                        sc, instr, mode, level, 'wf', lower(xyz{iComp}), fstart, fend);

			% Save
			print(fig, fullfile(plot_dir, fnames{iComp}), '-dpng');
			disp( ['File saved to: ' fnames{iComp} ] );
		end
	end % iComp
end