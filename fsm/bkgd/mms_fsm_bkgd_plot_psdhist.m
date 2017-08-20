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
%   2016-09-09      Written by Matthew Argall
%
%***************************************************************************
function fnames = mms_fsm_bkgd_plot_psdhist(sc, mode, optdesc, tstart, tend, plot_dir)

	% Declare global attributes
	global cal_path_root data_path_root dropbox_root hk_root log_path_root unh_data_root

	% Defaults
	if nargin() < 6
		plot_dir = '';
	end
	
	% Initialize
	mms_fsm_init();
	
	% Constants
	instr    = 'fsm';
	level    = 'l2plus';
	
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
% LPP Noise Floor                    %
%------------------------------------%
	lpp_noise = mms_lpp_read_floor();
	
%------------------------------------%
% Variable Names                     %
%------------------------------------%

	% Parse file name to get variable name prefix and suffix
	prefix    = [sc '_' instr '_'];
	suffix    = ['_' mode '_' level];
	mag_instr = regexp(optdesc, '-', 'split');
	mag_instr = mag_instr{2};

	hist_vname  = [prefix 'psd' '_hist'   suffix];
	floor_vname = [prefix 'psd' '_floor'  suffix];

%------------------------------------%
% Loop Over Each File                %
%------------------------------------%
		
	% Figures
	%   - Orbit plots
	%   - Combined noise floor plots
	f1 = figure('Position', [ 50, 50,1200,600], 'Units', 'pixels');

	% Allocate memory
	nFiles = length(files);
	fnames = cell(1, nFiles);
	xyz    = {'X', 'Y', 'Z'};

	% Loop over files
	for ii = 1 : nFiles
	%------------------------------------%
	% Read Data                          %
	%------------------------------------%
		% Read data
		[hist, f, bins, flag, comp] = MrCDF_Read(files{ii}, hist_vname, 'RowMajor', true);
		noiseFloor                  = MrCDF_Read(files{ii}, floor_vname);

	%------------------------------------%
	% Setup Figure                       %
	%------------------------------------%
		% Number of columns
		nComp = length(comp);
		nFlag = length(flag);
	
		% Create figure
		xyz    = {'X', 'Y', 'Z'};
		range  = [min( hist(:) ) max( hist(:) )];
		xrange = [min(f) max(f)];
		yrange = [min(bins) max(bins)];
		
		% Create plot positions
		[inPosHist, outPosHist] = MrLayout( [nFlag, nComp],     ...
		                                    'Figure',   f1,      ...
		                                    'IXMargin', [0,0],   ...
		                                    'IYMargin', [0,0],   ...
		                                    'OXMargin', [12,15], ...
		                                    'OYMargin', [4,2],   ...
		                                    'XGap',     2,       ...
		                                    'YGap',     2 );

	%------------------------------------%
	% Loop Over Component                %
	%------------------------------------%
		for iComp = 1 : nComp

		%------------------------------------%
		% Loop Over Flags                    %
		%------------------------------------%
			for iFlag = 1 : nFlag
		
			%------------------------------------%
			% Create Surface Plot                %
			%------------------------------------%
				% Scale the histogram data
				sclHist = MrRescale( squeeze(hist(:,iFlag,iComp,:)), 1, 64, ...
				                     'MinValue', range(1),            ...
				                     'MaxValue', range(2) );
				
				% Current position
				iPlot = iComp + (iFlag-1)*nComp;
				subplot( nFlag, nComp, iPlot, ...
				         'OuterPosition', outPosHist(iPlot,:), ...
				         'Position',      inPosHist(iPlot,:) );
				
				% Create the surface plot
				%   - Orient so the xy-plane is flat on the screen
				surf( f, bins, sclHist, 'EdgeColor', 'none' );
				ax = gca();
				ax.View = [0, 90];
				ax.XLim = xrange;
				ax.YLim = yrange;
%				ax.CLim = range;

				% Overplot the noise floor
				hold on
				plot3(f, noiseFloor(:,iFlag,iComp), repmat( range(2), 1, length(f) ), '--k', 'LineWidth', 2 );
				hold off
		
			%------------------------------------%
			% Overplot LLP Noise Floor           %
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
				% LEFT Edge
				if iComp == 1
					ax.YLabel.String = {'PSD', 'log_{10}(nT^{2}/Hz)'};
				else
					ax.YTickLabel = [];
				end
				
				% RIGHT Edge
				if iComp == nComp
					cb = colorbar(ax);
					cb.Label.String = 'Counts';
					ax.Position = inPosHist(iPlot,:);
				end
				
				% BOTTOM row
				if iFlag == nFlag
					ax.XLabel.String = 'f (Hz)';
				else
					ax.XTickLabel = [];
				end
				
				% TITLE
				%   1 = slow/fast/brst
				%   2 = lorange
				%   3 = hirange
				%   4 = pre-perigee
				%   5 = post-perigee
				%   6 = deck 32/64
				theTitle = [ upper(mag_instr) ' ' xyz{comp(iComp)+1} ];
				tf_flag  = bitget(flag(iFlag), [1, 2, 3, 4, 5, 6]);
				
				% Mode
				if tf_flag(1)
					if strcmp(mode, 'brst')
						theTitle = [theTitle ' brst'];
					else
						theTitle = [theTitle ' fast'];
					end
				else
					theTitle = [theTitle ' slow'];
				end
				if tf_flag(2)
					theTitle = [theTitle ' lo'];
				end
				if tf_flag(3)
					theTitle = [theTitle ' hi'];
				end
				if tf_flag(4)
					theTitle = [theTitle ' pre'];
				end
				if tf_flag(5)
					theTitle = [theTitle ' post'];
				end
				if tf_flag(6)
					theTitle = [theTitle ' deck32'];
				else
					theTitle = [theTitle ' deck64'];
				end
				title( ax, theTitle );
			end % flag
		end % comp

	%------------------------------------%
	% Plot Each Chart                    %
	%------------------------------------%
		% Save each figure
		if ~isempty(plot_dir)
			[~, ~, ~, ~, fstart] = mms_dissect_filename(files{ii});
			fnames{ii} = sprintf('%s_%s_%s_%s_%s-psd-hist_%s.png', sc, instr, mode, level, optdesc, fstart);
			print(f1, fullfile(plot_dir, fnames{ii}), '-dpng');
			
			% Clear figure1
			clf(f1);
		end

	end % files
	
	% Close the figure
	if nargout() == 0
		clear fnames
	else
		close(f1);
	end
end
