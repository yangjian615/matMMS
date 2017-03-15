%
% Name
%   mms_fsm_test_fgm_fit
%
% Purpose
%   Test different methods for fitting SCM data. Two methods are currently chosen:
%     1) fit bi-gaussian via the curve fitting toolbox
%     2) minimize least squared error to bi-gaussian
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

function [data, h, t, f, bins] = mms_fsm_test_histfit(file, varargin)

	% Data file to plot
	file        = '/nfs/fsm/temp/mms4_fsm_srvy_l2plus_cal-dfg_20151101153733_v1.0.0.cdf';
%	file      = '/nfs/fsm/temp/mms4_fsm_srvy_l2plus_cal-scm_20151103160519_v1.0.0.cdf';
	theComp     = 'z';
	theFlag     = 1;
	bin_range   = [];
	bin_size    = 0.1;
	datatype    = 'psd';
	f0          = 1.0;
	fc          = 0.5;
	f_range     = [];
	f_bin_size  = [];
	fit_method  = 'savgol';
	
	% Optional Arguments
	nOptArgs = length( varargin );
	for ii = 1 : 2 : nOptArgs
		switch varargin{ii}
			case 'BinRange'
				bin_range = varargin{ii+1};
			case 'BinSize'
				bin_size = varargin{ii+1};
			case 'DataType'
				datatype = varargin{ii+1};
			case 'FBinSize'
				f_bin_size = varargin{ii+1};
			case 'FDiagnostic'
				fc = varargin{ii+1};
			case 'FInit'
				f0 = varargin{ii+1};
			case 'FRange'
				f_range = varargin{ii+1};
			case 'FitMethod'
				fit_method = varargin{ii+1};
			otherwise
				error( ['Optional argument not recognized: "' varargin{ii} '".'] )
		end
	end

%------------------------------------%
% Variable Names                     %
%------------------------------------%

	% Parse file name to get variable name prefix and suffix
	[sc, instr, mode, level, ~, ~, optdesc] = mms_dissect_filename(file);
	prefix    = [sc '_' instr '_'];
	suffix    = ['_' mode '_' level];
	mag_instr = optdesc(5:end);

	% Create Variable Names
	flag_vname  = [prefix 'flag'         suffix];
	amp_vname   = [prefix 'amp'   '_omb' suffix];
	phase_vname = [prefix 'phase' '_omb' suffix];
	psd_vname   = [prefix 'psd'   '_omb' suffix];

%------------------------------------%
% Read Data                          %
%------------------------------------%
	% Flag
	flag = MrCDF_Read(file, flag_vname, 'RowMajor', true);

	% Spectra
	switch datatype
		case 'psd'
			[data, t, f, comp] = MrCDF_Read(file, psd_vname, 'RowMajor', true);
		case 'amp'
			[data, t, f, comp] = MrCDF_Read(file, amp_vname, 'RowMajor', true);
		case 'phase'
			[data, t, f, comp] = MrCDF_Read(file, phase_vname, 'RowMajor', true);
		otherwise
			error( 'DATATYPE must be {"amp" | "phase" | "psd"}.' )
	end

	% Transpose data to be
	%   [ component, time, frequency ]
	data = permute( data, [2,3,1] );

%------------------------------------%
% Pick Indices                       %
%------------------------------------%

	% Component
	switch theComp
		case 'x'
			iComp = 1;
		case 'y'
			iComp = 2;
		case 'z'
			iComp = 3;
		otherwise
			error( ['Invalid component: "' theComp '".'] );
	end

	% Flag
	iFlag = find(flag == theFlag);
	assert( ~isempty(iFlag), ['Invalid flag provided (' sprintf('%0i', theFlag) ').'] );
	
	% Prune data
	t    = t(iFlag);
	data = data(:,iFlag,:);

%------------------------------------%
% Datatype-Dependent Defaults        %
%------------------------------------%
	
	% Amplitude
	if strcmp( datatype, 'amp' )
		data      = alog10( data );
		units     = 'nT';
		if isempty(bin_range)
			bin_range = [-4.0, 4.0];
		end
	
	% Phase
	elseif strcmp( datatype, 'phase' )
		units     = 'Degrees';
		if isempty(bin_range)
			bin_range = [-180.0, 180.0];
		end
	
	% PSD
	elseif strcmp( datatype, 'psd' )
		data      = log10( data );
		units     = 'log_{10}(nT^2/Hz)';
		if isempty(bin_range)
			bin_range = [-10, 5];
		end
	
	% Other
	else
		error( ['Invalid value for DATATYPE: "' datatype '".'] );
	end

%------------------------------------%
% Histogram Data                     %
%------------------------------------%
	% Create bins
	nTimes  = length(t);
	nFreqs  = length(f);
	nBins   = 1 + floor( ( bin_range(2) - bin_range(1) ) ./ bin_size );
	bins    = linspace( bin_range(1), bin_range(2), nBins);
	
	
	if isempty(f_range)
		f_range = [ min(f) max(f) ];
	end
	if isempty(f_bin_size)
		nFHist = length(f);
		fHist  = f;
	else
		nFHist = 1 + floor( ( f_range(2) - f_range(1) ) ./ f_bin_size );
		fHist  = linspace( f_range(1), f_range(2), nFHist);
	end

	% Allocate memory
	h = zeros( 3, nBins, nFreqs );
	
	% Loop over all frequencies
	for jj = 1 : nFreqs
		h(1,:,jj) = histc( data(1,:,jj), bins );
		h(2,:,jj) = histc( data(2,:,jj), bins );
		h(3,:,jj) = histc( data(3,:,jj), bins );
	end
	
	if nFHist ~= nFreqs
		% Locate frequencies F within the new grid FREQS
		%   - Find unique frequencies by their index
		%   - Sum non-unique indices together
		ifreq  = MrValue_Locate(fHist, f, 'RoundUp', true);
		iufreq = unique( ifreq );
		nufreq = length( iufreq );
		temp   = zeros( 3, nBins, nFHist );

		for jj = 1 : nFHist
			ifnew             = iufreq(jj);
			ifsum             = find( ifreq == iufreq(jj) );
			temp( :,:,ifnew ) = sum( h(:,:,ifsum), 3 );
		end
		
		h = temp;
		clear temp
	end

%------------------------------------%
% Determine Noise Floor              %
%------------------------------------%
	% Prepare data
	%   [ component, flag, bin, freq ]
	dims = size(h);
	bins = double(bins);
	h    = reshape( double(h), dims(1), 1, dims(2), dims(3) );

	% Function used for data processing
	switch fit_method
		case 'gauss'
			[floorFit, fit0] = mms_fsm_bkgd_fit_gauss( h, fHist, bins, fc, f0 );
		case 'bigauss'
			[floorFit, sigFit, fit0] = mms_fsm_bkgd_fit_bigauss( h, fHist, bins, fc, f0 );
		case 'gauss12'
			[floorFit, sigFit, fit0] = mms_fsm_bkgd_fit_gauss12( h, fHist, bins, fc, f0 );
		case 'bigauss noise'
			[floorFit, sigFit, fit0] = mms_fsm_bkgd_fit_bigauss_noise( h, fHist, bins, fc, f0 );
		case 'trigauss'
			[floorFit, sigFit, noiseFit, fit0] = mms_fsm_bkgd_fit_trigauss( h, fHist, bins, fc, f0 );
		case 'savgol'
			result = mms_fsm_bkgd_fit_savgol( h, fHist, bins, fc, f0 );
		otherwise
			error( ['Invalid fitting method: "' fit_method '".'] )
	end
	nGauss = numcoeffs( fit0{1} ) / 3;

	% Initial fit index
	ifc = MrValue_Locate(fHist, fc);

%------------------------------------%
% Setup Figure                       %
%------------------------------------%
	strComp = {'X', 'Y', 'Z'};

	% Create a figure
	fig = figure( 'OuterPosition', [ 5, 5, 1100, 700] );
	
	% Create plot positions
	[inPos, outPos] = MrLayout( [3, 3],              ...
	                            'Figure',   fig,     ...
	                            'IXMargin', [0,0],   ...
	                            'IYMargin', [0,0],   ...
	                            'OXMargin', [12,15], ...
	                            'OYMargin', [4,2],   ...
	                            'XGap',     3,       ...
	                            'YGap',     4 );

%------------------------------------%
% Plot Results                       %
%------------------------------------%
	nCol = length(comp);
	for iCol = 1 : nCol

	%------------------------------------%
	% Plot Spectra                       %
	%------------------------------------%
		% Plot Data
		iRow  = 1;
		iPlot = (iRow-1)*nCol + iCol;
		subplot( 3, 3, iPlot, ...
		         'OuterPosition', outPos(iPlot,:), ...
		         'Position',      inPos(iPlot,:) );

		pSpec = pcolor( MrCDF_epoch2ssm(t), f, squeeze( data(iCol,:,:) )' );
		pSpec.EdgeColor = 'none';
		title( [ upper(sc) ' ' upper(mag_instr) ' ' strComp{iCol} ] );
		axSpec               = gca();
		axSpec.XLabel.String = 'Time (sec)';

	%------------------------------------%
	% Plot Histogram                     %
	%------------------------------------%
	
		% Plot the histogram
		iRow  = 2;
		iPlot = (iRow-1)*nCol + iCol;
		subplot( 3, 3, iPlot, ...
		         'OuterPosition', outPos(iPlot,:), ...
		         'Position',      inPos(iPlot,:) );

		pHist = pcolor( fHist, bins, squeeze( h(iCol,1,:,:) ) );
		pHist.EdgeColor      = 'none';
		axHist               = gca();
		axHist.XLabel.String = 'f (Hz)';
		colorbar
		
		% Overplot the resulting fit
		hold on
		if nGauss == 1
			plot(fHist, squeeze( floorFit(iCol,1,:) ), '--g');
		end
		if nGauss == 2
			plot(fHist, squeeze( floorFit(iCol,1,:) ), '--g');
			plot(fHist, squeeze( sigFit(iCol,1,:) ),   '--m');
		end
		if nGauss == 3
			plot(fHist, squeeze( floorFit(iCol,1,:) ), '--g');
			plot(fHist, squeeze( sigFit(iCol,1,:) ),   '--m');
			plot(fHist, squeeze( noiseFit(iCol,1,:) ), '--r');
		end
		line([f0, f0], [bins(1) bins(end)], 'LineStyle', '--', 'Color', 'k');
		line([fc, fc], [bins(1) bins(end)], 'LineStyle', '--', 'Color', 'c');
		hold off

	%------------------------------------%
	% Plot Initial Fit                   %
	%------------------------------------%
		
		% Plot the Initial fit
		iRow  = 3;
		iPlot = (iRow-1)*nCol + iCol;
		subplot( 3, 3, iPlot, ...
		         'OuterPosition', outPos(iPlot,:), ...
		         'Position',      inPos(iPlot,:) );
		theFit = fit0{iComp,1};
		plot( theFit, bins, squeeze( h(iCol,1,:,ifc) ) );

		hold on
		if nGauss == 3
			gauss1 = theFit.a1 * exp( -( ( bins - theFit.b1 ) / theFit.c1 ).^2.0 );
			plot( bins, gauss1, '--k' );
			legend( {'Data', 'Fit', 'Gauss1'} );
		end
		if nGauss == 2
			gauss1 = theFit.a1 * exp( -( ( bins - theFit.b1 ) / theFit.c1 ).^2.0 );
			gauss2 = theFit.a2 * exp( -( ( bins - theFit.b2 ) / theFit.c2 ).^2.0 );
			plot( bins, gauss1, '--g', bins, gauss2, '--r', bins, gauss1 + gauss2 , '--k' );
			legend( {'Data', 'Fit', 'Gauss1', 'Gauss2', 'Tot'} );
		end
		if nGauss == 3
			gauss1 = theFit.a1 * exp( -( ( bins - theFit.b1 ) / theFit.c1 ).^2.0 );
			gauss2 = theFit.a2 * exp( -( ( bins - theFit.b2 ) / theFit.c2 ).^2.0 );
			gauss3 = theFit.a3 * exp( -( ( bins - theFit.b3 ) / theFit.c3 ).^2.0 );
			plot( bins, gauss1, '--g', bins, gauss2, '--m', bins, gauss3, '--r', bins, gauss1 + gauss2 + gauss3, '--k' );
			legend( {'Data', 'Fit', 'Gauss1', 'Gauss2', 'Gauss3', 'Tot'} );
		end
		hold off
		
		axFit               = gca();
		axFit.XLabel.String = [ upper( datatype) ' (' units ')' ];

	%------------------------------------%
	% Make Pretty                        %
	%------------------------------------%
		if iCol == 1
			axSpec.YLabel.String = {'f', '(Hz)'};
			axHist.YLabel.String = { upper(datatype), ['(' units ')'] };
			axFit.YLabel.String  = '#';
		else
			axSpec.YTickLabel = [];
			axHist.YTickLabel = [];
			axFit.YTickLabel  = [];
			axFit.YLabel.String = '';
		end
	end

%------------------------------------%
% Clean Up                           %
%------------------------------------%
	if nargout() == 0
		clear data h t f bins
	end
end