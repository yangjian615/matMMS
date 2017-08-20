%
% Name
%   mms_fsm_test_fit
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

% Data file to plot
% file = '/nfs/fsm/temp/mms4_fsm_srvy_l2plus_cal-dfg_20151101153733_v0.0.1.cdf';
%file    = '/nfs/fsm/temp/mms4_fsm_srvy_l2plus_cal-scm_20151103160519_v0.1.0.cdf';
file    = '/nfs/fsm/temp/mms4_fsm_brst_l2plus_cal-scm-day_20151101000000_v0.0.0.cdf';
theComp = 'z';
theFlag = 33;
theFreq = 12.0;
fc      = 0.5;            % Corner frequency of high-pass filter
fitmethod = 'process';    % curvefit | fminsearch | adaptive | process

%------------------------------------%
% Variable Names                     %
%------------------------------------%

% Parse file name to get variable name prefix and suffix
[sc, instr, mode, level, ~, ~, optdesc] = mms_dissect_filename(file);
prefix    = [sc '_' instr '_'];
suffix    = ['_' mode '_' level];
mag_instr = optdesc(5:end);

% Create Variable Names
amp_vname         = [prefix 'amp'   '_omb'   suffix];
phase_vname       = [prefix 'phase' '_omb'   suffix];
psd_vname         = [prefix 'psd'   '_omb'   suffix];
amp_hist_vname    = [prefix 'amp'   '_hist'  suffix];
phase_hist_vname  = [prefix 'phase' '_hist'  suffix];
psd_hist_vname    = [prefix 'psd'   '_hist'  suffix];
amp_floor_vname   = [prefix 'amp'   '_floor' suffix];
phase_floor_vname = [prefix 'phase' '_floor' suffix];
psd_floor_vname   = [prefix 'psd'   '_floor' suffix];

%------------------------------------%
% Read Data                          %
%------------------------------------%

% Spectra
% [amp, t, f] = MrCDF_Read(file, amp_vname,   'RowMajor', true);
% phase       = MrCDF_Read(file, phase_vname, 'RowMajor', true);
% [psd, t, f] = MrCDF_Read(file, psd_vname,   'RowMajor', true);

% Histogram
% [amp_hist,   ~, amp_bins, flag, comp] = MrCDF_Read(file, amp_hist_vname,   'RowMajor', true);
% [phase_hist, ~, phase_bins]           = MrCDF_Read(file, phase_hist_vname, 'RowMajor', true);
[psd_hist, f, psd_bins, flag, comp]     = MrCDF_Read(file, psd_hist_vname,   'RowMajor', true);

% Noise floor
% amp_floor   = MrCDF_Read(file, amp_floor_vname);
% phase_floor = MrCDF_Read(file, phase_floor_vname);
psd_floor   = MrCDF_Read(file, psd_floor_vname);

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
assert( ~isempty(iFlag), 'Invalid flag provided.' );

%------------------------------------%
% Create 2D Coordinate Grids         %
%------------------------------------%
%nTime     = length(t);
nFreq     = length(f);
nFlags    = length(flag);
nPSD_bins = length(psd_bins);

%psd      = squeeze(psd(:,iComp,:));
psd_hist = squeeze(psd_hist(:,iFlag,iComp,:));

% Convert to a 2D coordinate grid
%t_psd_grid = repmat( transpose( MrCDF_epoch2ssm(t) ), nFreq, 1 );
%f_psd_grid = repmat( transpose(f), 1, nTime );

f_hist_grid   = repmat( f, nPSD_bins, 1 );
bin_hist_grid = repmat( transpose(psd_bins), 1, nFreq );

%------------------------------------%
% Scale Data                         %
%------------------------------------%

% Scale the PSD
%sclPSD = MrRescale( log10( psd ), 1, 64, ...
%                    'MinValue', -6, ...
%                    'MaxValue', -1, ...
%                    'Class',    'uint8' );

% Scale the histogram
sclPSD_hist = MrRescale( psd_hist, 1, 64, ...
                         'Class', 'uint8' );

%------------------------------------%
% Pure MATLAB Fit                    %
%------------------------------------%
if strcmp(fitmethod, 'fminsearch')
	% Gaussian distribution
	%   - f = p1 exp( -(x - p2)^2 / (2 p3^2) ) + 
	%         p3 exp( -(x - p4)^2 / (2 p5^2) )
	%   - p1 = amplitude
	%   - p2 = mean
	%   - p3 = standard deviation
	% The goal is to minimize the sum of the least squares error.
	iFreq   = MrValue_Locate(f, theFreq);
	x       = psd_bins;
	y       = double(psd_hist(:,iFreq))';
	oneDist = y;
	fgauss  = @(p) sum( ( y - ( p(1) * exp( -(x - p(2)).^2 / (2*p(3)^2) ) + ...
	                            p(4) * exp( -(x - p(5)).^2 / (2*p(6)^2) ) ) ).^2 );
	
	% Allocate memory
	oneGuess = [50, -4.75, 0.5, 250, -1.5, 0.75];
	p        = fminsearch( fgauss, oneGuess );
	oneFit   = p(1) * exp( -(x - p(2)).^2 / (2*p(3)^2) ) + ...
	           p(4) * exp( -(x - p(5)).^2 / (2*p(6)^2) );

	% Loop through all points
	fitLo    = zeros(1,nFreq);
	fitHi    = zeros(1,nFreq);
	theGuess = [50, -4.75, 0.5, 250, -1.5, 0.75];
	iFreq    = find(f > 0.5, 1, 'first');
	for ii = iFreq : nFreq
		% Create a new function handle
		%   - Must do so to update the Y-data
		x      = psd_bins;
		y      = double(psd_hist(:,ii))';
		fgauss = @(p) sum( ( y - ( p(1) * exp( -(x - p(2)).^2 / (2*p(3)^2) ) + ...
		                           p(4) * exp( -(x - p(5)).^2 / (2*p(6)^2) ) ) ).^2 );
		
		% Get the fit parameters
		temp      = fminsearch( fgauss, theGuess );
		fitLo(ii) = temp(2);
		fitHi(ii) = temp(5);
		
		% Next iteration
		theGuess   = temp;
	end

%------------------------------------%
% Curve Fitting Toolbox              %
%------------------------------------%
elseif strcmp(fitmethod, 'curvefit')
	% Prepare data
	%   - Prevent warnings from fit().
	psd_bins = double(psd_bins');
	psd_hist = double(psd_hist);
	
	% Allocate memory
	iFreq   = MrValue_Locate(f, theFreq);
	oneDist = psd_hist(:,iFreq);
	oneFit  = fit( psd_bins, oneDist, 'gauss2', ...
	               'StartPoint', [50, -4.75, 0.5, 250, -1.5, 0.75] );
	
	% Loop through all points
	allFit = zeros(1,nFreq);
	fitOptions = fitoptions('gauss2', ...
	                        'StartPoint', [20, -6.5, 0.5, 400, -3.5, 1.0]);
	
	for ii = 1 : nFreq
		fprintf('Fit %i of %i.\n', ii, nFreq);
		temp       = fit( psd_bins, psd_hist(:,ii), 'gauss2', fitOptions );
		allFit(ii) = temp.b1;
		fitOptions.StartPoint = [temp.a1 temp.b1 temp.c1 temp.a2 temp.b2 temp.c2];
	end

%------------------------------------%
% Adaptive Method                    %
%------------------------------------%
elseif strcmp(fitmethod, 'adaptive')
	% Prepare data
	%   - Prevent warnings from fit().
	psd_bins = double(psd_bins');
	psd_hist = double(psd_hist);
	
	% Example fit.
	iFreq    = MrValue_Locate(f, theFreq);
	oneDist  = psd_hist(:,iFreq);
	oneFit   = fit( psd_bins, oneDist, 'gauss2' );

	% Initial guess
	iFreq     = MrValue_Locate(f, fc);
	initDist  = psd_hist(:,iFreq);
	initFit   = fit( psd_bins, initDist, 'gauss2' );
	initGuess = [ initFit.a1, initFit.b1, initFit.c1, initFit.a2, initFit.b2, initFit.c2 ];

	% Proceed with fminsearch
	fitLo    = zeros(1,nFreq);
	fitHi    = zeros(1,nFreq);
	ibin     = zeros(1,nFreq);
	theGuess = initGuess;
	for ii = iFreq : nFreq
		% Create a new function handle
		%   - Must do so to update the Y-data
		x      = psd_bins;
		y      = double(psd_hist(:,ii));
		fgauss = @(p) sum( ( y - ( p(1) * exp( -(x - p(2)).^2 / (2*p(3)^2) ) + ...
		                           p(4) * exp( -(x - p(5)).^2 / (2*p(6)^2) ) ) ).^2 );
		
		% Get the fit parameters
		temp      = fminsearch( fgauss, theGuess );
		fitLo(ii) = temp(2);
		fitHi(ii) = temp(5);

		% Nearest of fit
		[~, ibin(ii)] = min( ( psd_bins-temp(5) ).^2 );
		
		% Next iteration
		theGuess   = temp;
	end

%------------------------------------%
% Processing Method                  %
%------------------------------------%
elseif strcmp(fitmethod, 'process')
	% Prepare data
	dims     = size(psd_hist);
	psd_bins = double(psd_bins);
	psd_hist = reshape( double(psd_hist), 1, 1, dims(1), dims(2) );

	% Function used for data processing
	allFit = mms_fsm_bkgd_fit_bigauss( psd_hist, f, psd_bins, fc );
else
	error( ['Unknown fit method: "' fitmethod '".'] )
end

%------------------------------------%
% Plot Data                          %
%------------------------------------%
% Create a figure
fig = figure();

% Plot the PSD
subplot(3,1,1);
%h = pcolor(t_psd_grid, f_psd_grid, log10( psd ) );
%h.EdgeColor = 'none';
%colorbar

% Plot the histogram
subplot(3,1,2);
h = pcolor(f_hist_grid, bin_hist_grid, sclPSD_hist);
h.EdgeColor = 'none';
colorbar
hold on
h = plot(f, fitLo, '--k');
h = plot(f, fitHi, '--k');
line([theFreq, theFreq], [psd_bins(1) psd_bins(end)], 'LineStyle', '--', 'Color', [0,0,0])
hold off

% Plot a vertical slice of the histogram
subplot(3,1,3);
if strcmp(fitmethod, 'curvefit')
	h = plot(oneFit, psd_bins, oneDist);
elseif strcmp(fitmethod, 'adaptive')
	h = plot( oneFit, psd_bins, oneDist);
elseif strcmp(fitmethod, 'process')
	h = plot( psd_bins, allFit );
elseif strcmp(fitmethod, 'fminsearch')
	h = plot( psd_bins, oneDist, psd_bins, oneFit );
	legend({'data', 'fit'});
else

end