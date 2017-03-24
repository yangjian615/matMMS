%
% Name
%   mms_fsm_bkgd_fit_savgol
%
% Purpose
%   Find the noise floor of a distribution by fitting it to a gaussian function.
%
%   Gaussian distribution
%     - f = p1 exp( -(x - p2)^2 / (2 p3^2) )
%     - p1 = amplitude
%     - p2 = mean
%     - p3 = standard deviation
%   The goal is to minimize the sum of the least squares error.
%
% Calling Sequence
%   H_FLOOR = mms_fsm_bkgd_compute(F, H, BINS)
%     Fit each histogram distribution in H at each frequency F to a gaussian
%     distribution. Evaluate at each histogram bin location BINS. The peak value
%     constitutes the noise floor of the distribution, H_FLOOR.
%
%   H_FLOOR = mms_fsm_bkgd_compute(__, P)
%     Supply the initial guess for the gaussian as a vector of three elements,
%     P = [amplitude, mean, deviation].
%
% Parameters
%   F               in, required, type = float
%   H               in, required, type = integer
%   BINS            in, required, type = double
%   P               in, optional, type = double, default = [30.0, -5.0, 1.0]
%
% Returns
%   H_FLOOR         out, required, type=single
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2016-09-14      Written by Matthew Argall
%
function [h_floor, h_sig, fit0] = mms_fsm_bkgd_fit_savgol(h, f, bins, fc, f0)

	% Allocate memory
	%   - h = [comp, flag, bin, freq]
	dims    = size(h);
	nComp   = dims(1);
	nFlag   = dims(2);
	nBins   = dims(3);
	nFreq   = dims(4);
	h_floor = zeros(nComp, nFlag, nFreq);
	
	
	ifc = MrValue_Locate(f, fc);
	if0 = MrValue_Locate(f, f0);
	
	
	% Loop through all points
	for iComp = 1 : nComp
	for iFlag = 1 : nFlag
		result = sgolayfilt( squeeze(h(iComp,iFlag,:,:)), 6, 31, [], 2 );
	end
	end

%------------------------------------%
% Plot Results                       %
%------------------------------------%
	binRange = [min(bins), max(bins)]

	f1 = figure();
	
	% Create plot positions
	[inPosHist, outPosHist] = MrLayout( [3, 1],              ...
	                                    'Figure',   f1,      ...
	                                    'IXMargin', [0,0],   ...
	                                    'IYMargin', [0,0],   ...
	                                    'OXMargin', [12,15], ...
	                                    'OYMargin', [4,2],   ...
	                                    'XGap',     2,       ...
	                                    'YGap',     2 );

	% Scale the histogram data
	sclIn  = MrRescale( squeeze( h(3,1,:,:) ), 1, 64 );
	sclOut = MrRescale( result, 1, 64 );
	
	subplot( 3, 1, 1, ...
	         'OuterPosition', outPosHist(1,:), ...
	         'Position',      inPosHist(1,:) );
	
	% Input histogram
	%   - Orient so the xy-plane is flat on the screen
	surf( f, bins, sclIn, 'EdgeColor', 'none' );
	ax = gca();
	ax.View = [0, 90];
	ax.Title.String = 'In-Situ';
	ax.YLim = binRange;
	ax.YLabel.String = {'PSD', 'log_{10}(nT^2/Hz)'};
	
	
	
	subplot( 3, 1, 2, ...
	         'OuterPosition', outPosHist(2,:), ...
	         'Position',      inPosHist(2,:) );
	
	% Output histogram
	%   - Orient so the xy-plane is flat on the screen
	surf( f, bins, sclIn, 'EdgeColor', 'none' );
	ax = gca();
	ax.View = [0, 90];
	ax.Title.String = 'Savitky-Golay Smoothing';
	ax.XLabel.String = {'f (Hz)'};
	ax.YLim = binRange;
	ax.YLabel.String = {'PSD', 'log_{10}(nT^2/Hz)'};
	line( [f0,f0], binRange, [64,64], 'LineStyle', '--', 'LineWidth', 2, 'Color', 'r' );
	
	
	
	subplot( 3, 1, 3, ...
	         'OuterPosition', outPosHist(3,:), ...
	         'Position',      inPosHist(3,:) );

	% Fit
	plot(bins, sclIn(:,if0), bins, sclOut(:,if0) );
	ax = gca()
	ax.Title.String = '1D Cut';
	ax.XLabel.String = {'PSD', 'log_{10}(nT^2/Hz)'};
	ax.YLabel.String = 'Occurrence';
	legend( {'Data', 'Savgol'} )
	
	
keyboard
	
end