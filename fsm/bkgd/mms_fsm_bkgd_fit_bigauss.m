%
% Name
%   mms_fsm_bkgd_fit_bigauss
%
% Purpose
%   Find the noise floor of a distribution by fitting it to a bi-gaussian function.
%
%   Gaussian distribution
%     - f = p1 exp( -(x - p2)^2 / (2 p3^2) ) + 
%           p4 exp( -(x - p5)^2 / (2 p6^2) )
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
%   H_FLOOR = mms_fsm_bkgd_compute(__, FC)
%     Give the corner frequency at which to begin fitting, FC. This is usually the
%     corner frequency of a high-pass filter.
%
%   H_FLOOR = mms_fsm_bkgd_compute(__, FC, F0)
%     Give the frequency of the initial fit, F0. The initial fit is performed using
%     the curve fitting toolbox. Subsequent fits are performed iteratively using
%     fminsearch.
%
%   [..., H_SIG] = mms_fsm_bkgd_compute(__)
%     Return the peak value of the second gaussian distribution. This is interpreted
%     as the mean signal strength.
%
%   [..., FIT0] = mms_fsm_bkgd_compute(__)
%     Return the solutions to the initial fits in FIT0.
%
% Parameters
%   F               in, required, type = float
%   H               in, required, type = integer
%   BINS            in, required, type = double
%   FC              in, optional, type = double, default = 0.5
%
% Returns
%   H_FLOOR         out, required, type=single
%   H_SIG           out, optional, type=single
%   FIT0            out, optional, type=cfit object
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2016-09-14      Written by Matthew Argall
%   2016-09-26      Create an initial guess with the curve fitting toolbox. - MRA
%   2016-10-04      Turn off output from fminsearch. Report warning at end. - MRA
%   2016-10-19      Added the FIT0 output parameter. - MRA
%
function [h_floor, h_sig, fit0] = mms_fsm_bkgd_fit_bigauss(h, f, bins, fc, f0)

	% Initial guess
	if nargin < 5 || isempty(f0)
		f0 = 1.0;
	end
	if nargin < 4 || isempty(fc)
		fc = 0.5;
	end

%------------------------------------%
% Allocate Memory to Results         %
%------------------------------------%
	
	% Allocate memory
	%   - h = [comp, flag, bin, freq]
	dims    = size(h);
	nComp   = dims(1);
	nFlag   = dims(2);
	nBins   = dims(3);
	nFreq   = dims(4);
	h_floor = zeros(nComp, nFlag, nFreq, 'single');
	h_sig   = zeros(nComp, nFlag, nFreq, 'single');
	fit0    = cell(nComp, nFlag);
	
	% Convert to double precision to prevent warnings
	if ~isa(h, 'double')
		h = double(h);
	end
	if ~isa(bins, 'double')
		bins = double(bins);
	end

%------------------------------------%
% Configure Loop                     %
%------------------------------------%
	
	% Suppress output of fminsearch
	options  = optimset('Display', 'off');
	nMaxIter = 0;

	% Loop through components and flags
	ifc = MrValue_Locate(f, fc);
	if0 = MrValue_Locate(f, f0);
	for iComp = 1 : nComp
	for iFlag = 1 : nFlag

	%------------------------------------%
	% Initial Fit                        %
	%------------------------------------%
	
		% Initial guess.
		initDist = squeeze( h(iComp,iFlag,:,if0) );
		[amp,imax] = max( initDist );

		% An initial guess using the curve fitting toolbox
		initDist = squeeze( h(iComp,iFlag,:,if0) );
		initFit  = fit( bins', initDist, 'gauss2' );
		p        = [ initFit.a1, initFit.b1, initFit.c1, initFit.a2, initFit.b2, initFit.c2 ];
	%------------------------------------%
	% Fit Forward from Fc                %
	%------------------------------------%

		% Loop forward to last frequency
		for ii = if0 : nFreq
			% Create a new function handle
			%   - Must do so to update the Y-data
			x        = bins;
			y        = double( squeeze(h(iComp,iFlag,:,ii))' );
			fbigauss = @(p) sum( ( y - ( p(1) * exp( -(x - p(2)).^2 / (2*p(3)^2) ) + ...
			                             p(4) * exp( -(x - p(5)).^2 / (2*p(6)^2) ) ) ).^2 );
			
			% Get the fit parameters
			[temp, ~, flag] = fminsearch( fbigauss, p, options );
	
			% The peaks occur at the mean values of the distribution
			h_floor(iComp,iFlag,ii) = min( temp( [2,5] ) );
			h_sig(iComp,iFlag,ii)   = max( temp( [2,5] ) );
			
			% Keep track of the number of times the max iteration was exceeded
			if flag == 0
				nMaxIter = nMaxIter + 1;
			end
		end
		
	%------------------------------------%
	% Fit Backward from Fc               %
	%------------------------------------%
		
		% Loop backward to first frequency
		p = [ initFit.a1, initFit.b1, initFit.c1, initFit.a2, initFit.b2, initFit.c2 ];
		for ii = if0 : -1 : ifc
			% Create a new function handle
			%   - Must do so to update the Y-data
			x        = bins;
			y        = double( squeeze(h(iComp,iFlag,:,ii))' );
			fbigauss = @(p) sum( ( y - ( p(1) * exp( -(x - p(2)).^2 / (2*p(3)^2) ) + ...
			                             p(4) * exp( -(x - p(5)).^2 / (2*p(6)^2) ) ) ).^2 );
			
			% Get the fit parameters
			[temp, ~, flag] = fminsearch( fbigauss, p, options );
	
			% The peaks occur at the mean values of the distribution
			h_floor(iComp,iFlag,ii) = min( temp( [2,5] ) );
			h_sig(iComp,iFlag,ii)   = max( temp( [2,5] ) );
			
			% Keep track of the number of times the max iteration was exceeded
			if flag == 0
				nMaxIter = nMaxIter + 1;
			end
			
			% Next iteration
			p = temp;
		end
		
		% Store the fit
		fit0{iComp, iFlag} = initFit;
	end
	end
	
	% Report
	if nMaxIter > 0
		mrfprintf( 'logwarn', 'Maximum number of iterations exceeded %i of %i times.', ...
		                      nMaxIter, nComp*nFlag*nFreq );
	end
end