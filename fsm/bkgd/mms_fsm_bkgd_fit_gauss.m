%
% Name
%   mms_fsm_bkgd_fit_gauss
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
function h_floor = mms_fsm_bkgd_fit_gauss(h, f, bins, p)

	% Initial guess
	if nargin < 4 || isempty(p)
		p = [ 30.0, -5.0, 1.0];
	end
	
	% Allocate memory
	%   - h = [comp, flag, bin, freq]
	dims    = size(h);
	nComp   = dims(1);
	nFlag   = dims(2);
	nBins   = dims(3);
	nFreq   = dims(4);
	h_floor = zeros(nComp, nFlag, nFreq);

	% Loop through all points
	for iComp = 1 : nComp
	for iFlag = 1 : nFlag
	for iFreq = 1 : nFreq
		% Create a new function handle
		%   - Must do so to update the Y-data
		x      = bins;
		y      = double( squeeze(h(iComp,iFlag,:,iFreq))' );
		fgauss = @(p) sum( ( y - ( p(1) * exp( -(x - p(2)).^2 / (2*p(3)^3) ) ) ).^2 );

		% Get the fit parameters
		temp = fminsearch( fgauss, p );

		% Peaks occur at the mean value of the distribution
		h_floor(iComp,iFlag,iFreq) = temp(2);
		
		% Next iteration
		p = temp;
	end
	end
	end
end