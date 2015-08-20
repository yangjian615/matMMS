%
% Name
%   mms_edi_cost_function
%
% Purpose
%   Compute the cost function of a set of EDI beams.
%
% Calling Sequence
%   COSTFN = mms_edi_cost_function(FV_BPP, POS_BPP, BEAM_WIDTH, GRID)
%     Compute the cost function COSTFN of a set of EDI firing vectors
%     FV_BPP with width BEAM_WIDTH fired from guns at positions POS_BPP,
%     all in the plane perpendicular to the magnetic field (BPP). The
%     cost function is with respect to a set of gridded points GRID,
%     also in BPP.
%
% Parameters
%   FV_BPP          in, required, type=1xN float
%   POS_BPP         in, required, type=1xN float
%   BEAM_WIDTH      in, required, type=3xN float
%   GRID            in, required, type=3xN float
%
% Returns
%   COSTFN          out, required, type=NxM float
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-08-18      Written by Matthew Argall
%
function costFn = mms_edi_cost_function(fv_bpp, pos_bpp, beam_width, grid)
	
	% Check sizes
	N     = length(beam_width);
	szFV  = size(fv_bpp);
	szPos = size(pos_bpp);
	assert( szFV(2) == N && szPos(2) == N, ...
	        'Incompatible number of firing vectors, positions, and beam widths.')

	% Get grid dimensions and allocate memory
	%   - First dimension is radius, second is azimuth
	dims   = size(grid);
	costFn = zeros( dims );
	
	% Calculate firing angles in BPP
	firing_angle = atan2( fv_bpp(2,:), fv_bpp(1,:) );
	
	% Loop through all beams, calculating their contribution
	% to the cost function at each grid point.
	for ii = 1 : N
		% Angle from gun position to each grid point
		grid_angle = atan2( grid.y - pos_bpp(2,ii), grid.x - pos_bpp(1,ii) );
		
		% Difference between firing angle and angle to grid point
		dAngle = grid_angle - firing_angle(ii);
	
		%
		% Re-normalize angular difference
		% Reasons (numbers are given in degrees not radians, for clarity)
		% 1) the input range of -360 to +360 degrees needs to be reduced
		%    to -180 to 180 (add 360 if below -180, subtract 360 if above 180)
		% 2) the angular deviation between two directions cannot be larger
		%    than 90 degrees (in other words: no distinction between parallel
		%    and anti-parallel, since we are only interested in orientation,
		%    not in direction)
		% The two transformations combined are summarized in the following
		% table:
		%       IN                  OUT
		%   -360 ... -270         0 ... 90
		%   -270 ... -180       -90 ...  0
		%   -180 ...  -90         0 ... 90
		%    -90 ...    0       -90 ...  0
		%      0 ...   90         0 ... 90
		%     90 ...  180       -90 ...  0
		%    180 ...  270         0 ... 90
		%    270 ...  360       -90 ...  0
		%
		% The formula below achieves this transformation
		%
		dAngle = mod( dAngle + 2.5*pi, pi) - pi/2;
	
		% Add contribution of current beam to cost function
		costFn = costFn + ( dAngle ./ beam_width(ii) ).^2;
	end
end