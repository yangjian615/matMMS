%
% Name
%   mms_edi_polar_grid
%
% Purpose
%   Create a polar grid and return the cartesian coordinate of each
%   grid point.
%
% Calling Sequence
%   GRID = mms_edi_polar_grid(R, PHI)
%     Create a structure array GRID containing the x- and y-
%     coordinates of each grid point on a polar grid with 
%     radial extent R and polar angle extent PHI. R and PHI
%     should be 3-element vectors containing the minimum,
%     maximum, and spacing of the grid in each dimension:
%     [r_min, r_step, r_max], [phi_min, phi_step, phi_max].
%
% Parameters
%   R               in, required, type=1x3 double
%   PHI             in, required, type=1X3 double
%
% Returns
%   GRID            out, required, type=structure
%                   Fields are:
%                     'x' - x-coordinates of grid points
%                     'y' - y-coordinates of grid points
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-08-14      Written by Matthew Argall
%
function grid = mms_edi_polar_grid(r, phi)
	% Create the grid points
	phi_grid = phi(1) : phi(2) : phi(3);
	r_grid   =   r(1) :   r(2) :   r(3);
	
	% Calculate the cartesian coordinates throughout the grid
	%   - Create the grid by taking the inner product
	%     along the shallow dimension.
	deg2rad = pi / 180.0;
	x = r_grid' * cos( phi_grid * deg2rad );
	y = r_grid' * sin( phi_grid * deg2rad );
	
	% Create a structure array of grid points.
	grid = struct( 'x', num2cell(x), ...
	               'y', num2cell(y) );
end