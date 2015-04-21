%
% Name
%   mms_edi_bavg
%
% Purpose
%   Create a transformation matrix that will transform a coordinate system
%   into one in which the z-axis is along +/-B, oriented so that the new z-axis
%   points toward +z in the original coordinate system.
%
% Calling Sequence
%   xyz2bpp = mms_edi_xxyz2bpp(B)
%     Use magnetic field vectors to define a coordinate system transformation
%     defined by the b-hat normal direction.
%
% Parameters
%   B               in, required, type=3xN
%
% Returns
%   XYZ2BPP         out, required, type=3x3xN
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-20      Written by Matthew Argall
%
function xyz2bpp = mms_edi_xxyz2bpp(b)

	% Coordinate axes directions of bpp relative to xyz
	%   - Normalize the magnetic field to create z-hat
	z_hat = mrvector_normalize(b);

	% Ensure z-hat points up with respect to xyz
	ineg = find( z_hat(3,:) < 0 );
	if ~isempty(ineg)
		z_hat(:,ineg) = -z_hat(:,ineg);
	end
	
	% Create X and Y
	x_hat = mrvector_cross([0 1 0], z_hat);
	x_hat = mrvector_normalize(x_hat);
	y_hat = mrvector_cross(z_hat, x_hat);

	% Build matrix
	xyz2bpp        = zeros( [3 size(b)] );
	xyz2bpp(1,:,:) = x_hat;
	xyz2bpp(2,:,:) = y_hat;
	xyz2bpp(3,:,:) = z_hat;
end