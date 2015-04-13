%
% Name
%   mms_fg_calparams2matrix
%
% Purpose
%   Create an orthogonalization matrix from the fluxgate calibration
%   parameters.
%
% Calling Sequence
%   [SENSOR_MATRIX, M] = mms_fg_calparams2matrix(G, DPHI, DTHETA, U3)
%     Use the instrument gain correction G, angular offset of 123 from the
%     XYZ Y-axis DPHI, elevation angle of 12 above XY, and the (x,y)
%     components of a unit vector that project 3 onto Z, U3, to create a
%     rotation matrix from the orthongal XYZ magnetometer frame to the
%     nonorthogonal 123 sensor frame, M, and its inverse, SENSOR_MATRIX.
%
% Examples
%  Create a set of test data:
%    G      = [0.965819,     0.998249, 0.993915]';
%    dPhi   = [-0.0760805,  -0.0660412]';
%    dTheta = [-0.241677,   -0.0810649]';
%    u3     = [-0.00187507,  0.00423245]';
%    [sensor_matrix, m] = mms_fg_calparams2matrix(G, dPhi, dTheta, u3)
%    sensor_matrix =
%                   1.0354     -0.0011907   0.0019465
%                   0.0013123   1.0017     -0.0042374
%                   0.0042458   0.0014186   1.0061
%    m =
%                   0.96581     0.0011506  -0.0018637
%                  -0.0012825   0.99825     0.0042067
%                  -0.0040739  -0.0014124   0.9939
%
% Parameters
%   G             in, required, type=3xN float
%   dPhi          in, required, type=2xN float
%   dTheta        in, required, type=2xN float
%   u3            in, required, type=3xN float
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-03-21      Written by Matthew Argall
%
function [sensor_matrix, m] = mms_fg_calparams2matrix(G, dPhi, dTheta, u3)

	% Convert to degrees
	deg2rad = pi / 180.0;
	
	% Theta
	theta1 = (90.0 - dTheta(1, :)) .* deg2rad;
	theta2 = (90.0 - dTheta(2, :)) .* deg2rad;
	
	% Phi
	phi1 = dPhi(1, :) .* deg2rad;
	phi2 = (90.0 + dPhi(2, :)) .* deg2rad;
	
	a = u3(1, :);
	b = u3(2, :);
	
	% Allocate memory to output
	N = int64(length(G(:)) ./ 3);
	m = zeros(3, 3, N);
	
	% Row1
	m(1, 1, :) = G(1, :) .* sin(theta1) .* cos(phi1);
	m(2, 1, :) = G(1, :) .* sin(theta1) .* sin(phi1);
	m(3, 1, :) = G(1, :) .* cos(theta1);
	
	% Row2
	m(1, 2, :) = G(2, :) .* sin(theta2) .* cos(phi2);
	m(2, 2, :) = G(2, :) .* sin(theta2) .* sin(phi2);
	m(3, 2, :) = G(2, :) .* cos(theta2);
	
	% Row3
	m(1, 3, :) = G(3, :) .* a;
	m(2, 3, :) = G(3, :) .* b;
	m(3, 3, :) = G(3, :) .* sqrt(1 - (a.^2 + b.^2));
	
	% Orthogonalization matrix
	%   - Only calculate if you are asked.
	sensor_matrix = zeros(3, 3, N);
	for ii = 1 : N
		sensor_matrix(:, :, ii) = inv( m(:, :, ii) );
	end
end