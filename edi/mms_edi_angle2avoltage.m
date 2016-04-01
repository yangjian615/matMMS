%
% Name
%   mms_edi_angle2avoltage
%
% Purpose
%   Convert EDI firing angles to analog voltages.
%
% Calling Sequence
%   AVOLTAGE = mms_edi_avoltage2angle(SOLID_ANGLE)
%     Convert a solid angle to an analog voltage. Solid
%     angle should be [azimuth, elevation]. The analog
%     voltage is the voltage that needs to be applied in
%     the [x, y] directions to achieve the desired firing
%     angles.
%
%   FIRE_DIR = mms_edi_avoltage2angle(SOLID_ANGLE, 'Degrees')
%     Indicate the SOLID_ANGLE is given in degrees.
%
% Parameters
%   SOLID_ANGLE     in, required, type=2xN double
%   'Degrees'       in, optional, type=char
%
% Returns
%   AVOLTAGE        out, required, type=2xN double
%
% Examples
%   Convert firing angles to analog voltages and back.
%     >> avoltage    = mms_edi_angle2avoltage([45; 75], 'Degrees');
%     >> solid_angle = mms_edi_avoltage2angle( avoltage(1), avoltage(2) );
%     >> solid_angle(1:2) * 180/pi
%        solid_angle = 45.00   74.9999
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-18      Written by Matthew Argall
%
function avoltage = mms_edi_angle2avoltage(solid_angle, arg2)
	
	% Degrees given?
	tf_degrees = false;
	if nargin > 1
		assert( ischar(arg2) && strcmp(arg2, 'Degrees'), 'Unknown parameter in 2nd position.');
		tf_degrees = true;
	end
	
	% Convert to radians?
	if tf_degrees
		theta = deg2rad( solid_angle(1,:) );
		phi   = deg2rad( solid_angle(2,:) );
	else
		theta = solid_angle(1,:);
		phi   = solid_angle(2,:);
	end

	% 
	%  MMS constants for calculation of angles from voltages
	%    vx    = aVx - 13824
	%    vy    = aVy - 13824
	% 
	%    vx = 14431.2 * cos(theta) * sin(0.77*phi)
	%    vy = 14431.2 * sin(theta) * sin(0.77*phi)
	% 
	%    r = vx^2 + vy^2
	%      = sqrt( 14431.2^2 * cos(theta)^2 * sin(0.77*phi)^2 + 
	%              14431.2^2 * sin(theta)^2 * sin(0.77*phi)^2 )
	%      = sqrt( 14431.2^2 * sin(0.77*phi)^2 * 
	%              (cos(theta)^2 + sin(theta)^2) )
	%      = sqrt( 14431.2^2 * sin(0.77*phi)^2 )
	%      = 14431.2 * sin(0.77*phi)
	% 
	%    phi   = 1/0.77 * asin(r / 14431.2)
	%    theta = atan(vy / vz)
	% 
	% --------------------
	% 
	%   Analytic voltages used in the flight software have this relationship with gun
	%   firing angles:
	% 
	%    vx = 14431.2 * cos(theta) * sin(0.77*phi)
	%    vy = 14431.2 * sin(theta) * sin(0.77*phi)
	% 
	%   o) The numbers 14431.2 and 0.77 are the coefficients of this model.
	%   o) 13824 (0x3600) is a constant that is used to shift the bipolar analytic
	%   voltages vx,vy into a strictly positive range. These positive numbers are
	%   what is transmitted in telemetry. Therefore, to regain the original bipolar
	%   analytic voltages, the constant must be subtracted on the ground.
	% 
	%   You may wonder why we are using a constant that is smaller than 14431.2 for
	%   shifting into a positive range. The reason this works is that we only go up
	%   to theta=95 degrees. With this restriction the sin() in the above equations
	%   does not get larger than sin(0.77 * 95) = 0.957, and therefore the max
	%   negative vx or vy we will encounter is -14431.2 * 0.957 = -13810.7
	% 
	%   The choice of 13824 (0x3600) is not arbitrary. It has to do with lookup table
	%   grid spacing.
	% 

	avx = 14431.2 .* cos(theta) .* sin(0.77 .* phi) + 13824;
	avy = 14431.2 .* sin(theta) .* sin(0.77 .* phi) + 13824;

	avoltage = [avx; avy];
end


