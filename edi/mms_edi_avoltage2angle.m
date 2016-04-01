%
% Name
%   mms_edi_avoltage2angle
%
% Purpose
%   Convert EDI analog voltages to firing directions.
%
% Calling Sequence
%   FIRE_DIR = mms_edi_avoltage2angle(AVX, AVY)
%     Convert x- and y-components of the analog voltages
%     to a firing direction in spherical coordinates, ordered
%     ordered as [azimuth, elevation, r], where r is
%     the total voltage applied.
%
%   FIRE_DIR = mms_edi_avoltage2angle(AVX, AVY, 'Cartesian')
%     Return the firing direction in cartesian coordaintes,
%     [x, y, z]. The firing vector is not normalized.
%
% Parameters
%   AVX             in, required, type=1xN double
%   AVY             in, required, type=1xN double
%   'Cartesian'     in, optional, type=char
%
% Returns
%   FIRE_DIR        out, required, type=3xN double
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
function fire_dir = mms_edi_avoltage2angle(aVx, aVy, arg3)

	% Default
	tf_cart = false;

	if nargin > 2
		assert( ischar(arg3) && strcmp(arg3, 'Cartesian'), 'Unknown parameter in 3rd position.');
		tf_cart = true;
	end
	
	% Analog voltages must be converted from int16 to float
	% before coming here. sqrt() and atan() do not like integers.
	if isa(aVx, 'int16') || isa(aVy, 'int16')
		error( 'Analog voltages must be converted to floats before being passed in.' )
	elseif ~isfloat(aVx) || ~isfloat(aVy)
		error( 'Analog voltages must be floats.' )
	end

	%
	% MMS constants for calculation of angles from voltages
	%   vx    = aVx - 13824
	%   theta = 1/0.77 * asin( r / 14431.2 )
	%
	%--------------------
	%
	%  Analytic voltages used in the flight software have this relationship with gun
	%  firing angles:
	%
	%  vx = 14431.2 * cos(phi) * sin(0.77*theta)
	%  vy = 14431.2 * sin(phi) * sin(0.77*theta)
	%
	%  o) The numbers 14431.2 and 0.77 are the coefficients of this model.
	%  o) 13824 (0x3600) is a constant that is used to shift the bipolar analytic
	%  voltages vx,vy into a strictly positive range. These positive numbers are
	%  what is transmitted in telemetry. Therefore, to regain the original bipolar
	%  analytic voltages, the constant must be subtracted on the ground.
	%
	%  You may wonder why we are using a constant that is smaller than 14431.2 for
	%  shifting into a positive range. The reason this works is that we only go up
	%  to theta=95 degrees. With this restriction the sin() in the above equations
	%  does not get larger than sin(0.77 * 95) = 0.957, and therefore the max
	%  negative vx or vy we will encounter is -14431.2 * 0.957 = -13810.7
	%
	%  The choice of 13824 (0x3600) is not arbitrary. It has to do with lookup table
	%  grid spacing.
	%

	% Shift analytic voltages back to +/- values.
	Vx  = aVx - 13824.0;
	Vy  = aVy - 13824.0;
	
	% Radius
	%   - As read from the CDF file, analog voltages are int16 
	%   - sqrt does not like int16.
	r = sqrt( Vx.^2 + Vy.^2 );
	
	% Azimuthal Angle (Radians)
	theta = atan(Vy ./ Vx);
	
	% Elevation Angle (Radians)
	phi = (1.0 ./ 0.77) .* asin( r ./ 14431.2 );

	% Cartesian
	if tf_cart
		[x, y, z] = cart2sph(theta, phi, r);
		fire_dir  = [x; y; z];
	else
		fire_dir  = [theta; phi; r];
	end
end