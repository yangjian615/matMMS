%
% Name
%   mms_edi_beam_width
%
% Purpose
%   Create a polar grid and return the cartesian coordinate of each
%   grid point.
%
% Method
%   We start by finding the semi-major and -minor axes of the uncertainty ellipse.
%   These are the theta-hat and phi-hat spherical coordinate system axes at the
%   location of the beam, with r-hat being in the direction of the beam.
%
%              | -sin(phi) |                | cos(theta) cos(phi) |
%      e_phi = |  cos(phi) |      e_theta = | cos(theta) sin(phi) |
%              |     0     |                |     -sin(theta)     |
%
%   Next, we find the B-perp direction by crossing the firing vector with the
%   magnetic field vector (assume B points into Quadrant II and the beam is
%   fired perpendicular to e_phi and e_theta -- the beam cross-section)::
%
%         B_perp = v x B
%
%   Finally, we find the beam uncertainty in BPP using the equation of an ellipse and
%   that of distance::
%
%         x^2 / dPhi^2  +  y^2 / dTheta^2 = 1
%         d = sqrt( x^2 + y^2 )
%
%   where dPhi and dTheta are the lengths of the semi-major and semi-minor axes, and
%   (x,y) is the location at which B_perp intersects the edge of the ellipse.
%   Solving the first equation for y^2, we get
%
%         y^2 = ( tan(alpha)/dPhi^2 + 1/dTheta^2 )^-1
%
%   where tan(alpha) = y / x, and alpha is the angle from the semi-major axes
%   to B_perp. Then, we note that the beam width, width = 2*d, factor out an y^2,
%   and make substitutions
%
%        width = 2 * y^2 * sqrt( 1 + x^2 / y^2 )
%        width = 2 * y^2 * sqrt( 1 + tan(alpha)^2 )
%        width = 2*sqrt( (1+tan(alpha)^2) / $
%                        ( (1/dTheta^2) + (tan(alpha)^2/dPhi^2) ) )
%
% Calling Sequence
%   SIGMA = mms_edi_beam_width(AZIMUTH, POLAR, B, GUN_ID)
%     Return the beam width (beam uncertainty) SIGMA along the direction
%     tangent to the plane perpendicular to the magnetic field B. The
%     beam is fired with an azimuthal and polar angle AZIMUTH and
%     POLAR, respectively. Each firing angle is associated with the gun
%     identified by its gun ID, GUN_ID (either 1 or 2). If GUN_ID is
%     a scalar, then all firing angles are assumed to be from the same
%     gun. There must be a 1-to-1 correspondence between firing angle
%     and magnetic field vector, and B must be in the GDU1 coordinate
%     system.
%
%   SIGMA = mms_edi_beam_width(FV, B, GUN_ID)
%     Provide the firing vectors FV in cartesian coordinates instead
%     of the gun firing angles. FV and B are expected to be in the
%     EDI1 instrument coordinate system.
%
%   [SIGMA, ALPHA] = mms_edi_beam_width(__)
%     Return the angle between the semi-major axis of the uncertainty
%     ellipse and unit vector perpendicular to the magnetic field.
%
%   [__] = mms_edi_beam_width(..., 'ParamName', ParamValue)
%     Provide any parameter name-value pair given below.
%
% Parameters
%   AZIMUTH         in, required, type=1xN float
%   POLAR           in, required, type=1xN float
%   FV              in, required, type=3xN float
%   B               in, required, type=3xN float
%   GUN_ID          in, required, type=1xN integer
%   'Radians'       in, optional, type=boolean, default=false
%                   Indicate that AZIMUTH and POLAR have units of
%                     radians, not degrees.
%   'SemiMajor'     in, optional, type=float
%                   The length of the semi-major axis of the beam
%                     uncertainty ellipse. If not provided, it is
%                     assumed to vary linarly from 0.5 to 2.0 as
%                     the azimuthal angle varies from 0 to 90 deg.
%   'SemiMinor'     in, optional, type=float
%                   The length of the semi-minor axis of the beam
%                     uncertainty ellipse. If not provided, it is
%                     assumed to vary linarly from 0.5 to 0.25 as
%                     the polar angle varies from 0 to 90 deg.
%
% Returns
%   SIGMA           out, required, type=1xN float
%   ALPHA           out, optional, type=1xN float
%
% Example
%  Compute the beam width for a given set of firing angles and magnetic field
%   >> azimuth = 35;
%   >> polar   = 52;
%   >> B       = [0; 1; 1];
%   >> gun_id  = 1;
%   >> [sigma, alpha] = mms_edi_beam_width(azimuth, polar, B, gun_id)
%     sigma = 0.023767838765260
%     alpha = 2.653537401507073
%
%  Compute the beam width for a given set of firing vectors and magnetic field
%   >> fv     = [0.6455;  0.4520;  0.6157];
%   >> B      = [0; 1; 1];
%   >> gun_id = 1;
%   >> [sigma, alpha] = mms_edi_beam_width(fv, B, gun_id)
%     sigma = 0.021347975231659
%     alpha = 2.619237582084609
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-08-18      Written by Matthew Argall
%   2015-08-28      If AZIMUTH and POLAR are given, rotate firing vectors
%                     from EDI2 to EDI1 coordinates. - MRA
%
function [sigma, alpha] = mms_edi_beam_width(arg1, arg2, arg3, arg4, varargin)
%------------------------------------%
% Required Parameters                %
%------------------------------------%

	% Figure out how the program was called
	szArg1 = size(arg1);
	szArg2 = size(arg2);
		
	% With FIRING VECTORS
	%   - It and the magnetic field should both be 3xN
	if szArg1(1) == 3 && szArg2(1) == 3
		fv      = arg1;
		azimuth = [];
		polar   = [];
		B       = arg2;
		gun_id  = arg3;
		
		if nargin > 3
			varargin = [ arg4 varargin ];
		end
	
	% With AZIMUTH and POLAR
	%   - Both should be column or row vectors
	elseif ( isrow(arg1) || iscolumn(arg1) ) && ( isrow(arg2) || iscolumn(arg2) )
		fv      = [];         % GDU firing vectors
		azimuth = arg1;       % Firing angle (azimuth)
		polar   = arg2;       % Firing angle (polar)
		B       = arg3;       % Magnetic field
		gun_id  = arg4;       % Gun ID
	else
		error( 'Unexpected calling sequence. Provide the firing angles or firing vectors.' )
	end

%------------------------------------%
% Defaults & Optional Parameters     %
%------------------------------------%
	% Defaults
	semi_major = [];
	semi_minor = [];
	tf_bcs     = false;
	tf_radians = false;
	
	% Optional arguments
	nOptArgs = length(varargin);
	for ii = 1 : 2 : nOptArgs
		switch varargin{ii}
			case 'Radians'
				tf_radians = varargin{ii+1};
			case 'SemiMajor'
				semi_major = varargin{ii+1};
			case 'SemiMinor'
				semi_minor = varargin{ii+1};
			otherwise
				if ischar(varargin{ii})
					error( ['Unexpected parameter name "' varargin{ii} '".'] );
				else
					error( ['Parameter name expected, not value of class ' class(varargin{ii})] )
				end
		end
	end
	
%------------------------------------%
% Check Inputs                       %
%------------------------------------%
	% Check sizes
	nGunID = length(gun_id);
	nAz    = length(azimuth);
	nPol   = length(polar);
	szB    = size(B);
	szFV   = size(fv);
	N      = max( [nAz nPol szFV(2) szB(2)] );
	
	% All must be same size
	assert( szFV(1) == 0 || szFV(1) == 3, 'FV must be 3xN.' );
	assert( isequal( szB, [3, N] ), 'B must be 3xN.' );
	assert( nAz == 0     || nAz     == N, 'AZIMUTH has incorrect number of elements.' );
	assert( nPol == 0    || nPol    == N, 'POLAR has incorrect number of elements.' );
	assert( nGunID == 1  || nGunID  == N, 'Incorrect number of Gun IDs.' );
	
	% Gun ID must be 1 || 2
	assert( isempty( setdiff( gun_id, [1 2] ) ), 'GUN_ID must be 1 or 2.' )
	
	% Expand the GunID
	if nGunID == 1
		gun_id = repmat(gun_id, 1, N);
	end
	
%------------------------------------%
% Firing Vectors & Angles            %
%------------------------------------%
	
	% Compute the firing vectors
	if isempty(fv)
		% Convert to radians
		if ~tf_radians
			deg2rad = pi / 180.0;
			azimuth = azimuth * deg2rad;
			polar   = polar   * deg2rad;
		end
		
		% Sines and Cosines
		cosAz  = cos(azimuth);
		cosPol = cos(polar);
		sinAz  = sin(azimuth);
		sinPol = sin(polar);
	
		% Cartesian components of firing direction
		fv      = zeros(3, N);
		fv(1,:) = sinPol .* cosAz;
		fv(2,:) = sinPol .* sinAz;
		fv(3,:) = cosPol;
		
		% Rotate Gun2 Firing vectors from EDI2 to EDI1 coordinates
		%          | 1  0  0 |
		%   edi1 = | 0 -1  0 | edi2
		%          | 0  0 -1 |
		fv(2:3, gun_id==2) = -fv(2:3, gun_id==2);

	% Compute the firing angles
	else
		% Normalize firing vector
		fv = mrvector_normalize(fv);

		% Azimuth firing angle
		%   - Careful of the poles
		azimuth  = zeros(1, N);
		iNotPole = find( abs( fv(1,:) ) >= 1.0d-8 | abs( fv(1,:) )  < 1.0d-8 );
		if ~isempty(iNotPole)
			azimuth(iNotPole) = atan2( fv(2, iNotPole),  fv(1, iNotPole) );
		end

		% Polar firing angle
		polar = acos( fv(2,:) );
		
		% Sines and Cosines
		cosAz  = cos(azimuth);
		cosPol = cos(polar);
		sinAz  = sin(azimuth);
		sinPol = sin(polar);
	end
	
%------------------------------------%
% Semi-Major and -Minor Axes         %
%------------------------------------%
	% Semi-minor axis in azimuth direction.
	%   - Varies from 0.5 to 0.25 as the azimuth angle increases from 0 to 90 degree
	if isempty(semi_minor)
		semi_minor = 0.5 * (1.0 - polar / pi );
	else
		semi_minor = repmat(semi_minor, N);
	end
	
	% Semi-major axis in polar direction
	%   - Varies from 0.5 to 2.0 as polar angle increases from 0 to 90 degrees
	if isempty(semi_major)
		semi_major = 0.5 * ( 1.0 + 6.0 * polar / pi );
	else
		semi_major = repmat( semi_major, N );
	end
	
%------------------------------------%
% Magnetic Field                     %
%------------------------------------%
	
	% Normalize the magnetic field
	b_hat = mrvector_normalize(B);
	
%------------------------------------%
% Beam Unit Vectors                  %
%------------------------------------%
	
	% Calculate the unit vector perpendicular to B and the gun firing vector.
	%   - e_perp = fv x B
	e_perp = mrvector_cross( fv, b_hat );
	e_perp = mrvector_normalize(e_perp);
	
	% Calculate the spherical polar unit vector measured along the beam
	e_pol = zeros(3, N);
	e_pol(1,:) =  cosPol .* cosAz;
	e_pol(2,:) =  cosPol .* sinAz;
	e_pol(3,:) = -sinPol;

	% Angle between the perpendicular and polar directions
	%   - Take the inner product of the two unit vectors
	alpha = acos( mrvector_dot(e_perp, e_pol) );
	
%------------------------------------%
% Beam Width                         %
%------------------------------------%
	
	% Calculate the full width of the beam along the
	% B-perp line that cuts through the beam profile.
	enum  = 1 + tan(alpha).^2;
	denom = (1.0 ./ semi_major).^2 + ( tan(alpha) ./ semi_minor ).^2;
	
	sigma = (pi/180.0 * 2.0) * sqrt( enum ./ denom );
end