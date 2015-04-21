%
% Name
%   mms_fg_xomb2smpa
%
% Purpose
%   Create a coordinate system transformation matrix that rotates
%   the spinning major principal axis (SMPA) system into the
%   orthogonal magnetometer system (OMB).
%
% Calling Sequence
%   OMB2SMPA = mms_fg_xomb2smpa()
%     Return the transformation from OMB to SMPA.
%
% Returns
%   OMB2SMPA        out, required, type=3x3 double
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-13      Written by Matthew Argall
%   2015-04-21      Needed to rotate cw (-225), not ccw. - MRA
%
function omb2smpa = mms_fg_xomb2smpa()

	%
	% From 7.2.2.5 of the MMS coordinate system handbook
	%
	%   The OMB coordinate system is defined to be a fixed 225°
	%   rotation from the SMPA coordinate system, as in Equation
	%   7.2-11. As described in Table 7.2-6, the OMB Z-axis is
	%   aligned with the MPA, while the OMB X-axis is closely 
	%   aligned with the AFG Boom X-axis. Because the OMB X, Y,
	%   and Z-axes are closely aligned with sensor axes 1, 2, and
	%   3, respectively, this coordinate system is useful as an 
	%   intermediate state of calibrated, orthogonalized AFG & DFG
	%   data before rotation into body coordinates.
	% 

	% Constants
	theta    = -225.0 * pi / 180.0;
	cosTheta = cos(theta);
	sinTheta = sin(theta);
	
	% Rotation matrix
	omb2smpa = [  cosTheta  sinTheta  0;
	             -sinTheta  cosTheta  0;
	                 0         0      1 ];
end