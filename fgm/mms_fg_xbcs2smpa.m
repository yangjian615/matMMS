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
% See Also
%   mrvector_cross.m
%   mrvector_normalize.m
%
% History:
%   2015-04-13      Written by Matthew Argall
%
function bcs2smpa = mms_fg_xbcs2smpa(mpa)

	% Create unit vectors
	smpaz = mrvector_normalize(mpa);
	smpay = mrvector_cross(smpaz, [1 0 0]);
	smpay = mrvector_normalize(smpay);
	smpax = mrvector_cross(smpay, smpaz);

	% Create the rotation matrix
	bcs2smpa        = zeros( [3 size(mpa)] );
	bcs2smpa(:,1,:) = smpax;
	bcs2smpa(:,2,:) = smpay;
	bcs2smpa(:,3,:) = smpaz;
end