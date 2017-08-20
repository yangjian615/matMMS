%
% Name
%   mms_fdoa_xbcs2smpa
%
% Purpose
%   Create a coordinate transformation from the spinning spacecraft
%   body coordinate system (BCS) to the spinning major priciple axis
%   coordinate system (SMPA).
%
% Calling Sequence
%   BCS2SMPA = mms_fdoa_xbcs2smpa(MPA)
%     Use the z-MPA axis extracted from the definitive attitude data
%     header to create a rotation matrix from SMPA to BCS.
%
% Parameters
%   MPA             in, required, type=1x3 double
%
% Returns
%   BCS2SMPA        out, required, type=3x3 double
%
% See Also
%   mms_fg_interp_cal.m    --   FGM reports & interpolates its own zMPA
%   mms_fdoa_read_defatt.m
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-11      Written by Matthew Argall
%
function bcs2smpa = mms_fdoa_xbcs2smpa(mpa)

	% Each FDOA file has a single z-MPA vector.
	%   - We expect it to be constant.
	%   - If not, we should interpolate (slerp)
	assert( ( isrow(mpa) || iscolumn(mpa) ) && length(mpa) == 3, ...
	        'DefAtt MPA vector must be 1x3 or 3x1.' )
	
	% Create unit vectors
	%   - To deal with arrays of vectors, see:
	%       mrvector_cross.m
	%       mrvector_magnitude.m
	%       mrvector_normalize.m
	smpaz = mpa / norm(mpa, 2);
	smpay = cross(smpaz, [1 0 0]);
	smpay = smpay / norm(smpay, 2);
	smpax = cross(smpay, smpaz);
	
	% Create the rotation matrix
	bcs2smpa      = zeros(3, 3);
	bcs2smpa(:,1) = smpax;
	bcs2smpa(:,2) = smpay;
	bcs2smpa(:,3) = smpaz;
end