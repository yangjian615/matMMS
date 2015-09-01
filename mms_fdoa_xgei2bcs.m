%
% Name
%   mms_fdoa_xgei2bcs
%
% Purpose
%   Interpolate the quaternion rotation from GEI2000 to BCS.
%
% Calling Sequence
%   QBSCS2GEI = mms_fdoa_xgei2bcs(T_ATT, Q, TIME)
%     Interpolate the quaternians Q at temporal grid points T_ATT
%     onto the gridpoints TIME. Quaternions are measured in the
%     GEI inertial frame and the vector to be rotation should be
%     in the spacecraft body coordinate system (BCS). Rotatetion
%     by Q, then, rotates from GEI2000 to BCS.
%
%   QGEI2BCS = mms_fdoa_xgei2bcs(T_ATT, Q, TIME, 'inverse')
%     Return the quaterions that perform the inverse transformation,
%     from BCS to GEI2000.
%
% Parameters
%   T_ATT           in, required, type=int64 (cdf_time_tt2000)
%   Q               in, required, type=4xN float
%   TYPE            in, required, type=int64 (cdf_time_tt2000)
%   'inverse'       in, optional, type=char
%
% Returns
%   Q               out, required, type=4xN float
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-08-13      Written by Matthew Argall
%
function Q = mms_fdoa_xgei2bcs(t_att, q, time, arg4)

	% Defaults
	tf_inverse = false;

	% Inverse transformation
	if nargin == 4
		assert( ischar(arg4) && strcmp(arg4, 'inverse'), 'Unexpected fourth parameter.' );
		tf_inverse = true;
	end

	% Spacecraft attitude (quaternion). 
	% NOTE: this has noise at level of
	% 0.023 deg for spin phase
	% 0.008 deg for x-y axis

	% Interpolate quaternions and rotate input vectors to J2000 
	% Quaternions are J2000 to BCS; qtvrot convention is to use invert=1 to transform J2000 to BSC.
	% DON'T FORGET TO USE /SLERP!!
	% TODO: watch out for data gaps > 0.5 spin period.
	
	% Times cannot be int64, so convert to SSE
	%   - Reference epoch is the first attitude time stamp
	t0 = t_att(1);
	
	% Inteprolate the quaternions
	%   - Use SLERP for smoother interpolation.
	Q = mrquaternion_interp( double(t_att - t0), q, double(time - t0), 'slerp' );

	% Invert
	%   - DEFATT quaternions transform from BCS to J2000
	if tf_inverse
		Q(1:3,:) = -Q(1:3,:);
	end
end