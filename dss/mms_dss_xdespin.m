%
% Name
%   mms_dss_xdespin
%
% Purpose
%   Crate an array of matrices that will transform from a spinning to a
%   not-spinning coordinate system.
%
% Calling Sequence
%   SPIN2DESPUN = mms_dss_xdespin( SUNPULSE, TIMES )
%     Take a structure returned by mms_dss_read_sunpulse, determine
%     the spin phase, interpolate to times TIMES, and create a
%     set of transformation matrices to go from the spinning
%     coordinate system to a despun coordinate system.
%
% Parameters
%   SUNPULSE        in, required, type=struct
%   TIMES           in, required, type=int64 (cdf_time_tt2000)
%
% Returns
%   SPIN2DESPUN     out, required, type=3x3xN double
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-19      Written by Matthew Argall
%
function spin2despun = mms_dss_xdespin( sunpulse, times, tf_spinup )

	if nargin < 3
		tf_spinup = false;
	end

	% Build matrix
	spin_phase = mms_dss_sunpulse2phase( sunpulse, times ) * (pi/180.0);

	% The spin phase is the angle of rotation away from the s/c-sun
	% line. To despin, we want to rotate against the phase. To spin-
	% up, we rotate with the phase.
	if ~tf_spinup
		spin_phase = -spin_phase;
	end

	% The DSS requires a -76 degree rotation about the z-axis to align
	% with BCS. Aligning BCS with the sun sensor requires a +76 degree
	% rotation.
	offset = 76 * pi/180;

	% Sine and Cosine of phase
	sinPhase = sin(spin_phase + offset);
	cosPhase = cos(spin_phase + offset);
	
	% Create the rotation matrix
	%                 |  cos  sin  0  |
	%   spin2despun = | -sin  cos  0  |
	%                 |   0    0   1  |
	spin2despun        =  zeros(3, 3, length(times));
	spin2despun(1,1,:) =  cosPhase;
	spin2despun(1,2,:) =  sinPhase;
	spin2despun(2,1,:) = -sinPhase;
	spin2despun(2,2,:) =  cosPhase;
	spin2despun(3,3,:) =  1;
end