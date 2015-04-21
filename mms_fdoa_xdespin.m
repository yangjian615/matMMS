%
% Name
%   mms_fdoa_xbcs2despun
%
% Purpose
%   Create a coordinate transformation from the spinning spacecraft
%   body coordinate system (BCS) to a despun coordinate system.
%
%   Process:
%     1. Find and read attitude data
%     2. Interpolate the phase
%     3. Build transformation matrix.
%
% Calling Sequence
%   BSCS2DESPUN = mms_fdoa_xbcs2despun(ATTITUDE, TIME)
%     Create a transformation matrix BSCS2DESPUN using an attitude
%     data structure returned by mms_fdoa_read_defatt.m ATTITUDE.
%     Despin using the phase angle of the angular momentum system, 'L'.
%
%   BSCS2DESPUN = mms_fdoa_xbcs2despun(..., TYPE)
%     Despin using the phase angle of type TYPE. Options are:
%       'L' - angular momentum
%       'P' - major principal axis of inertia
%       'w' - spin vector
%       'z' - body z-axis
%
%   BSCS2DESPUN = mms_fdoa_xbcs2despun(..., 'SpinUp')
%     Spin up the data instead of despinning.
%
% Parameters
%   ATTITUDE        in, required, type=struct
%   TIME            in, required, type=1xN int64 (cdf_tt2000)
%   TYPE            in, optional, type=char, default='L'
%   INSTR           in, optional, type=char, default='BCS'
%   'SpinUp',       in, optionanl, type=char
%
% Returns
%   BCS2DESPUN      out, required, type=3x3xN double
%
% See Also
%   mms_fdoa_read_defatt.m
%   mms_fdoa_interp_phase.m
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-11      Written by Matthew Argall
%   2015-04-14      Do not read data here. - MRA
%   2015-04-18      Added the 'SpinUp' parameter. - MRA
%
function bcs2despun = mms_fdoa_xdespin(attitude, time, type, instr, arg5)

	% Inputs
	tf_spinup = false;
	if nargin < 3
		type = 'L';
	end
	if nargin < 4
		instr = 'BCS';
	end
	if nargin > 4
		assert( ischar(arg5) && strcmp(arg5, 'SpinUp'), 'Unexpected parameter in position 5');
		tf_spinup = true;
	end
	DEG2RAD = pi / 180.0;
	
	%
	% The spin phase reported by FDOA is relative to the physical
	% BCS x-axis. In order to despin instrument data properly, we
	% must offset the FDOA phase by the phase difference between
	% when the instrument x-axis points along the s/c-sun line and
	% when x-BCS points along the s/c-sun line.
	%
	% If the instrument is in BCS already, no phase offset needs to
	% be applied. The offset is the angle between the original CS
	% x-axis and the x-axis of BCS.
	%
	offset = 0;

	% Interpolate the phase
	%   - Returned in degrees
	phase   = mms_fdoa_interp_phase(attitude, time, type) * DEG2RAD;
	
	% If we are spinning up the data, we have to invert the phase
	if tf_spinup
		phase = -phase;
	end
	
	% Sine and cosine of phase;
	cosPhase = cos(phase + offset);
	sinPhase = sin(phase + offset);
	
	% Create the rotation matrix
	bcs2despun        =  zeros(3, 3, length(time));
	bcs2despun(1,1,:) =  cosPhase;
	bcs2despun(2,1,:) =  sinPhase;
	bcs2despun(1,2,:) = -sinPhase;
	bcs2despun(2,2,:) =  cosPhase;
	bcs2despun(3,3,:) =  1;
end