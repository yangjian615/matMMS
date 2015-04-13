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
%   BSCS2DESPUN = mms_fdoa_xbcs2despun(SC, TSTART, TEND, TYPE, TIME, ATT_DIR)
%     Create a transformation matrix BSCS2DESPUN using attitude
%     data from MMS satellite SC in the time interval [TSTART,TEND].
%     Despin using the phase angle of type TYPE and interpolate to
%     times TIME. Attitude data can be found in directory ATT_DIR.
%
%   [..., PHASE] = mms_fdoa_xbcs2despun(__)
%     Also return the interpolated, unwrapped phase angle PHASE.
%
%   [..., ATTITUDE] = mms_fdoa_xbcs2despun(__)
%     Also return all of the attitude data.
%
% Parameters
%   SC              in, required, type=char
%   TSTART          in, required, type=char
%   TEND            in, required, type=char
%   TYPE            in, required, type=char
%   TIME            in, required, type=1xN int64 (cdf_tt2000)
%   ATT_DIR         in, required, type=char
%
% Returns
%   BCS2DESPUN      out, required, type=3x3xN double
%   PHASE           out, required, type=1xN double
%   ATTITUDE        out, required, type=struct
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
%
function [bcs2despun, phase, attitude] = mms_fdoa_xdespin(sc, tstart, tend, type, time, att_dir)
	
	% Read the FDOA data
	attitude = mms_fdoa_read_defatt(sc, tstart, tend, att_dir);
	
	% Interpolate the phase
	%   - Returned in degrees
	DEG2RAD = pi / 180.0;
	phase   = mms_fdoa_interp_phase(attitude, time, type) * DEG2RAD;
	
	% Sine and cosine of phase;
	cosPhase = cos(phase);
	sinPhase = sin(phase);
	
	% Create the rotation matrix
	bcs2despun        =  zeros(3, 3, length(time));
	bcs2despun(1,1,:) =  cosPhase;
	bcs2despun(2,1,:) =  sinPhase;
	bcs2despun(1,2,:) = -sinPhase;
	bcs2despun(2,2,:) =  cosPhase;
	bcs2despun(3,3,:) =  1;
end