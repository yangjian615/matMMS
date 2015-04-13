%
% Name
%   mms_fdoa_interp_phase
%
% Purpose
%   Interpolate phase.
%
% Calling Sequence
%   PHASEOUT = mms_fdoa_interp_phase(DEFATT, TIMEOUT)
%     Interpolate the phase of the Sun-to-body-X dihedral angle about
%     major principal axis P. DEFATT is a structure returned by
%     mms_fdoa_read_defatt.m. TIMEOUT are CDF_TIME_TT2000 values at
%     which to interpolate the phase.
%
%   PHASEOUT = mms_fdoa_interp_phase(DEFATT, TIMEOUT, TYPE)
%     Specify the phase type: 'L', 'P', 'w', or 'Z'.
%       'L' = Sun-to-body-X dihedral angle about angular momentum vector L.
%       'P' = Sun-to-body-X dihedral angle about major principal axis P.
%       'w' = Sun-to-body-X dihedral angle about rotation rate vector.
%       'z' = Sun-to-body-X dihedral angle about body Z-axis
%
% Parameters
%   DEFATT          in, required, type=structure
%   TIMEOUT         in, optional, type=1xN int64 (cdf_tt2000)
%   TYPE            in, optional, type=char
%
% Returns
%   PHASEOUT        out, required, type=1xM double
%
% See Also
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-11      Written by Matthew Argall
%
function phaseOut = mms_fdoa_interp_phase(defatt, timeOut, type)
	
	if nargin < 3
		type = 'L';
	else
		assert( ismember(type, {'L', 'P', 'w', 'z'}), 'Type must be {"L" | "P" | "w" | "z"}.' )
	end
	
	% Convert TAI to TT2000
	att_tt2000 = mms_fdoa_epoch2tt2000(defatt.('TAI'), 'AttTAI', true);
	
	% Check if there are leaps in the phase.
	%   - w(:,4) is the spin phase in radians/second
	%   - If the spin phase changes by more than 180 degrees, MrPhaseUnwrap will not know
	%     which direction the phase jumped.
%	nGaps = sum( diff(att_tt2000) > ( 175.0 * 1.0d9 / defatt.('w')(:,4) ) )
%	if nGaps > 0
%		warning('MMS_FDOA_InterpPhase:Gaps', 'Gaps found while interpolating phase.')
%	end

	% Start by unwrapping the phase
	phunwrap = MrPhaseUnwrap(defatt.(type)(:,3)', 360.0);
	
	% Inteprolate the phase
	%   - First, convert to double
	t0 = att_tt2000(1);
	phaseOut = interp1( double(att_tt2000 - t0), phunwrap, double(timeOut - t0) );
end