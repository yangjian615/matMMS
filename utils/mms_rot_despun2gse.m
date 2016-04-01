%
% Name
%   mms_rot_despun2gse
%
% Purpose
%   Transform data from a spinning coordinate system to GSM.
%
% Calling Sequence
%   V_GSE = mms_rot_despun2gsm(T, V, DEFATT)
%       Rotate the vector(s) V measured at times T from DMPA to GSE
%       coordinates using the major principal axis found in the
%       structure DEFATT. DEFATT is returned by the program
%       mms_fdoa_read_defatt.m
%
%   V_GSE = mms_rot_despun2gsm(T, V, DEFATT, TYPE)
%       If V is not in DMPA coordinates, but some other despun spacecraft
%       coordinate system, then specify which axis is the despun z-axis
%       with the TYPE parameter. Choices are "L", "P", "W", or "Z".
%
%   [..., V_GEI] = mms_rot_despun2gsm(__)
%       Also return the vector in GEI coordinates.
%
% Parameters
%   T            in, required, type=int64 (cdf_time_tt2000)
%   V            in, required, type=3xN float
%   DEFATT       in, required, type=struct
%   TYPE         in, optional, type=char, default='P'
%
% See Also
%    mms_fdoa_read_defatt
%    mms_fdoa_xgei2despun
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-10-08      Written by Matthew Argall
%
function [v_gse, v_gei] = mms_rot_despun2gse(t, v, defatt, type)
	
	% Default
	if nargin < 4
		type = 'P';
	end
	
	% Transformation matrix
	gei2despun = mms_fdoa_xgei2despun(defatt, t, type)
	dmpa2gei   = permute(gei2dmpa, [2,1,3]);
	v_gei      = mrvector_rotate(despun2gei, v);

	% Rotate to GSE
	[mjd, utc] = tt2000toMJDutc(t');
	xgei2gse   = gei2gse(mjd, utc);
	v_gse      = mrvector_rotate(xgei2gse, b_gei);
end

