%
% Name
%   mms_rot_despun2gsm
%
% Purpose
%   Transform data from a spinning coordinate system to GSM.
%
%   IGRF coordinates of Earth's dipole have been hard-coded into the program
%   and are valid between the dates 2015 and 2019.
%
% Calling Sequence
%   V_GSM = mms_rot_despun2gsm(T, V, DEFATT)
%       Rotate the vector(s) V measured at times T from DMPA to GSM
%       coordinates using the major principal axis found in the
%       structure DEFATT. DEFATT is returned by the program
%       mms_fdoa_read_defatt.m
%
%   V_GSM = mms_rot_despun2gsm(T, V, DEFATT, TYPE)
%       If V is not in DMPA coordinates, but some other despun spacecraft
%       coordinate system, then specify which axis is the despun z-axis
%       with the TYPE parameter. Choices are "L", "P", "W", or "Z".
%
%   [..., V_GSE] = mms_rot_despun2gsm(__)
%       Also return the vector in GSE coordinates.
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
function [v_gsm, v_gse, v_gei] = mms_rot_despun2gsm(t, v, defatt, type)
	
	% Default
	if nargin < 4
		type = 'P';
	end
	
	% Transformation matrix
	gei2despun = mms_fdoa_xgei2despun(defatt, t, type);
	despun2gei = permute(gei2despun, [2,1,3]);
	v_gei      = mrvector_rotate(despun2gei, v);

	% Rotate to GSE
	[mjd, utc] = tt2000toMJDutc(t');
	xgei2gse   = gei2gse(mjd, utc);
	v_gse      = mrvector_rotate(xgei2gse, v_gei);

	% Make sure all years are 2015
	year = MrCDF_Epoch_Breakdown(t);
	assert( sum(year(:,1) < 2015 & year(:,1) >= 2020) == 0, 'Update IGRF coefficients. Time must be from 2015 - 2019.');
	
	% GSE --> GSM
	g10        = -29442.0;
	g11        =  -1501.0;
	h11        =   4797.1;
	xgse2gsm   = gse2gsm(g10, g11, h11, mjd, utc);
	v_gsm      = mrvector_rotate(xgse2gsm, v_gse);
end

