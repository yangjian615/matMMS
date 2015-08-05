%
% Name
%   mms_fg_h_stemp_etemp
%
% Purpose
%   Calculate the forward gain H(stemp, etemp) from counts to nT, given ground
%   calibration constants from AFG or DFG. Applied via multiplication after offsets
%   have been removed.
%
% Calling Sequence
%   FGAIN = mms_fg_calparams2matrix(STEMP, ETEMP)
%     Calculate the forward gain FGAIN of the instrument from the sensor
%     and electronics temperatures, STEMP and ETEMP. The resulting matrix
%     is applied via multiplication in the sensor coordinate system after
%     offsets have been removed.
%
%   CDATA = mms_fg_calparams2matrix(..., DATA)
%     If given uncorrected data, DATA, the temperature correction factors will
%     be calculated and applied to them. The corrected data, CDATA, is returned.
%
% Parameters
%   STEMP         in, required, type=scalar or 1xN float
%   ETEMP         in, required, type=scalar or 1xN float
%   DATA          in, optional, type=3x1 or 3xN float
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-07-27      Written by Matthew Argall
%
function fgain = mms_fg_h_stemp_etemp(stemp, etemp, cal_const, data)

	% Defaults
	if nargin < 4
		data = [];
	end
	
	% Check Inputs
	ns = length(stemp);
	ne = length(etemp);
	n  = max( [ns ne] );
	assert( ns == n | ns == 1, 'STEMP must be scalar or same length as ETEMP.');
	assert( ne == n | ne == 1, 'ETEMP must be scalar or same length as ETEMP.');

	% Allocate memory
	if isempty( data )
		fgain = ones(3, n);
	else
		sz = size( data );
		assert( isequal( sz, [3 1] ) | isequal( sz, [1 n] ), 'DATA must be 3x1 or 3xN, where N is the length of STEMP.' );
		fgain = double( data );
	end

	% Calculate the forward gain
	for ii = 1 : 3
		fgain(ii,:) = fgain(ii,:) * cal_const.rsf * cal_const.psf(ii)           .* ...
		              ( cal_const.h_stemp(ii) + cal_const.m_stemp(ii) * stemp ) .* ...
		              ( 1 + cal_const.m_etemp(ii) * etemp );
	end
end