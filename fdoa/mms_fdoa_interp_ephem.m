%
% Name
%   mms_fdoa_interp_ephem
%
% Purpose
%   Interpolate ephemeris positions.
%
% Calling Sequence
%   R_OUT = mms_fdoa_interp_ephem( T, R, T_OUT )
%     Interpolate ephemeris positions R measured at time T. The
%     result is an array of position vectors R_OUT at times T_OUT.
%
%   R_OUT = mms_fdoa_interp_ephem( ..., V )
%     Provide the velocity vectors V at times T. Used during
%     interpolation. If not given, V will be calculated directly
%     from R and T.
%
%   [R_OUT, V_OUT] = mms_fdoa_interp_ephem( __ )
%     Also return the velocity vectors at interpolated time points.
%
% Parameters
%   T               in, required, type = 1xN int64 (cdf_time_tt2000)
%   R               in, required, type = 3xN float
%   T_OUT           in, required, type = 1xM int64 (cdf_time_tt2000)
%   V               in, required, type = 1xN float
%
% Returns
%   R_OUT           out, required, type = 3xM float
%   V_OUT           out, optional, type = 3xM float
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-08-21      Written by Matthew Argall
%
function [r_out, v_out] = mms_fdoa_interp_ephem(t, r, t_out, v)
	if nargin < 4
		v = zeros(3, 0);
	end
	
	% Convert from tt2000 to seconds for interpolation
	t0    = min( [t(1) t_out(1)] );
	t     = double(t     - t0) * 1e-9;
	t_out = double(t_out - t0) * 1e-9;
	
	% Allocate memory
	r_out = zeros(3, length(t_out));
	v_out = zeros(3, length(t_out));

	% Interpolate
	[r_out(1,:), v_out(1,:)] = mrhermite( t, r(1,:), t_out, v(1,:), 'extrap' );
	[r_out(2,:), v_out(2,:)] = mrhermite( t, r(2,:), t_out, v(2,:), 'extrap' );
	[r_out(3,:), v_out(3,:)] = mrhermite( t, r(3,:), t_out, v(3,:), 'extrap' );
end