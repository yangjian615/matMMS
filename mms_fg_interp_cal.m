%
% Name
%   mms_fg_interp_cal
%
% Purpose
%   Calibrate MMS fluxgate data.
%
% Calling Sequence
%   CAL = mms_fg_interp_cal(CAL_PARAMS, TIME)
%     Interpolate fgm calibration prameters CAL_PARAMS returned
%     from a call to mms_fg_read_cal.m to the times TIME, and
%     return an updated structure, CAL.
%
%   CAL = mms_fg_interp_cal(CAL_PARAMS, TIME, LASTVAL)
%     Instead of interpolating, use the last calibration parameter.
%
% Parameters
%   CAL_PARAMS      in, required, type = 3xN double
%   TIME            in, required, type = 1xN int64 (cdf_time_tt2000)
%   LASTVAL         in, required, type = boolean, default = false
%
% Returns
%   CAL             out, required, type=struct
%                   Fields are:
%                     'Epoch'  -  CDF_TIME_TT2000 time tags.
%                     'Gain'   -  Instrument gain for each axis.
%                     'dTheta' -  Orthogonalization angle.
%                     'dPhi'   -  Orthogonalization angle.
%                     'U3'     -  Orthogonalization angle.
%                     'Offset' -  DC offset for each axis.
%                     'MPA'    -  Major principal axis of inertia.
%
% See Also
%   mms_fg_calibrate.m
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-12      Written by Matthew Argall
%   2015-04-17      Finding exprapoled indices incorrectly. Now use
%                     MrValue_Locate to do so. - MRA.
%   2015-07-27      Update to support new version of calibration files. - MRA
%
function cal = mms_fg_interp_cal(cal_params, time, lastval)

	if nargin < 3
		tf_lastval = false;
	end

	% Number of points in the output
	npts = length(time);
	
	% Allocate memory
	cal = struct( 'Epoch',   zeros(1, npts), ...
	              'Gain',    zeros(3, npts), ...
	              'dPhi',    zeros(2, npts), ...
	              'dTheta',  zeros(2, npts), ...
	              'U3',      zeros(2, npts), ...
	              'Offset',  zeros(3, npts), ...
	              'MPA',     zeros(3, npts), ...
	              'gPrime',  zeros(3, npts), ...
	              'stemp',   zeros(1, npts), ...
	              'etemp',   zeros(1, npts), ...
	              'stemp_r', zeros(1, npts), ...
	              'etemp_r', zeros(1, npts), ...
	              'alpha_s', zeros(3, npts), ...
	              'alpha_e', zeros(3, npts) ...
	            );
	
	% Extract the epoch times for easy access
	cal_epoch = cal_params.('Epoch');
	
	% Extrapolate or interpolate?
	tf_before = time <  cal_epoch(1);
	tf_interp = time >= cal_epoch(1)   & time <= cal_epoch(end);
	tf_extrap = time >  cal_epoch(end);
	
	% Should not have values before the first calibration time
	if sum(tf_before) > 0
		error( 'Times must be after the first calibration period.' );
	end
	
	% histc() does not like int64.
	time_sse  = MrCDF_epoch2sse(time,                 cal_params.('Epoch')(1));
	t_cal_sse = MrCDF_epoch2sse(cal_params.('Epoch'), cal_params.('Epoch')(1));
	
	% Map times to closest calibration time.
	ilast_val = MrValue_Locate(t_cal_sse, time_sse);

	% Extrapolation
	%   - ILAST_VAL is the location of each output point within the cal data
	%   - If we are also interpolating, those values will be overwritten.
	if sum(tf_extrap) > 0 || tf_lastval
		% Use the last value
		for ii = 1 : 3
			% 3-vectors
			cal.('Gain')(ii, :)   = cal_params.('Gain')(ii, ilast_val);
			cal.('Offset')(ii, :) = cal_params.('Offset')(ii, ilast_val);
			cal.('MPA')(ii, :)    = cal_params.('MPA')(ii, ilast_val);
			cal.('gPrime')(ii, :) = cal_params.('gPrime')(ii, ilast_val);
			cal.('alpha_s')(ii, :) = cal_params.('alpha_s')(ii, ilast_val);
			cal.('alpha_e')(ii, :) = cal_params.('alpha_e')(ii, ilast_val);

			% 2-vectors
			if ii <= 2
				cal.('dPhi')(ii, :)   = cal_params.('dPhi')(ii, ilast_val);
				cal.('dTheta')(ii, :) = cal_params.('dTheta')(ii, ilast_val);
				cal.('U3')(ii, :)     = cal_params.('U3')(ii, ilast_val);
			end
		end
			
		% Scalars
		cal.('stemp')   = cal_params.('stemp')(ilast_val);
		cal.('etemp')   = cal_params.('etemp')(ilast_val);
		cal.('stemp_r') = cal_params.('stemp_r')(ilast_val);
		cal.('etemp_r') = cal_params.('etemp_r')(ilast_val);
	end
	
	% Interpolation
	if sum(tf_interp) > 0 && ~tf_lastval
		t_interp = times(tf_interp);

		% Use the last value
		for ii = 1 : 3
			% 3-Vectors
			cal.('Gain')(ii, tf_interp)   = interp1(cal_epoch, cal_params.('Gain')(ii, :), t_inerp);
			cal.('Offset')(ii, tf_interp) = interp1(cal_epoch, cal_params.('Offset')(ii, :), t_inerp);
			cal.('MPA')(ii, tf_interp)    = interp1(cal_epoch, cal_params.('MPA')(ii, :), t_inerp);
			% Slerp the MPA... what is slerp?
			cal.('gPrime')(ii, tf_interp)  = interp1(cal_epoch, cal_params.('gPrime')(ii, :), t_inerp);
			cal.('alpha_s')(ii, tf_interp) = interp1(cal_epoch, cal_params.('alpha_s')(ii, :), t_inerp);
			cal.('alpha_e')(ii, tf_interp) = interp1(cal_epoch, cal_params.('alpha_e')(ii, :), t_inerp);

			% 2-Vectors
			if ii <= 2 ii = 1 : 2
				cal.('dPhi')(ii, tf_interp)   = interp1(cal_epoch, cal_params.('dPhi')(ii, :), t_inerp);
				cal.('dTheta')(ii, tf_intepr) = interp1(cal_epoch, cal_params.('dTheta')(ii, :), t_inerp);
				cal.('U3')(ii, tf_interp)     = interp1(cal_epoch, cal_params.('U3')(ii, :), t_inerp);
			end
		end
		
		% Scalars
		cal.('stemp')(ii, tf_interp)   = interp1(cal_epoch, cal_params.('stemp')(ii, :), t_inerp);
		cal.('etemp')(ii, tf_interp)   = interp1(cal_epoch, cal_params.('etemp')(ii, :), t_inerp);
		cal.('stemp_r')(ii, tf_interp) = interp1(cal_epoch, cal_params.('stemp_r')(ii, :), t_inerp);
		cal.('etemp_r')(ii, tf_interp) = interp1(cal_epoch, cal_params.('etemp_r')(ii, :), t_inerp);
	end
end