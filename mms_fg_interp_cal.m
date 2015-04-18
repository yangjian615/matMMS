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
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-12      Written by Matthew Argall
%   2015-04-17      Finding exprapoled indices incorrectly. Now use
%                     MrValue_Locate to do so. - MRA.
%
function cal = mms_fg_interp_cal(cal_params, time, lastval)

	if nargin < 3
		tf_lastval = false;
	end

	% Number of points in the output
	npts = length(time);
	
	% Allocate memory
	cal = struct( 'Epoch',  zeros(1, npts), ...
	              'Gain',   zeros(3, npts), ...
	              'dPhi',   zeros(2, npts), ...
	              'dTheta', zeros(2, npts), ...
	              'U3',     zeros(2, npts), ...
	              'Offset', zeros(3, npts), ...
	              'MPA',    zeros(3, npts) );
	
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
	if sum(tf_extrap) > 0 || tf_lastval
		% Use the last value
		for ii = 1 : 3
			cal.('Gain')(ii, :)   = cal_params.('Gain')(ii, ilast_val);
			cal.('Offset')(ii, :) = cal_params.('Offset')(ii, ilast_val);
			cal.('MPA')(ii, :)    = cal_params.('MPA')(ii, ilast_val);
		end
		for ii = 1 : 2
			cal.('dPhi')(ii, :)   = cal_params.('dPhi')(ii, ilast_val);
			cal.('dTheta')(ii, :) = cal_params.('dTheta')(ii, ilast_val);
			cal.('U3')(ii, :)     = cal_params.('U3')(ii, ilast_val);
		end
	end
	
	% Interpolation
	if sum(tf_interp) > 0 && ~tf_lastval
		t_interp = times(tf_interp);

		% Use the last value
		for ii = 1 : 3
			cal.('Gain')(ii, tf_interp)   = interp1(cal_epoch, cal_params.('Gain')(ii, :), t_inerp);
			cal.('Offset')(ii, tf_interp) = interp1(cal_epoch, cal_params.('Offset')(ii, :), t_inerp);
			cal.('MPA')(ii, tf_interp)    = interp1(cal_epoch, cal_params.('MPA')(ii, :), t_inerp);
			% Slerp the MPA... what is slerp?
		end
		for ii = 1 : 2
			cal.('dPhi')(ii, tf_interp)   = interp1(cal_epoch, cal_params.('dPhi')(ii, :), t_inerp);
			cal.('dTheta')(ii, tf_intepr) = interp1(cal_epoch, cal_params.('dTheta')(ii, :), t_inerp);
			cal.('U3')(ii, tf_interp)     = interp1(cal_epoch, cal_params.('U3')(ii, :), t_inerp);
		end
	end
end