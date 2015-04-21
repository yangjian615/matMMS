%
% Name
%   mms_fg_calibrate
%
% Purpose
%   Calibrate MMS fluxgate data. This transforms magnetic field from
%   the uncalibrated sensor frame (123) into the orthogonalized
%   magnetometer from (OMB).
%
%   Process:
%     1. Separate hi- and lo-range data
%     2. Calibrate hi- and lo-range data separately
%       a) Map data to nearest calibration time
%       b) Apply offsets
%       c) Create orthogonalization matrix
%       d) Apply orthogonalization
%
% Calling Sequence
%   B_OMB = mms_fg_calibrate(B_123, T, RANGE, T_RANGE, HICAL, LOCAL)
%     Calibrate and orthogonalize magnetometer data in the non-orthogonal
%     123 coordinate system B_123. B_123 has cdf_time_tt2000 time tags T
%     and is in range RANGE at times T_RANGE. The calibration parameters
%     for hi- and lo-range data are HICAL and LOCAL, respectively, and can
%     be obtained from mms_fg_read_cal.m. Return results in B_OMB.
%
%   [B_OMB, MPA] = mms_fg_calibrate(B_123, T, RANGE, T_RANGE, HICAL, LOCAL)
%     Also return the z-MPA axis as viewed from BCS.
%
% Parameters
%   B_123           in, required, type = 3xN double
%   T               in, required, type = 1xN int64 (cdf_time_tt2000)
%   RANGE           in, required, type = 1xN logical
%   T_RANGE         in, required, type = 1xN int64 (cdf_time_tt2000)
%   HICAL           in, required, type = struct
%   LOCAL           in, required, type = struct
%
% Returns
%   B_OMB           out, required, type=3xN double
%   MPA             out, optional, type=3xN double
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-13      Written by Matthew Argall
%
function [B_omb, mpa] = mms_fg_calibrate(B_123, t, range, t_range, hiCal, loCal)

	%
	% Histogram ranges
	%   - T_RANGE are reported once per packet
	%   - B_FG times occur multiple times per packet. 
	%   - T_RANGE will serve as buckets into which
	%     magnetic field times will fall -- essentially
	%     rounding T_FG down to the nearest T_RANGE
	%     and using its range flag.
	%   - histc() reqires T_FG to fall completely within
	%     the range of T_RANGE or invalid indices will
	%     be returned. Extend T_RANGE if need be.
	%
	% However, the packet time marks the beginning of the
	% data capture, so there should be no problem.
	%
	
	% histc() does not take int64.
	t_sse       = MrCDF_epoch2sse(t, t_range(1));
	t_range_sse = MrCDF_epoch2sse(t_range, t_range(1));
	
	% Map each T_FG onto the values of RANGE.
	[~, range_inds] = histc(t_sse, t_range_sse);
	
	%
	% Turns out, it is possible for data points to not
	% have a range flag. At the beginning of the array,
	% the first point may not. At the end of the array,
	% all points in the last packet will not because
	% histc() does not have a right-edge for the last
	% bin.
	%
	if t(1) < t_range(1)
		% Find the first point with range information.
		iInRange = find( t >= t_range(1), 1, 'first');
		wrn_msg  = sprintf('First %d points lack range info.', iInRange-1);
		warning('mms_fg_calibrate:range', wrn_msg);
		
		% Set out-of-range values to the first in-range value.
		range_inds(1:iInRange) = range_inds(iInRange);
	end
	if t(end) > t_range(end)
		% Find the last point with range information.
		iInRange = find( t <= t_range(end), 1, 'last');
		wrn_msg  = sprintf('Last %d points may lack range info.', length(t) - iInRange);
		warning('mms_fg_calibrate:range', wrn_msg);
		
		% Set out-of-range values to the last in-range value.
		range_inds(iInRange+1:end) = range_inds(iInRange);
	end
	
	% hirange data -- logical array
	%   1 if hirange
	%   0 if lorange
	hirange = range(range_inds) == 1;
	lorange = ~hirange;
	nHi     = sum( hirange );
	nLo     = sum( lorange );

%------------------------------------%
% Apply Calibration                  %
%------------------------------------%
	% Allocate memory to output
	B_omb = B_123;
	
	% hirange data
	if nHi > 0
		[B_omb(:, hirange), zmpaHi] = mms_fg_calibrate_apply( B_123(:, hirange), t(hirange), hiCal );
	end
	
	% lorange data
	if nLo > 0
		[B_omb(:, lorange), zmpaLo] = mms_fg_calibrate_apply( B_123(:, lorange), t(lorange), loCal );
	end

%------------------------------------%
% MPA z-Axis                         %
%------------------------------------%
	if nargout > 1
		mpa = zeros(size(B_omb));
		if nHi > 0
			mpa(:, hirange ) = zmpaHi;
		end
		if nLo > 0
			mpa(:, lorange ) = zmpaLo;
		end
	end
end


%
% Name
%   mms_fg_calibrate_apply
%
% Purpose
%   A helper routine for the mms_fg_calibrate function. Will calibrate
%   either hi-range or lo-range data.
%
%   Process:
%     1. Map data to nearest calibration time
%     2. Apply offsets
%     3. Create orthogonalization matrix
%     4. Apply orthogonalization
%
% Calling Sequence
%   B_CAL = mms_fg_calibrate_apply(B, TIME, CAL_PARAMS)
%     Calibrate MMS fgm data B with time stamps TIME using the
%     calibration parameters stored in CAL_PARAMS. Return the
%     calibrated magnetic field B_CAL. CAL_PARAMS can be obtained
%     from mms_fg_read_cal.m.
%
%   [B_CAL, CAL] = mms_fg_calibrate_apply(__)
%     Also return the interpolated calibration matrix.
%
% Parameters
%   B               in, required, type=3xN double
%   TIME            in, required, type=int64 (cdf_time_tt2000)
%   CAL_PARAMS      in, required, type=struct
%
% Returns
%   B_CAL           out, required, type=3xN double
%   CAL             out, optional, type=struct
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-13      Written by Matthew Argall
%   2015-04-21      Avoid interpolating calibration parameters. - MRA
%
function [b_cal, zmpa] = mms_fg_calibrate_apply(B_123, time, cal_params)

%------------------------------------%
% Map Cal Params to Values of B      %
%------------------------------------%
	
	%
	% Interpolate the calibration parameters
	%   - Always take the last calibrated point (no interpolation).
	%   - Ken: cal values will always be constant for L1B,
	%          regardless of how new the cal file is.
	%
	% mms_fg_calparams2matrix is super slow because it has to find the
	% inverse of all of the orthogonalization matrices. Since all we are
	% doing is mapping each point in B_123 onto the last calibration time,
	% there is no need to obtain calibration parameters for every value
	% of B_123. We need only to compute the few matrices within our time
	% of interest and apply the same matrix repeatedly to the appropriate
	% values oB_123
	%
	
	%
	% Do not find a calibration parameter at each point
	%
	%     cal = mms_fg_interp_cal(cal_params, time, true);
	%
	% Instead, map each value of B_123 onto the appropriate set of
	% calibration parameters.
	%
	
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
	
	% MrValue_Locate() cannot handle int64.
	time_sse  = MrCDF_epoch2sse(time,      cal_epoch(1));
	t_cal_sse = MrCDF_epoch2sse(cal_epoch, cal_epoch(1));
	
	% Map times to closest calibration time.
	cal_inds = MrValue_Locate(t_cal_sse, time_sse);

%------------------------------------%
% Apply Calibration Parameters       %
%------------------------------------%
	% Subtract DC offsets
	b_cal      = B_123;
	b_cal(1,:) = B_123(1,:) - cal_params.('Offset')(1,cal_inds);
	b_cal(2,:) = B_123(2,:) - cal_params.('Offset')(2,cal_inds);
	b_cal(3,:) = B_123(3,:) - cal_params.('Offset')(3,cal_inds);
	
	% Create orthogonalization matrix
	orthog_mat = mms_fg_calparams2matrix( cal_params.('Gain'), cal_params.('dTheta'), ...
	                                      cal_params.('dPhi'), cal_params.('U3') );
	
	% Orthogonalize the data
	b_cal = mrvector_rotate( orthog_mat(:,:,cal_inds), b_cal );

	% Major principle axis
	zmpa = cal_params.('MPA')(:,cal_inds);
end