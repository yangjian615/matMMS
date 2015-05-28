%
% Name
%   mms_sc_calibrate
%
% Purpose
%   Calibrate MMS search coil data and transform into BCS.
%
%   Process:
%     1. Find and read calibration file
%     2. Find and read search coil data
%     3. Create compensation array for transfer function
%     4. Convert numbers to nano-Tesla
%     5. Calibrate
%       a) Fourier transform intervals
%       b) Apply compensate array
%       c) Inverse fourier transform
%
% Calling Sequence
%   [T, B_BCS] = mms_sc_bcs(SC, INSTR, MODE, TSTART, TEND)
%     Find and read search coil magnetometer and transfer function
%     data from MMS scacecraft SC (e.g. 'mms1') from instrument
%     INSTR (i.e. 'scm') in telemetry mode MODE (e.g. 'comm') during
%     the time interval [TSTART, TEND]. Calibrate the data and
%     transform into the spacecraft body frame (BCS). Return the data
%     and its time tags B_BCS and T.
%
%   [__] = mms_sc_bcs(__, 'ParamName', ParamValue)
%     Any parameter name-value pair found below.
%
% Parameters
%   SC_FILES        in, required, type = char/cell
%   CAL_FILE        in, required, type = char
%   TSTART          in, optional, type = char
%   TEND            in, optional, type = char
%   DURATION        in, optional, type = char default = 600.0
%
% Returns
%   T               out, optional, type=1xN int64 (cdf_time_tt2000)
%   B_BCS           out, required, type=3xN double
%   B_OMB           out, required, type=3xN double
%   B_123           out, required, type=3xN double
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-14      Written by Matthew Argall
%   2015-05-11      Take file names as inputs. Return data in OMB and 123. - MRA
%
function [t, b_bcs, b_omb, b_123] = mms_sc_bcs(sc_files, cal_file, tstart, tend, duration)

	if nargin < 5
		duration = 600.0;
	end

%------------------------------------%
% Read Data                          %
%------------------------------------%

	% Read the L1A files
	sc_l1a = mms_sc_read_l1a(sc_files, tstart, tend);

	% Read the cal file
	[transfr_fn, f] = mms_sc_read_caltab(cal_file);

%------------------------------------%
% Calibrate the Data                 %
%------------------------------------%
	% Number of points to calibrate.
	nPts = length(sc_l1a.tt2000);
	
	% SR is uint16. MATLAB does not like to multiply integers
	% by arrays of doubles. Convert to single.
	sr = single(sc_l1a.sample_rate);
	
	% Find the locations of TT2000 within TT2000_PACKET.
	%   - Must convert to double first
	t_sse  = MrCDF_epoch2sse(sc_l1a.tt2000, sc_l1a.tt2000_packet(1));
	sr_sse = MrCDF_epoch2sse(sc_l1a.tt2000_packet);

	% Find points in T before T_SR(1) and after T_SR(end);
	%   - T_SR are given at the beginning of every packet.
	%   - Make sure to capture all data in the last packet.
	deltat  = median( diff(t_sse) );
	iBefore = find( sc_l1a.tt2000 >= sc_l1a.tt2000_packet(1),            1, 'first' ) - 1;
	iAfter  = find( sc_l1a.tt2000 <= sc_l1a.tt2000_packet(end) + deltat, 1, 'last'  ) + 1;
	
	% Number of points before the first packet.
	if isempty(iBefore)
		nBefore = 0;
	else
		nBefore = iBefore;
	end
	
	% Number of points after the last packet
	if isempty(iAfter)
		nAfter = 0;
	else
		nAfter = nPts - iAfter + 1;
	end

	% Data starts before sample rate.
	%   - Issue a warning if the sample rate is different
	%     from the first reported sample rate.
	%   - If nBefore == 1, assume it is with the same
	%     sampling rate as the first packet.
	if nBefore > 1
		% Number of nominal samples between points
		nsr = round( diff( t_sse(1:iBefore) ) * sr(1) );
		if sum(nsr) ~= nBefore - 1
			msg = sprintf( ['%d points out of cal interval. ' ...
			                'Mean sample rate %f * expected.'], ...
			               iBefore, mean(nsr) )
			warning('MMS_SC_BSC:Calibrate', msg)
		end
	end
	
	% Data ends after sample rate.
	%   - Issue a warning if the sample rate is different
	%     from the last reported sample rate.
	%   - If nAfter == 1, assume it is with the same
	%     sampling rate as the last packet.
	if nAfter > 1
		% Number of nominal samples between points
		nsr = round( diff( t_sse(iAfter:end) ) * sr(end) );
		if sum(nsr) ~= nAfter - 1
			msg = sprintf( ['%d points out of cal interval. ' ...
			                'Mean sample rate %f * expected.'], ...
			               nAfter, mean(nsr) )
			warning('MMS_SC_BSC:Calibrate', msg)
		end
	end

	% Find locations of T_SSE in SR_SSE
	inds   = MrValue_Locate(sr_sse, t_sse, 'RoundUp', true);

	% Calibrate each interval with a different sample rate
	dsr    = diff(sr);
	istart = [ 1 find( dsr ~= 0 ) ];
	iend   = [ istart(2:end)-1 nPts ];
	nCal   = length(istart);

	% Calibrate
	b_omb = zeros(3, nPts);
	for ii = 1 : nCal
		b_omb(:, istart(ii):iend(ii)) ...
			= mms_sc_calibrate(sc_l1a.b_123(:, istart(ii):iend(ii)), sr(inds(ii)), transfr_fn, f, duration);
	end
	
	% Extract data
	if nargout > 3
		b_123 = sc_l1a.b_123;
	end
	t = sc_l1a.tt2000;
	
	% Clear the data structure
	clear sc_l1a

%------------------------------------%
% Rotate to BCS                      %
%------------------------------------%

	% BCS is the same as OCS, except for a translation along z,
	% which plays no role in a coordinate transformation. Could
	% otherwise use mms_instr_xxyz2instr('SCM_123', 'BCS')
	sc2bcs = mms_instr_xxyz2ocs('SCM_123');
	
	% Rotate
	b_bcs = mrvector_rotate(sc2bcs, b_omb);
end