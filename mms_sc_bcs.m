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
%   [B_BCS, T] = mms_sc_bcs(SC, INSTR, MODE, TSTART, TEND)
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
%   SC              in, required, type = char
%   INSTR           in, required, type = char
%   MODE            in, required, type = char
%   TSTART          in, optional, type = char
%   TEND            in, optional, type = char
%   'CalDir'        in, optional, type = char default = 'DataDir'
%                   Directory in which to find calibration data.
%   'DataDir'       in, optional, type = char default = pwd()
%                   Directory in which to find search coil data.
%   'Duration'      in, optional, type = char default = 600.0
%                   Time duration over which to calibrate data.
%
% Returns
%   B_BCS           out, required, type=3xN double
%   T               out, optional, type=1xN int64 (cdf_time_tt2000)
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-14      Written by Matthew Argall
%
function [b_bcs, t] = mms_sc_bcs(sc, instr, mode, tstart, tend, varargin)

	% Defaults
	cal_dir   = '';
	data_dir  = pwd();
	duration  = 600.0;

	nOptArgs = length(varargin);
	for ii = 1 : 2 : nOptArgs
		switch varargin{ii}
			case 'CalDir'
				cal_dir   = varargin{ii+1};
			case 'DataDir'
				data_dir  = varargin{ii+1};
			case 'Duration'
				duration  = varargin{ii+1};
			otherwise
				error( ['Parameter name not recognized: "' varargin{ii} '".'] );
		end
	end

%------------------------------------%
% Find and Read B SCM                %
%------------------------------------%
	% Create a file name
	fpattern = mms_construct_filename(sc, instr, mode, 'l1a_sc128', ...
	                                  'Directory', data_dir, ...
	                                  'Tokens',    true);
	
	% Find the files
	[files, nFiles] = MrFile_Search(fpattern,              ...
	                                'TStart',    tstart,   ...
	                                'TEnd',      tend,     ...
	                                'TimeOrder', '%Y%M%d', ...
	                                'Closest',   true);
	assert(nFiles > 0, ['No SCM files found matching "' fpattern '".'] )

	% Create variable names.
	%   - Magnetic field and sample rate.
	b_name  = mms_construct_varname(sc, instr, 'sc128_123');
	sr_name = mms_construct_varname(sc, instr, 'samplerate_sc128');
	
	% Read the data
	[b_123, t] = MrCDF_nRead(files, b_name,  'sTime', tstart, 'eTime', tend);
	[sr, t_sr] = MrCDF_nRead(files, sr_name, 'sTime', tstart, 'eTime', tend);

	% Make row vectors
	b_123 = b_123';
	sr    = sr';
	t     = t';
	t_sr  = t_sr';

%------------------------------------%
% Read the Transfer Function         %
%------------------------------------%
	
	% Find the calibration file
	sc_cal_ftest = fullfile(cal_dir, [sc '_' instr sc(4) '_caltab_*.txt']);
	sc_cal_file  = dir(sc_cal_ftest);
	sc_cal_file  = fullfile(cal_dir, sc_cal_file.name);

	% Read the cal file
	[transfr_fn, f] = mms_sc_read_caltab(sc_cal_file);

%------------------------------------%
% Calibrate the Data                 %
%------------------------------------%
	% Number of points to calibrate.
	nPts = length(t);
	
	% SR is uint16. MATLAB does not like to multiply integers
	% by arrays of doubles. Convert to single.
	sr = single(sr);
	
	% Find the locations of T within T_SR.
	%   - Must convert to double first
	t_sse  = MrCDF_epoch2sse(t, t_sr(1));
	sr_sse = MrCDF_epoch2sse(t_sr);

	% Find points in T before T_SR(1) and after T_SR(end);
	%   - T_SR are given at the beginning of every packet.
	%   - Make sure to capture all data in the last packet.
	deltat  = median( diff(t_sse) );
	iBefore = find( t >= t_sr(1),            1, 'first' ) - 1;
	iAfter  = find( t <= t_sr(end) + deltat, 1, 'last'  ) + 1;
	
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
	b_cal = zeros(3, nPts);
	for ii = 1 : nCal
		b_cal(:, istart(ii):iend(ii)) ...
			= mms_sc_calibrate(b_123(:, istart(ii):iend(ii)), sr(inds(ii)), transfr_fn, f, duration);
	end

%------------------------------------%
% Rotate to BCS                      %
%------------------------------------%

	% BCS is the same as OCS, except for a translation along z,
	% which plays no role in a coordinate transformation. Could
	% otherwise use mms_instr_xxyz2instr('SCM_123', 'BCS')
	sc2bcs = mms_instr_xxyz2ocs('SCM_123');
	
	% Rotate
	b_bcs = mrvector_rotate(sc2bcs, b_cal);
end