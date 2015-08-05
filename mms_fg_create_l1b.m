%
% Name
%   mms_fg_bcs
%
%   Process
%     1. Find and read calibration files
%     2. Find and read data files
%     3. Calibrate data
%     4. Transform from OMB to SMPA
%.    5. Transform from SMPA to BCS (equivalent to OCS).
%
% Purpose
%   Read digital and analog fluxgate (DFG & AFG) magnetometer level
%   1A data and turn it into level 1B data. L1B implies calibrated
%   data in the spinning spacecraft body coordinate system.
%
% Calling Sequence
%   [T, B_BCS] = mms_fg_create_l1b(FILES, HICAL_FILE, LOCAL_FILE)
%     Read, calibrate and rotate into BCS level 1A fluxgate data
%     from files with names FILES, using hi- and lo-range
%     calibration data found in files HICAL_FILE and LOCAL_FILE.
%
%   [T, B_BCS] = mms_fg_create_l1b(..., TSTART, TEND)
%     Return data within the time interval beginning at TSTART and
%     ending at TEND.
%
%   [__, B_SMPA, B_OMB, B_123] = mms_fg_create_l1b(__)
%     Also return data in the SMPA and 123 coordinate systems.
%
% Parameters
%   FILES           in, required, type = char/cell
%   HICAL_FILE      in, required, type = char
%   LOCAL_FILE      in, required, type = char
%   TSTART          in, required, type = char
%   TEND            in, required, type = char
%
% Returns
%   T               out, optional, type=1xN int64 (cdf_time_tt2000)
%   B_BCS           out, required, type=3xN double
%   B_SMPA          out, required, type=3xN double
%   B_OMB           out, required, type=3xN double
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-13      Written by Matthew Argall
%   2015-05-21      Take filenames as input, not pieces of filenames. - MRA
%   2015-06-21      Renamed from mms_fg_bcs to mms_fg_create_l1b. - MRA
%
function [t, b_bcs, b_smpa, b_omb, b_123] = mms_fg_create_l1b(files, hiCal_file, loCal_file, tstart, tend, hk_file)

	% Defaults
	if nargin() < 4
		tstart = '';
	end
	if nargin() < 5
		tend = '';
	end
	if nargin() < 6
		hk_file = '';
	end
	
	%
	% Right now 2015-07-29, Ken is using the reference temperature to calibrate
	% DFG/AFG data. In-fight temperature measurements from the 0x10e housekeeping
	% files are ignored.
	%

%------------------------------------%
% Read Cal Data                      %
%------------------------------------%
	% Check calibration file names
	assert( ischar(hiCal_file) && isrow(hiCal_file), 'HiCal_File must be a single file name.' );
	assert( ischar(loCal_file) && isrow(loCal_file), 'LoCal_File must be a single file name.' );

	% Calibration data
	[hiCal, hiConst] = mms_fg_read_cal(hiCal_file, tstart, tend);
	[loCal, loConst] = mms_fg_read_cal(loCal_file, tstart, tend);
	
	% Sensor temperature data
	if isempty(hk_file)
		hk_0x10e = [];
	else
		hk_0x10e = mms_hk_read_0x10e(hk_file, tstart, tend);
	end

%------------------------------------%
% Read Mag Data                      %
%------------------------------------%

	% Read files
	fg_l1a = mms_fg_read_l1a(files, tstart, tend);

%------------------------------------%
% Calibrate Mag Data                 %
%------------------------------------%
	% Calibrate
	[b_omb, mpa] = mms_fg_calibrate( fg_l1a.b_123, fg_l1a.tt2000,    ...
	                                 fg_l1a.range, fg_l1a.tt2000_ts, ...
	                                 hiCal, loCal, hiConst, loConst, ...
	                                 hk_0x10e );

	% Extract data
	if nargout() > 5
		b_123 = fg_l1a.b_123;
	end
	t = fg_l1a.tt2000;
	
	% Clear the data structure
	clear fg_l1a

%------------------------------------%
% Transform from OMB to SMPA         %
%------------------------------------%

	% OMB -> SMPA
	omb2smpa = mms_fg_xomb2smpa();
	b_smpa   = mrvector_rotate(omb2smpa, b_omb);

%------------------------------------%
% Transform from SMPA to BCS         %
%------------------------------------%

	% Build transformation from SMPA to BCS
	bcs2smpa = mms_fg_xbcs2smpa(mpa);
	smpa2bcs = permute(bcs2smpa, [2 1 3]);

	% Rotate to SMPA
	b_bcs = mrvector_rotate(smpa2bcs, b_smpa);
end