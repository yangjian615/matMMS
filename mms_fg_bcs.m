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
%   Calibrate the fluxgate data and transform into the spacecraft
%   body coordinate system (BCS).
%
% Calling Sequence
%   [B_BCS, T] = mms_fg_calibrate(SC, INSTR, TSTART, TEND)
%     Look in the present working directory for fgm magnetic field
%     and calibration data associated with MMS spacecraft SC,
%     instrument INSTR, and time interval [TSTART, TEND]. Read and
%     calibrate the data, then return the results, B_BCS and T_OUT.
%     TSTART and TEND should be ISO-8601 strings.
%
%   [__] = mms_sc_calibrate(__, 'ParamName', ParamValue)
%     Use any of the parameter name-value pairs below.
%
% Parameters
%   SC              in, required, type = char
%   INSTR           in, required, type = char
%   TSTART          in, required, type = char
%   TEND            in, required, type = char
%   'CalDir'        in, optional, type = char, defualt = 'DataDir';
%                   Directory in which to find calibration files.
%   'DataDir'       in, optional, type = char, defualt = pwd();
%                   Directory in which to find data files.
%
% Returns
%   B_BCS           out, required, type=3xN double
%   T               out, optional, type=1xN int64 (cdf_time_tt2000)
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-13      Written by Matthew Argall
%
function [b_bcs, t, b_smpa, b_omb] = mms_fg_bcs(sc, instr, mode, tstart, tend, varargin)

	% Defaults
	cal_dir   = '';
	data_dir  = pwd();

	nOptArgs = length(varargin);
	for ii = 1 : 2 : nOptArgs
		switch varargin{ii}
			case 'CalDir'
				cal_dir   = varargin{ii+1};
			case 'DataDir'
				data_dir  = varargin{ii+1};
			otherwise
				error( ['Parameter name not recognized: "' varargin{ii} '".'] );
		end
	end
	
	% Assume the calibration files are stored with the data.
	if isempty(cal_dir)
		cal_dir = data_dir;
	end

%------------------------------------%
% Read Cal Data                      %
%------------------------------------%

	% File name patterns
	hical_fname = mms_construct_filename(sc, instr, 'hirangecal', 'l2pre', ...
	                                     'Directory', cal_dir, ...
	                                     'Tokens',    true);
	local_fname = mms_construct_filename(sc, instr, 'lorangecal', 'l2pre', ...
	                                     'Directory', cal_dir, ...
	                                     'Tokens',    true);

	% Find the calibration file
	[hiCal_file, nHiCal] = MrFile_Search(hical_fname);
	[loCal_file, nLoCal] = MrFile_Search(local_fname);

	% Make sure the files exist
	assert( nHiCal > 0, 'HiCal file not found.' );
	assert( nLoCal > 0, 'LoCal file not found.' );

	% Calibrate hi-range
	hiCal = mms_fg_read_cal(hiCal_file, tstart, tend);
	loCal = mms_fg_read_cal(loCal_file, tstart, tend);

%------------------------------------%
% Read Mag Data                      %
%------------------------------------%

	% Create the file names
	fpattern = mms_construct_filename(sc, instr, mode, 'l1a', ...
	                                  'Directory', data_dir, ...
	                                  'Tokens',    true);
	
	% Find the files
	[files, nFiles] = MrFile_Search(fpattern,              ...
	                                'TStart',    tstart,   ...
	                                'TEnd',      tend,     ...
	                                'TimeOrder', '%Y%M%d', ...
	                                'Closest',   true);
	if nFiles == 0
		error( ['No files found matching "' fpattern '".'] );
	end

	% Create variable names
	b_name     = mms_construct_varname(sc, instr, '123');
	if strcmp(instr, 'afg')
		range_name = mms_construct_varname(sc, instr, 'hirange');
	else
		range_name = mms_construct_varname(sc, instr, 'range');
	end

	% Read the magnetometer data
	[b_123,  t]      = MrCDF_Read(files, b_name,     'sTime', tstart, 'eTime', tend);
	[range, t_range] = MrCDF_Read(files, range_name, 'sTime', tstart, 'eTime', tend);

	% Transpose the data to be row vectors.
	t       = t';
	b_123   = b_123';
	range   = range';
	t_range = t_range';

%------------------------------------%
% Calibrate Mag Data                 %
%------------------------------------%
	% Calibrate
	[b_omb, mpa] = mms_fg_calibrate(b_123, t, range, t_range, hiCal, loCal);

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