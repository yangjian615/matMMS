%
% Name
%   mms_fg_read_cal
%
% Purpose
%   Read MMS fluxgate calibration data.
%
% Calling Sequence
%   CAL_PARAMS = mms_fg_calibrate(FILENAME, TSTART, TSTOP)
%     Read calibration parameters CAL_PARAMS from calibration file(s)
%     FILENAME within the time interval [TSTART, TEND].
%
% Parameters
%   FILENAME        in, required, type=char/cell
%   TSTART          in, required, type=char
%   TSTOP           in, required, type=char
%
% Returns
%   CAL_PARAMS      out, required, type=structure
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
%   2015-03-22      Written by Matthew Argall
%   2015-05-21      Assert that the file name exists.
%
function cal_params = mms_fg_read_cal(filename, tstart, tstop)

	% Make sure the file exist
	assert( exist(filename, 'file') == 2, ['Calibration file does not exist: "' filename '".'] );

	% Dissect the file name
	[sc, instr, mode, level] = mms_dissect_filename(filename);

	% Construct the param name for the cal variables
	if strcmp(mode, 'lorangecal')
		name = ['lo_' level];
	else
		name = ['hi_' level];
	end
	
	% Construct the variable names
	epoch_name  = mms_construct_varname(sc, instr, name, 'Epoch');
	gain_name   = mms_construct_varname(sc, instr, name, 'G');
	dTheta_name = mms_construct_varname(sc, instr, name, 'dTheta');
	dPhi_name   = mms_construct_varname(sc, instr, name, 'dPhi');
	u3_name     = mms_construct_varname(sc, instr, name, 'U3');
	offset_name = mms_construct_varname(sc, instr, name, 'O');
	mpa_name    = mms_construct_varname(sc, instr, name, 'MPA');
	
	% Read data calibration data
	%   - Read all of it because there may not be calibration
	%     parameters within the time range yet.
	%   - Also, dTheta and dPhi do not have DEPEND_0 variables.
	%     However, their values appear to change. I have only
	%     pre-flight calibration files so cannot say.
	[gain_data, t_cal] = MrCDF_nRead(filename, gain_name,   'ColumnMajor', true);
	dTheta_data        = MrCDF_nRead(filename, dTheta_name, 'ColumnMajor', true);
	dPhi_data          = MrCDF_nRead(filename, dPhi_name,   'ColumnMajor', true);
	u3_data            = MrCDF_nRead(filename, u3_name,     'ColumnMajor', true);
	offset_data        = MrCDF_nRead(filename, offset_name, 'ColumnMajor', true);
	mpa_data           = MrCDF_nRead(filename, mpa_name,    'ColumnMajor', true);
	
	% Find the closest calibration time and pick those calibration params.
	%   - Transpose to return row vectors instead of column vectors.
	cal_params = struct( 'Epoch',  t_cal,       ...
	                     'Gain',   gain_data,   ...
	                     'dTheta', dTheta_data, ...
	                     'dPhi',   dPhi_data,   ...
	                     'U3',     u3_data,     ...
	                     'Offset', offset_data, ...
	                     'MPA',    mpa_data);
end