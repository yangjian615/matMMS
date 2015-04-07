%
% Name
%   mms_fg_calibrate
%
% Purpose
%   Calibrate fluxgate magnetometer data from the MMS mission.
%     1) Read B, T, and range from mag file
%     2) Find appropriate calibration file (hi/lo range)
%     3) Split B into hi- and low-range elements
%       a. Subtract offset
%       b. Create othogonalization matrix from calibration parameters.
%       c. Orthogonalize B
%
% Calling Sequence
%   [B_DATA] = mms_fg_calibrate(FILENAME)
%     Read magnetometer data from the file named FILENAME, calibrate it
%     by applying DC offsets, gain, and orthogonalization matrices, and
%     return the result as B_DATA.
%
% Parameters
%   FILENAME        in, required, type=char
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-03-22      Written by Matthew Argall
%
function [b_fg, t_fg] = mms_fg_calibrate(filename, cal_dir)

	% Create file names
	assert(exist(filename, 'file') == 2, ...
		     ['Flux gate data file does not exist: "' filename '".']);
	assert(exist(cal_dir, 'dir') == 7, ...
		     ['Calibration file directory does not exist: "' cal_dir '".']);

	% Dissect the FG data file name
	[sc, instr] = mms_dissect_filename(filename);
	
%------------------------------------%
% Read Magnetometer Data             %
%------------------------------------%

	% Create variable names
	b_name     = mms_construct_varname(sc, instr, '123');
	range_name = mms_construct_varname(sc, instr, 'range');
	
	
	% Read the magnetometer data
	fg_data = spdfcdfread(filename, ...
		                    'Variables',     {'Epoch' b_name range_name}, ...
                        'CombineRecords', true, ...
												'KeepEpochAsIs',  true);
								
	% Extract the data
	t_fg       = fg_data{1}';
	b_fg       = fg_data{2}';
	range_data = fg_data{3}';
	clear fg_data
	
	% Hi- and lo-range data
	iHi = find(range_data);
	iLo = find(~range_data);

%------------------------------------%
% Calibrate Hi-Range Data            %
%------------------------------------%
	if ~isempty(iHi)
		% Find the calibration file
		hiCal_file = mms_file_search(sc, instr, 'hirangecal', 'l2pre', ...
			                           'Directory', cal_dir);
		
		% Can handle only one file
		assert( length(hiCal_file) == 1, 'More than one cal file found. Cannot proceed.' );
		
		% Calibrate hi-range
		b_fg(:, iHi) = mms_fg_calibrate_hilo(hiCal_file{1}, b_fg(:, iHi), t_fg(1));
	end
	
%------------------------------------%
% Calibrate Lo-Range Data            %
%------------------------------------%
	if ~isempty(iLo)
		% Find the calibration file
		loCal_file = mms_file_search(sc, instr, 'lorangecal', 'l2pre', ...
			                           'Directory', cal_dir);
		
		% Can handle only one file
		assert( length(loCal_file) == 1, 'More than one cal file found. Cannot proceed.' );
		
		% Calibrate hi-range
		b_fg(:, iLo) = mms_fg_calibrate_hilo(loCal_file{1}, b_fg(:, iLo), t_fg(1));
	end
end


%
% Name
%   mms_fg_calibrate_hilo
%
% Purpose
%   A helper routine for the mms_fg_calibrate function. Will calibrate
%   either hi-range or lo-range data.
%
% Calling Sequence
%   [B_DATA] = mms_fg_calibrate(FILENAME)
%     Read magnetometer data from the file named FILENAME, calibrate it
%     by applying DC offsets, gain, and orthogonalization matrices, and
%     return the result as B_DATA.
%
% Parameters
%   FILENAME        in, required, type=char
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-03-22      Written by Matthew Argall
%
function [b_cal] = mms_fg_calibrate_hilo(filename, b_data, t_fg)

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
	theta_name  = mms_construct_varname(sc, instr, name, 'dTheta');
	phi_name    = mms_construct_varname(sc, instr, name, 'dPhi');
	u3_name     = mms_construct_varname(sc, instr, name, 'U3');
	offset_name = mms_construct_varname(sc, instr, name, 'O');
	
	% Read data calibration data 
	cal_data = spdfcdfread(filename, ...
		                     'Variables',     {epoch_name gain_name theta_name phi_name u3_name offset_name}, ...
                         'CombineRecords', true, ...
												 'KeepEpochAsIs',  true);
	
	% Extract the data
	%   - mms_fg_calparams2matrix expects 3xN or 2xN.
	%   - spdfcdfread reaturns Nx3 or Nx2
	%  => transpose data
	t_cal       = cal_data{1}';
	gain_data   = cal_data{2}';
	theta_data  = cal_data{3}';
	phi_data    = cal_data{4}';
	u3_data     = cal_data{5}';
	offset_data = cal_data{6}';
	clear cal_data
	
	% Find the closest calibration time and pick those calibration params.
	iCal       = find(t_cal < t_fg, 1, 'last');
	gain_data  = gain_data(:, iCal);
	theta_data = theta_data(:, iCal);
	phi_data   = phi_data(:, iCal);
	u3_data    = u3_data(:, iCal);

	% Subtract DC offsets
	b_cal      = b_data;
	b_cal(1,:) = b_data(1,:) - offset_data(1);
	b_cal(2,:) = b_data(2,:) - offset_data(2);
	b_cal(3,:) = b_data(3,:) - offset_data(3);
	
	% Create orthogonalization matrix
	[~, orthog_mat] = mms_fg_calparams2matrix(gain_data, theta_data, phi_data, u3_data);
	
	% Orthogonalize the data
	b_cal = mrvector_rotate(orthog_mat, b_cal);
end