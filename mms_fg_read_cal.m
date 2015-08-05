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
%                     'Epoch'   -  CDF_TIME_TT2000 time tags.
%                     'Gain'    -  Instrument gain for each axis.
%                     'dTheta'  -  Orthogonalization angle.
%                     'dPhi'    -  Orthogonalization angle.
%                     'U3'      -  Orthogonalization angle.
%                     'Offset'  -  DC offset for each axis.
%                     'MPA'     -  Major principal axis of inertia.
%                     'gprime'  -  Unitless on-orbit gain correction, relative to temperature-corrected gain
%                     'alpha_s' -  linear parameter for sensor temperature correction, normalized to reference temperature, of form (1+alpha_s(stemp-stemp_r))
%                     'alpha_e' -  linear parameter for electronics temperature correction, normalized to reference temperature, of form (1+alpha_e(etemp-etemp_r))
%                     'stemp'   -  effective sensor temperature during calibration period.
%                     'etemp'   -  effective electronics temperature during calibration period.
%                     'stemp_r' -  reference sensor temperature for on-orbit calibrations.
%                     'etemp_r' -  reference electronics temperature for on-orbit calibrations.
%                     'c'       -  coupling matrix from sensor 123 to orthoganalized magnetometer boom XYZ
%                     'oxyz'    -  offsets in nT, in orthogonalized magnetometer boom (OMB) system.
%                     'ix_12'   -  orthogonality angle between sensor axes 1 and 2
%                     'ix_13'   -  orthogonality angle between sensor axes 1 and 3
%                     'ix_23'   -  orthogonality angle between sensor axes 2 and 3
%
%   CAL_CONST       out, required, type=structure
%                   Parameters determined at Braunschweig. Fields are:
%                     'psf'      - preliminary scale factor used to multiply raw counts to yield nT
%                     'rsf'      - round scale factor used to divide raw counts to yield L1A pseudo-nT.
%                     'h_stemp'  - gain correction (forward gain) at 0 degrees C
%                     'm_stemp'  - slope of gain correction with respect to sensor temperature
%                     'etempp_c' - electronics temperature during calibrations
%                     'm_etemp'  - slope of gain correction with respect to electronics temperature
%                     'nl_a'     - non-linearity correction. linear polynomial coefficient for each sensor. to be applied after temperature correction. x' = nl_a*x + nl_b*x^2 + nl_c*x^3
%                     'nl_b'     - non-linearity correction.  quadratic polynomial coefficient for each sensor.  to be applied after temperature correction. x' = nl_a*x + nl_b*x^2 + nl_c*x^3
%                     'nl_c'     - non-linearity correction.  cubic polynomial coefficient for each sensor: to be applied after temperature correction. x' = nl_a*x + nl_b*x^2 + nl_c*x^3
%                     'dt13'     - relative time delay between sensor axis 1 and sensor axis 3. (axis 3 is the reference for L1A time tag)
%                     'dt23'     - relative time delay between sensor axis 2 and sensor axis 3. (axis 3 is the reference for L1A time tag)
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-03-22      Written by Matthew Argall
%   2015-05-21      Assert that the file name exists.
%   2015-07-27      Update to support newer version of cal files.
%
function [cal_params, cal_const] = mms_fg_read_cal(filename, tstart, tstop)

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
	
%------------------------------------%
% Variable Names                     %
%------------------------------------%
	% Construct the variable names
	epoch_name  = mms_construct_varname(sc, instr, name, 'Epoch'   );
	gain_name   = mms_construct_varname(sc, instr, name, 'G'       );
	dTheta_name = mms_construct_varname(sc, instr, name, 'dTheta'  );
	dPhi_name   = mms_construct_varname(sc, instr, name, 'dPhi'    );
	u3_name     = mms_construct_varname(sc, instr, name, 'U3'      );
	offset_name = mms_construct_varname(sc, instr, name, 'O'       );
	mpa_name    = mms_construct_varname(sc, instr, name, 'MPA'     );
	gprime_name = mms_construct_varname(sc, instr, name, 'gprime'  );
	alphas_name = mms_construct_varname(sc, instr, name, 'alpha_s' );
	alphae_name = mms_construct_varname(sc, instr, name, 'alpha_e' );
	stemp_name  = mms_construct_varname(sc, instr, name, 'stemp'   );
	etemp_name  = mms_construct_varname(sc, instr, name, 'etemp'   );
	stempr_name = mms_construct_varname(sc, instr, name, 'stemp_r' );
	etempr_name = mms_construct_varname(sc, instr, name, 'etemp_r' );
	c_name      = mms_construct_varname(sc, instr, name, 'c'       );
	oxyz_name   = mms_construct_varname(sc, instr, name, 'oxyz'    );
	ix12_name   = mms_construct_varname(sc, instr, name, 'xi_12'   );
	ix13_name   = mms_construct_varname(sc, instr, name, 'xi_13'   );
	ix23_name   = mms_construct_varname(sc, instr, name, 'xi_23'   );
	
	psf_name    = mms_construct_varname(sc, instr, name, 'psf'     );
	rsf_name    = mms_construct_varname(sc, instr, name, 'rsf'     );
	hstemp_name = mms_construct_varname(sc, instr, name, 'h_stemp' );
	mstemp_name = mms_construct_varname(sc, instr, name, 'm_stemp' );
	etempc_name = mms_construct_varname(sc, instr, name, 'etemp_c' );
	metemp_name = mms_construct_varname(sc, instr, name, 'm_etemp' );
	nla_name    = mms_construct_varname(sc, instr, name, 'nl_a'    );
	nlb_name    = mms_construct_varname(sc, instr, name, 'nl_b'    );
	nlc_name    = mms_construct_varname(sc, instr, name, 'nl_c'    );
	dt13_name   = mms_construct_varname(sc, instr, name, 'dt13'   );
	dt23_name   = mms_construct_varname(sc, instr, name, 'dt23'   );
	
%------------------------------------%
% Read Data                          %
%------------------------------------%
	% Read data calibration data
	%   - Read all of it because there may not be calibration
	%     parameters within the time range yet.
	%   - Also, dTheta and dPhi do not have DEPEND_0 variables.
	%     However, their values appear to change. I have only
	%     pre-flight calibration files so cannot say.
	[gain_data, t_cal] = MrCDF_nRead(filename, gain_name );
	dTheta_data        = MrCDF_nRead(filename, dTheta_name );
	dPhi_data          = MrCDF_nRead(filename, dPhi_name );
	u3_data            = MrCDF_nRead(filename, u3_name );
	offset_data        = MrCDF_nRead(filename, offset_name );
	mpa_data           = MrCDF_nRead(filename, mpa_name );
	gprime_data        = MrCDF_nRead(filename, gprime_name );
	alphas_data        = MrCDF_nRead(filename, alphas_name );
	alphae_data        = MrCDF_nRead(filename, alphae_name );
	stemp_data         = MrCDF_nRead(filename, stemp_name );
	etemp_data         = MrCDF_nRead(filename, etemp_name );
	stempr_data        = MrCDF_nRead(filename, stempr_name );
	etempr_data        = MrCDF_nRead(filename, etempr_name );
	c_data             = MrCDF_nRead(filename, c_name );
	oxyz_data          = MrCDF_nRead(filename, oxyz_name );
	xi12_data          = MrCDF_nRead(filename, ix12_name );
	xi13_data          = MrCDF_nRead(filename, ix13_name );
	xi23_data          = MrCDF_nRead(filename, ix23_name );
	
	% Constants
	psf_data    = MrCDF_nRead(filename, psf_name );
	rsf_data    = MrCDF_nRead(filename, rsf_name );
	hstemp_data = MrCDF_nRead(filename, hstemp_name );
	mstemp_data = MrCDF_nRead(filename, mstemp_name );
	etempc_data = MrCDF_nRead(filename, etempc_name );
	metemp_data = MrCDF_nRead(filename, metemp_name );
	nla_data    = MrCDF_nRead(filename, nla_name );
	nlb_data    = MrCDF_nRead(filename, nlb_name );
	nlc_data    = MrCDF_nRead(filename, nlc_name );
	dt13_data   = MrCDF_nRead(filename, dt13_name );
	dt23_data   = MrCDF_nRead(filename, dt23_name );

%------------------------------------%
% Eliminate Duplicates               %
%------------------------------------%
	% Check for unique times
	[uniq, iuniq] = unique(t_cal);
	if length(uniq) ~= length(t_cal)
		% Warn about eliminating values.
		warning('FG:Read_Cal', 'Duplicate calibration parameters. Taking unique values.')
		
		% Eliminate data
		t_cal       = t_cal(iuniq);
		gain_data   = gain_data(:, iuniq);
		dTheta_data = dTheta_data(:, iuniq);
		dPhi_data   = dPhi_data(:, iuniq);
		u3_data     = u3_data(:, iuniq);
		offset_data = offset_data(:, iuniq);
		mpa_data    = mpa_data(:, iuniq);
		gprime_data = gprime_data(:, iuniq);
		alphas_data = alphas_data(:, iuniq);
		alphae_data = alphae_data(:, iuniq);
		stemp_data  = stemp_data(iuniq);
		etemp_data  = etemp_data(iuniq);
		stempr_data = stempr_data(iuniq);
		etempr_data = etempr_data(iuniq);
		c_data      = c_data(:, iuniq);
		oxyz_data   = oxyz_data(:, iuniq);
		xi12_data   = xi12_data(iuniq);
		xi13_data   = xi13_data(iuniq);
		xi23_data   = xi23_data(iuniq);
	end
	
%------------------------------------%
% Create Output Structure            %
%------------------------------------%
	% Find the closest calibration time and pick those calibration params.
	%   - Transpose to return row vectors instead of column vectors.
	cal_params = struct( 'Epoch',   t_cal,       ...
	                     'Gain',    gain_data,   ...
	                     'dTheta',  dTheta_data, ...
	                     'dPhi',    dPhi_data,   ...
	                     'U3',      u3_data,     ...
	                     'Offset',  offset_data, ...
	                     'MPA',     mpa_data,    ...
	                     'gprime',  gprime_data, ...
	                     'alpha_s', alphas_data, ...
	                     'alpha_e', alphae_data, ...
	                     'stemp',   stemp_data,  ...
	                     'etemp',   etemp_data,  ...
	                     'stemp_r', stempr_data, ...
	                     'etemp_r', etempr_data, ...
	                     'c',       c_data,      ...
	                     'oxyz',    oxyz_data,   ...
	                     'xi_12',   xi12_data,   ...
	                     'xi_13',   xi13_data,   ...
	                     'xi_23',   xi23_data    ...
	                   );
	
	cal_const = struct( 'psf',      psf_data,    ...
	                    'rsf',      rsf_data,    ...
	                    'h_stemp',  hstemp_data, ...
	                    'm_stemp',  mstemp_data, ...
	                    'etempp_c', etempc_data, ...
	                    'm_etemp',  metemp_data, ...
	                    'nl_a',     nla_data,    ...
	                    'nl_b',     nlb_data,    ...
	                    'nl_c',     nlc_data,    ...
	                    'dt13',     dt13_data,   ...
	                    'dt23',     dt23_data    ...
	                  );
end