%
% Name
%   mms_edi_read_efield
%
% Purpose
%   Read MMS EDI electric field mode data files.
%
% Calling Sequence
%   EDI_L1A_EMODE = mms_edi_read_efield(FILENAMES)
%     Read EDI electric field mode data from files FILENAMES.
%
%   [__] = mms_edi_read_efield(__, TSTART, TEND)
%     Select data between the time interval beginning at TSTART and
%     ending at TEND. TSTART and TEND should be formatted as
%     'yyyy-mm-ddThh:mm:ss'.
%
%   [__] = mms_edi_read_efield(..., 'Quality', quality)
%     Select data between the time interval beginning at TSTART and
%     ending at TEND. TSTART and TEND should be formatted as
%     'yyyy-mm-ddThh:mm:ss'.
%
% Parameters
%   FILENAMES       in, required, type=char/cell
%   TSTART          in, required, type=char, default=''
%   TEND            in, required, type=char, default=''
%   'Quality'       in, optional, type=intarr, default=[]
%
% Returns
%   EDI_L1A_EMODE   out, required, type=structure
%                   Fields are:
%                     'epoch_gd12'      -  TT2000 Epoch time.
%                     'polar_gd12'      -  polar firing angle.
%                     'azimuth_gd12'    -  azimuth firing angles.
%                     'quality_gd12'    -  Quality flag.
%                     'tof_gd12'        -  Time of flight.
%                     'energy_gd12'     -  Beam energy.
%                     'num_chips_gd12'  -  Number of chips per code length.
%                     'm_gd12'          -  Correlator chip length.
%                     'n_gd12'          -  Correlator chip length.
%                     'max_addr_gd12'   -  Address of maximum beam correlation.
%
%                     'epoch_gd21'      -  TT2000 Epoch time.
%                     'polar_gd21'      -  polar firing angle.
%                     'azimuth_gd21'    -  azimuth firing angles.
%                     'q_gd21'          -  Quality flag.
%                     'tof_gd21'        -  Time of flight.
%                     'energy_gd21'     -  Beam energy.
%                     'num_chips_gd21'  -  Number of chips per code length.
%                     'm_gd21'          -  Correlator chip length.
%                     'n_gd21'          -  Correlator chip length.
%                     'max_addr_gd21'   -  Address of maximum beam correlation.
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-18      Written by Matthew Argall
%   2015-04-23      Cleaned up code. Return structure for each gun.
%   2015-05-25      Take file names as inputs. Read more data from files. Calculate
%                     firing vectors. Filter by quality. Return a single structure. - MRA
%   2015-06-03      spdfcdfread does not accept CDF IDs, so cannot pre-open files. - MRA
%
function edi_l1a_emode = mms_edi_read_efield(filenames, tstart, tend, varargin)

	% Defaults
	quality = [];
	if nargin() < 2
		tstart = '';
	end
	if nargin() < 3
		tend = '';
	end
	if nargin() == 5
		assert( strcmp(varargin{1}, 'Quality'), ['Unknown parameter name: "' varargin{1} '".'] );
		quality = varargin{2};
	end
	
	% Valid quality specification?
	if ~isempty(quality)
		if length(quality) > 4 || max(quality) > 3 || min(quality) < 0
			error('Quality can have only these unique values: 0, 1, 2, 3.');
		end
	end

%------------------------------------%
% Check Files                        %
%------------------------------------%
	% Number of files given
	if iscell(filenames)
		nFiles = length(filenames);
	else
		nFiles = 1;
	end

	% Dissect the file name
	[sc, instr, mode, level, ~, ~, optdesc] = mms_dissect_filename(filenames);

	% AFG or DFG
	assert( min( strcmp(instr, 'edi') ) == 1, 'Only EDI files are allowed.' );

	% Instr, Level, Mode
	if nFiles > 1
		assert( min( strcmp(mode, mode{1}) ) == 1, 'All files must have the same telemetry mode.' );
	else
	assert( min( strcmp(level,   'l1a') )    == 1, 'Only L1A files are allowed.' );
	assert( min( strcmp(optdesc, 'efield') ) == 1, 'Only E-field Mode files are allowed.' );

	% We now know all the files match, so keep on the the first value.
	if nFiles > 1
		sc      = sc{1};
		instr   = instr{1};
		mode    = mode{1};
		level   = level{1};
		optdesc = optdesc{1};
	end
	
%------------------------------------%
% File and Variable Names            %
%------------------------------------%
	
	% GDU1 variable names
	phi_gd12_vname      = mms_construct_varname(sc, instr, 'phi_gd12');
	theta_gd12_vname    = mms_construct_varname(sc, instr, 'theta_gd12');
	q_gd12_vname        = mms_construct_varname(sc, instr, 'sq_gd12');
	tof_gd12_vname      = mms_construct_varname(sc, instr, 'tof1_us');
	e_gd12_name         = mms_construct_varname(sc, instr, 'e_gd12');
	num_chips_gd12_name = mms_construct_varname(sc, instr, 'numchips_gd12');
	m_gd12_name         = mms_construct_varname(sc, instr, 'm_gd12');
	n_gd12_name         = mms_construct_varname(sc, instr, 'n_gd12');
	max_addr_gd12_name  = mms_construct_varname(sc, instr, 'max_addr_gd12');
	
	% GDU2 variable names
	phi_gd21_vname      = mms_construct_varname(sc, instr, 'phi_gd21');
	theta_gd21_vname    = mms_construct_varname(sc, instr, 'theta_gd21');
	q_gd21_vname        = mms_construct_varname(sc, instr, 'sq_gd21');
	tof_gd21_vname      = mms_construct_varname(sc, instr, 'tof2_us');
	e_gd21_name         = mms_construct_varname(sc, instr, 'e_gd21');
	num_chips_gd21_name = mms_construct_varname(sc, instr, 'numchips_gd21');
	m_gd21_name         = mms_construct_varname(sc, instr, 'm_gd21');
	n_gd21_name         = mms_construct_varname(sc, instr, 'n_gd21');
	max_addr_gd21_name  = mms_construct_varname(sc, instr, 'max_addr_gd21');
	
%------------------------------------%
% Read Data                          %
%------------------------------------%

	%
	% Firing vectors in [x,y,z] are not filled yet.
	% Convert analog voltages to firing vectors outside.
	%
	[theta_gd12, epoch_gd12] = MrCDF_nRead(filenames, theta_gd12_vname,    'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	phi_gd12                 = MrCDF_nRead(filenames, phi_gd12_vname,      'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	q_gd12                   = MrCDF_nRead(filenames, q_gd12_vname,        'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	tof_gd12                 = MrCDF_nRead(filenames, tof_gd12_vname,      'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	e_gd12                   = MrCDF_nRead(filenames, e_gd12_name,         'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	num_chips_gd12           = MrCDF_nRead(filenames, num_chips_gd12_name, 'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	m_gd12                   = MrCDF_nRead(filenames, m_gd12_name,         'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	n_gd12                   = MrCDF_nRead(filenames, n_gd12_name,         'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	max_addr_gd12            = MrCDF_nRead(filenames, max_addr_gd12_name,  'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);

	[theta_gd21, epoch_gd21] = MrCDF_nRead(filenames, theta_gd21_vname,    'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	phi_gd21                 = MrCDF_nRead(filenames, phi_gd21_vname,      'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	q_gd21                   = MrCDF_nRead(filenames, q_gd21_vname,        'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	tof_gd21                 = MrCDF_nRead(filenames, tof_gd21_vname,      'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	e_gd21                   = MrCDF_nRead(filenames, e_gd21_name,         'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	num_chips_gd21           = MrCDF_nRead(filenames, num_chips_gd21_name, 'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	m_gd21                   = MrCDF_nRead(filenames, m_gd21_name,         'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	n_gd21                   = MrCDF_nRead(filenames, n_gd21_name,         'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	max_addr_gd21            = MrCDF_nRead(filenames, max_addr_gd21_name,  'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);


	%
	% Firing angles
	%   - Theta is the polar angle
	%   - Phi is the azimuth angle
	%
	% From Hans Vaith
	%   gx = cos(phi) * sin(theta)
	%   gy = sin(phi) * sin(theta)
	%   gz = cos(theta)
	%
	
%------------------------------------%
% Filter by Quality                  %
%------------------------------------%
	if ~isempty(quality)
		if length(quality) > 4 || max(quality) > 3 || min(quality) < 0
			error('Quality can have only these unique values: 0, 1, 2, 3.');
		end
		
		% Find quality
		iq_gd12 = ismember(gd12.q_gd12, quality);
		iq_gd21 = ismember(gd21.q_gd21, quality);
		
		% Number of beams of the selected quality
		nq_gd12 = length(iq_gd12);
		nq_gd21 = length(iq_gd21);

		% Make sure something was found
		if nq_gd12 == 0 && nq_gd21 == 0
			error('EDI_BCS:Quality', 'No beams of the selected quality found.');
		elseif sum(iq_gd12) == 0
			warning('EDI_BCS:Quality', 'No beams of selected quality found for Gun1');
		elseif sum(iq_gd21) == 0
			warning('EDI_BCS:Quality', 'No beams of selected quality found for Gun2');
		end
		
		% Select data
		epoch_gd12     = epoch_gd12(iq_gd12);
		q_gd12         = q_gd12(iq_gd12);
		tof_gd12       = tof_gd12(iq_gd12);
		e_gd12         = e_gd12(iq_gd12);
		num_chips_gd12 = num_chips_gd12(iq_gd12);
		m_gd12         = m_gd12(iq_gd12);
		n_gd12         = n_gd12(iq_gd12);
		max_addr_gd12  = max_addr_gd12(iq_gd12);
		
		epoch_gd21     = epoch_gd21(iq_gd21);
		q_gd21         = q_gd21(iq_gd21);
		tof_gd21       = tof_gd21(iq_gd21);
		e_gd21         = e_gd21(iq_gd21);
		num_chips_gd21 = num_chips_gd21(iq_gd21);
		m_gd21         = m_gd21(iq_gd21);
		n_gd21         = n_gd21(iq_gd21);
		max_addr_gd21  = max_addr_gd21(iq_gd21);
	end

%------------------------------------%
% Firing Vectors                     %
%------------------------------------%
	% Convert to radians
	deg2rad      = pi / 180.0;
	polar_gd12   = theta_gd12 * deg2rad;
	azimuth_gd12 = phi_gd12   * deg2rad;
	polar_gd21   = theta_gd21 * deg2rad;
	azimuth_gd21 = phi_gd21   * deg2rad;

	% Convert to cartesian coordinates
	%   - sph2cart requires the elevation angle, measured up from the xy-plane,
	%     not down from the z-axis.
	fv_gd12      = zeros( [3, length(polar_gd12)] );
	fv_gd12(1,:) = sin( polar_gd12 ) .* cos( azimuth_gd12 );
	fv_gd12(2,:) = sin( polar_gd12 ) .* sin( azimuth_gd12 );
	fv_gd12(3,:) = cos( polar_gd12 );
	
	fv_gd21      = zeros( [3, length(polar_gd21)] );
	fv_gd21(1,:) = sin( polar_gd21 ) .* cos( azimuth_gd21 );
	fv_gd21(2,:) = sin( polar_gd21 ) .* sin( azimuth_gd21 );
	fv_gd21(3,:) = cos( polar_gd21 );

%------------------------------------%
% Create the Output Structure        %
%------------------------------------%
	edi_l1a_emode = struct( 'epoch_gd12',     epoch_gd12,     ...      GD12
	                        'polar_gd12',     theta_gd12,     ...
	                        'azimuth_gd12',   phi_gd12,       ...
	                        'fv_gd12_123',    fv_gd12,        ...
	                        'quality_gd12',   q_gd12,         ...
	                        'tof_gd12',       tof_gd12,       ...
	                        'energy_gd12',    e_gd12,         ...
	                        'num_chips_gd12', num_chips_gd12, ...
	                        'm_gd12',         m_gd12,         ...
	                        'n_gd12',         n_gd12,         ...
	                        'max_addr_gd12',  max_addr_gd12,  ...
	                        'epoch_gd21',     epoch_gd21,     ...      GD21
	                        'polar_gd21',     theta_gd21,     ...
	                        'azimuth_gd21',   phi_gd21,       ...
	                        'fv_gd21_123',    fv_gd21,        ...
	                        'quality_gd21',   q_gd21,         ...
	                        'tof_gd21',       tof_gd21,       ...
	                        'energy_gd21',    e_gd21,         ...
	                        'num_chips_gd21', num_chips_gd21, ...
	                        'm_gd21',         m_gd21,         ...
	                        'n_gd21',         n_gd21,         ...
	                        'max_addr_gd21',  max_addr_gd21 );
end