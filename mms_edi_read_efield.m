%
% Name
%   mms_edi_read_efield
%
% Purpose
%   Read MMS EDI electric field mode data files.
%
% Calling Sequence
%   [GD12, GD21] = mms_edi_read_efield(SC, INSTR, MODE, LEVEL, TSTART, TEND)
%     Read EDI electric field mode data captured by spacecraft SC
%     (e.g. 'mms3'), instrument INSTR (e.g. 'edi'), from telemetry
%     mode MODE and data product level LEVEL between the time
%     interval of [TSTART, TEND]. Data is search for in the present
%     working directory. Times should be provided in ISO format:
%     'yyyy-mm-ddThh:mm_ss'.
%
%   [__] = mms_edi_read_efield(__, edi_dir)
%     Specify the directory in which to look for EDI data.
%
% Parameters
%   SC              in, required, type=char/cell
%   INSTR           in, required, type=char
%   MODE            in, required, type=char
%   LEVEL           in, required, type=char
%   TSTART          in, required, type=char
%   TEND            in, required, type=char
%   EDI_DIR         in, optional, type=char
%
% Returns
%   gd12            out, required, type=structure
%                   Fields are:
%                     'epoch_gd12'    -  TT2000 Epoch time.
%                     'polar_gd12'    -  polar firing angle.
%                     'azimuth_gd12'  -  azimuth firing angles.
%                     'q_gd12'        -  Quality flag.
%                     'tof_gd12'      -  Time of flight.
%   gd12            out, required, type=structure
%                   Fields are:
%                     'epoch_gd21'    -  TT2000 Epoch time.
%                     'polar_gd21'    -  polar firing angle.
%                     'azimuth_gd21'  -  azimuth firing angles.
%                     'q_gd21'        -  Quality flag.
%                     'tof_gd21'      -  Time of flight.
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-18      Written by Matthew Argall
%   2015-04-23      Cleaned up code. Return structure for each gun.
%
function [gd12, gd21] = mms_edi_read_efield(sc, instr, mode, level, tstart, tend, edi_dir)
	
%------------------------------------%
% File and Variable Names            %
%------------------------------------%
	% Create the file name
	fpattern = mms_construct_filename(sc, instr, mode, level, ...
	                                  'Directory', edi_dir, ...
	                                  'Tokens',    true);

	% Search for the file
	fname = MrFile_Search(fpattern,              ...
	                      'TStart',    tstart,   ...
	                      'TEnd',      tend,     ...
	                      'TimeOrder', '%Y%M%d', ...
	                      'Closest',   true);
	
	% GDU1 variable names
	phi_gd12_vname   = mms_construct_varname(sc, instr, 'phi_gd12');
	theta_gd12_vname = mms_construct_varname(sc, instr, 'theta_gd12');
	q_gd12_vname     = mms_construct_varname(sc, instr, 'sq_gd12');
	tof_gd12_vname   = mms_construct_varname(sc, instr, 'tof1_us');
	
	% GDU2 variable names
	phi_gd21_vname   = mms_construct_varname(sc, instr, 'phi_gd21');
	theta_gd21_vname = mms_construct_varname(sc, instr, 'theta_gd21');
	q_gd21_vname     = mms_construct_varname(sc, instr, 'sq_gd21');
	tof_gd21_vname   = mms_construct_varname(sc, instr, 'tof2_us');
	
%------------------------------------%
% Read Data                          %
%------------------------------------%
	%
	% Firing vectors in [x,y,z] are not filled yet.
	% Convert analog voltages to firing vectors outside.
	%
	[theta_gd12, epoch_gd12] = MrCDF_Read(fname, theta_gd12_vname,  'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	phi_gd12                 = MrCDF_Read(fname, phi_gd12_vname,    'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	q_gd12                   = MrCDF_Read(fname, q_gd12_vname,      'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	tof_gd12                 = MrCDF_Read(fname, tof_gd12_vname,    'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);

	[theta_gd21, epoch_gd21] = MrCDF_Read(fname, theta_gd21_vname,  'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	phi_gd21                 = MrCDF_Read(fname, theta_gd21_vname,  'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	q_gd21                   = MrCDF_Read(fname, q_gd21_vname,      'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	tof_gd21                 = MrCDF_Read(fname, tof_gd21_vname,    'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);

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
% Create the Output Structure        %
%------------------------------------%
	gd12 = struct( 'epoch_gd12',    epoch_gd12, ...
	               'polar_gd12',    theta_gd12, ...
	               'azimuth_gd12',  phi_gd12,   ...
	               'q_gd12',        q_gd12,     ...
	               'tof_gd12',      tof_gd12 );
	
	gd21 = struct( 'epoch_gd21',    epoch_gd21, ...
	               'polar_gd21',    theta_gd21, ...
	               'azimuth_gd21',  phi_gd21,   ...
	               'q_gd21',        q_gd21,     ...
	               'tof_gd21',      tof_gd21 );
end



% Epoch
% epoch_timetag
% delta_t
% Epoch_beam_gd12
% Epoch_beam_gd21
% mms3_edi_t1
% mms3_edi_t1m2
% mms3_edi_tof_over
% mms3_edi_tof1_us
% mms3_edi_tof2_us
% mms3_edi_dtof_p_over
% mms3_edi_dtof_n_over
% mms3_edi_vax_gd12
% mms3_edi_vax_gd21
% mms3_edi_vay_gd12
% mms3_edi_vay_gd21
% mms3_edi_theta_gd12
% mms3_edi_theta_gd21
% mms3_edi_phi_gd12
% mms3_edi_phi_gd21
% mms3_edi_v1xg1
% mms3_edi_v1yg1
% mms3_edi_v1zg1
% mms3_edi_v2xg1
% mms3_edi_v2yg1
% mms3_edi_v2zg1
% mms3_edi_word14_gd12
% mms3_edi_word14_gd21
% mms3_edi_word15_gd12
% mms3_edi_word15_gd21
% mms3_edi_codelength_gd21
% mms3_edi_codelength_gd12
% mms3_edi_e
% mms3_edi_sq_gd12
% mms3_edi_sq_gd21
% mms3_edi_i1
% mms3_edi_i2
% mms3_edi_m_gd21
% mms3_edi_m_gd12
% mms3_edi_n_gd21
% mms3_edi_n_gd12
% mms3_edi_max_addr_gd12
% mms3_edi_max_addr_gd21
% epoch_crsf
% crsf_10min