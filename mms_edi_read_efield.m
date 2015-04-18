%
% Name
%   mms_edi_read_efield
%
% Purpose
%   Read MMS EDI electric field mode data files.
%
% Calling Sequence
%   EDI_STRUCT = mms_edi_read_efield(SC, INSTR, MODE, LEVEL, TSTART, TEND, EDI_DIR)
%     Read EDI electric field mode data captured by spacecraft SC
%     (e.g. 'mms3'), instrument INSTR (e.g. 'edi'), from telemetry
%     mode MODE and data product level LEVEL between the time
%     interval of [TSTART, TEND]. Data can be found in directory
%     EDI_DIR. Times should be provided in ISO format:
%     'yyyy-mm-ddThh:mm_ss'.
%
% Parameters
%   SC              in, required, type=char/cell
%   INSTR           in, required, type=char
%   MODE            in, required, type=char
%   LEVEL           in, required, type=char
%   TSTART          in, required, type=char
%   TEND            in, required, type=char
%   EDI_DIR         in, required, type=char
%
% Returns
%   EDI_STRUCT      out, required, type=structure
%                   Fields are:
%                     'epoch_gd12' -  TT2000 Epoch time for gun 1 and detector 2.
%                     'epoch_gd21' -  TT2000 Epoch time for gun 1 and detector 2.
%                     'g1aV'       -  Firing angle for gun 1.
%                     'g2aV'       -  Firing angle for gun 2.
%                     'g1pos'      -  Position of gun 1 in EDI_GUN system
%                     'g2pos'      -  Position of gun 2 in EDI_GUN system
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-18      Written by Matthew Argall
%
function [t_gd12, t_gd21, fa_gd12, fa_gd21] = mms_edi_read_efiled(sc, instr, mode, level, tstart, tend, edi_dir)
	
	if nargin < 7
		edi_dir = pwd();
	end
	
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
	
	% Gun positions with respect to GDU1 coordinate system
	g1pos_x_vname = mms_construct_varname(sc, instr, 'v1xg1');
	g1pos_y_vname = mms_construct_varname(sc, instr, 'v1xg1');
	g1pos_z_vname = mms_construct_varname(sc, instr, 'v1xg1');
	g2pos_x_vname = mms_construct_varname(sc, instr, 'v1xg1');
	g2pos_y_vname = mms_construct_varname(sc, instr, 'v1xg1');
	g2pos_z_vname = mms_construct_varname(sc, instr, 'v1xg1');
	g2pos_z_vname = mms_construct_varname(sc, instr, 'v1xg1');
	
	% Gun analog voltages
	%   - My version of mms_edi_avoltage2angles does not convert properly
	%   - Read Theta and Phi angles instead
%	g1aV_x_vname = mms_construct_varname(sc, instr, 'vax_gd12');
%	g1aV_y_vname = mms_construct_varname(sc, instr, 'vay_gd12');
%	g2aV_x_vname = mms_construct_varname(sc, instr, 'vax_gd21');
%	g2aV_y_vname = mms_construct_varname(sc, instr, 'vay_gd21');
	
	% Gun firing angles in EDI1 EDI2 coordinates
	g1_theta_vname = mms_construct_varname(sc, instr, 'theta_gd12');
	g2_theta_vname = mms_construct_varname(sc, instr, 'theta_gd21');
	g1_phi_vname   = mms_construct_varname(sc, instr, 'phi_gd12');
	g2_phi_vname   = mms_construct_varname(sc, instr, 'phi_gd21');
	
%------------------------------------%
% Read Data                          %
%------------------------------------%
	%
	% Firing vectors in [x,y,z] are not filled yet.
	% Convert analog voltages to firing vectors outside.
	%

	[g1_theta, t_gd12] = MrCDF_Read(fname, g1_theta_vname,  'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	g1_phi             = MrCDF_Read(fname, g1_phi_vname,  'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
%	g1pos_x          = MrCDF_Read(fname, g1pos_x_vname, 'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
%	g1pos_y          = MrCDF_Read(fname, g1pos_y_vname, 'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
%	g1pos_z          = MrCDF_Read(fname, g1pos_z_vname, 'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	[g2_theta, t_gd21] = MrCDF_Read(fname, g2_theta_vname,  'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	g2_phi             = MrCDF_Read(fname, g2_phi_vname,  'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
%	g2pos_x          = MrCDF_Read(fname, g2pos_x_vname, 'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
%	g2pos_y          = MrCDF_Read(fname, g2pos_y_vname, 'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
%	g2pos_z          = MrCDF_Read(fname, g2pos_z_vname, 'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	
	% Combine into 3- and 2-vectors to get firing vectors
%	g1pos = [g1pos_x; g1pos_y; g1pos_z];
%	g2pos = [g2pos_x; g2pos_y; g2pos_z];

	% Firing angles
	%   - Phi is the angle measured from x, in the xy-plane, positive toward y
	%   - Theta is the angle measured down from the z-axis
	%
	% From Hans Vaith
	%   gx = cos(phi) * sin(theta)
	%   gy = sin(phi) * sin(theta)
	%   gz = cos(theta)
	%
	% Order as [phi, theta] to pass directly to sph2cart
	fa_gd12 = [g1_phi; g1_theta];
	fa_gd21 = [g2_phi; g2_theta];
	
%------------------------------------%
% Create the Output Structure        %
%------------------------------------%
%	edi_struct = struct( 'epoch_gd12', epoch_gd12, ...
%	                     'epoch_gd21', epoch_gd21, ...
%	                     'g1aV',   g1aV,           ...
%	                     'g2aV',   g2aV            ...
%	                   );
%	                     'g1pos',  g1pos,      ...
%	                     'g2pos',  g2pos );
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