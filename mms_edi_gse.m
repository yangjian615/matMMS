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
%                     'g1fa'       -  Firing angle for gun 1.
%                     'g2fa'       -  Firing angle for gun 2.
%                     'g1pos'      -  Position of gun 1 in EDI_GUN system
%                     'g2pos'      -  Position of gun 2 in EDI_GUN system
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-18      Written by Matthew Argall
%
function edi_struct = mms_edi_bcs(sc, instr, mode, level, tstart, tend, edi_dir)

	% Inputs
	edi_data_dir = '/Users/argall/Documents/Work/Data/MMS/EDI/';
	dfg_data_dir = '/Users/argall/Documents/Work/Data/MMS/DFG/';
	fg_cal_dir   = '/Users/argall/Documents/Work/Data/MMS/FG_Cal/';
	att_dir      = '/Users/argall/Documents/Work/Data/MMS/Attitude/';
	sc           = 'mms3';
	instr        = 'edi';
	mode         = 'slow';
	level        = 'l1a_efield';
	tstart       = '2015-04-18T00:00:00';
	tend         = '2015-04-18T24:00:00';
	
%------------------------------------%
% Get EDI Data                       %
%------------------------------------%
	% EDI
	[edi1_bcs, edi2_bcs] = mms_edi_bcs(sc, instr, mode, level, tstart, tend, edi_data_dir);
	
	% Convert firing voltages to cartesian coordinates
	g1fa_edi1 = mms_edi_avoltage2angle( g1aV(1,:), g1aV(2,:), 'Cartesian' );
	g2fa_edi2 = mms_edi_avoltage2angle( g2aV(1,:), g2aV(2,:), 'Cartesian' );
	
%------------------------------------%
% Transform to SMPA                  %
%------------------------------------%
	% Read attitude data
	[attitude, att_hdr] = mms_fdoa_read_defatt(sc, tstart, tend, att_dir);

	% Build matrix
	bcs2smpa = mms_fdoa_xbcs2smpa( att_hdr.('zMPA')(:,1)' );

	% Transform
	fv_gd12_smpa = mrvector_rotate( bcs2smpa, edi1_bcs.fv_gd12_bcs );
	fv_gd21_smpa = mrvector_rotate( bcs2smpa, edi2_bcs.fv_gd21_bcs );
	pos_g1_smpa  = mrvector_rotate( bcs2smpa, edi1_bcs.pos_g1_bcs  );
	pos_g2_smpa  = mrvector_rotate( bcs2smpa, edi2_bcs.pos_g2_bcs  );
	pos_d1_smpa  = mrvector_rotate( bcs2smpa, edi1_bcs.pos_d1_bcs  );
	pos_d2_smpa  = mrvector_rotate( bcs2smpa, edi2_bcs.pos_d2_bcs  );

%------------------------------------%
% Despin Firing Vectors              %
%------------------------------------%

	%
	% Assume the principle axis of inertia (z-MPA)
	% is the same as the angular momentum vector (L)
	%
	smpa2dsl_g1 = mms_fdoa_xdespin(attitude, edi1_bcs.t_gd12, 'L', [ 'EDI1_GUN' ]);
	smpa2dsl_g2 = mms_fdoa_xdespin(attitude, edi1_bcs.t_gd21, 'L', [ 'EDI2_GUN' ]);

	% Rotate to DSL
	fa_gd12_dsl = mrvector_rotate(smpa2dsl_g1, fv_gd12_smpa);
	fa_gd21_dsl = mrvector_rotate(smpa2dsl_g2, fv_gd21_smpa);

%------------------------------------%
% Spin Up Firing Vectors             %
%------------------------------------%
	%
	% Assume the principle axis of inertia (z-MPA)
	% is the same as the angular momentum vector (L)
	%
	gun1_smpa2dsl = mms_fdoa_xdespin(attitude, t, 'L', [ 'EDI1_GUN' ],      'SpinUp');
	gun2_smpa2dsl = mms_fdoa_xdespin(attitude, t, 'L', [ 'EDI2_GUN' ],      'SpinUp');
	det1_smpa2dsl = mms_fdoa_xdespin(attitude, t, 'L', [ 'EDI1_DETECTOR' ], 'SpinUp');
	det2_smpa2dsl = mms_fdoa_xdespin(attitude, t, 'L', [ 'EDI2_DETECTOR' ], 'SpinUp');

	% Rotate to DSL
	gun1_dsl = mrvector_rotate(gun1_smpa2dsl, pos_g1_smpa);
	gun2_dsl = mrvector_rotate(gun2_smpa2dsl, pos_g2_smpa);
	det1_dsl = mrvector_rotate(det1_smpa2dsl, pos_d1_smpa);
	det2_dsl = mrvector_rotate(det2_smpa2dsl, pos_d2_smpa);

%------------------------------------%
% Rotate to GSE                      %
%------------------------------------%
	% Transform into GEI
	gei2despun = mms_fdoa_xgei2despun(attitude, t, 'L');
	despun2gei = permute(gei2despun, [2, 1, 3]);
	b_gei      = mrvector_rotate(despun2gei, b_dsl);

	% Transform matrix GEI -> GSE
	%   - Modified Julian Date (mjd).
	%   - UTC seconds since midnight (ssm).
	timevec = MrCDF_Epoch_Breakdown(t);
	mjd     = date2mjd(timevec(1,:), timevec(2,:), timevec(3,:));
	ssm = timevec(4,:) * 3600.0 + ...
		  timevec(5,:) * 60.0   + ...
		  timevec(6,:)          + ...
		  timevec(7,:) * 1e-3   + ...
		  timevec(8,:) * 1e-6   + ...
		  timevec(9,:) * 1e-9;
	GEI2GSE = gei2gse(mjd, ssm);

	% Transform to GSE
	b_gse = mrvector_rotate(GEI2GSE, b_gei);


end