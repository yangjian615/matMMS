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
function [gd12_dsl, gd21_dsl, gd12_bcs, gd21_bcs] = mms_edi_gse(sc, instr, mode, level, tstart, tend, varargin)

	% Defaults
	edi_dir      = pwd();
	sunpulse_dir = '';

	nOptArgs = length(varargin);
	for ii = 1 : 2 : nOptArgs
		switch varargin{ii}
			case 'SunPulseDir'
				sunpulse_dir   = varargin{ii+1};
			case 'DataDir'
				edi_dir  = varargin{ii+1};
			otherwise
				error( ['Parameter name not recognized: "' varargin{ii} '".'] );
		end
	end
	
	% Assume the calibration files are stored with the data.
	if isempty(sunpulse_dir)
		sunpulse_dir = edi_dir;
	end
	
%------------------------------------%
% Get EDI Data                       %
%------------------------------------%
	% EDI
	[gd12_bcs, gd21_bcs] = mms_edi_bcs(sc, instr, mode, level, tstart, tend, edi_dir);

%------------------------------------%
% Rotate to SMPA                     %
%------------------------------------%
	% Data not available yet


	% Read attitude data
%	[attitude, att_hdr] = mms_fdoa_read_defatt(sc, tstart, tend, att_dir);
%
%	% Build matrix
%	bcs2smpa = mms_fdoa_xbcs2smpa( att_hdr.('zMPA')(:,1)' );
%
%	% Transform
%	fv_gd12_smpa  = mrvector_rotate(bcs2smpa, gd12_bcs.fv_gd12_bcs);
%	gun_gd12_smpa = mrvector_rotate(bcs2smpa, gd12_bcs.gun_gd12_bcs);
%	det_gd12_smpa = mrvector_rotate(bcs2smpa, gd12_bcs.det_gd12_bcs);
%	
%	fv_gd21_smpa  = mrvector_rotate(bcs2smpa, gd12_bcs.fv_gd21_bcs);
%	gun_gd21_smpa = mrvector_rotate(bcs2smpa, gd12_bcs.gun_gd21_bcs);
%	det_gd21_smpa = mrvector_rotate(bcs2smpa, gd12_bcs.det_gd21_bcs);

	fv_gd12_smpa  = gd12_bcs.fv_gd12_bcs;
	gun_gd12_smpa = gd12_bcs.gun_gd12_bcs;
	det_gd12_smpa = gd12_bcs.det_gd12_bcs;

	fv_gd21_smpa  = gd21_bcs.fv_gd21_bcs;
	gun_gd21_smpa = gd21_bcs.gun_gd21_bcs;
	det_gd21_smpa = gd21_bcs.det_gd21_bcs;

%------------------------------------%
% Despin Firing Vectors              %
%------------------------------------%
	% Read sun pulse data
	sunpulse = mms_dss_read_sunpulse(sc, tstart, tend, ...
	                                 'UniquePulse', true, ...
	                                 'Directory', sunpulse_dir);

	% Build matrix
	spin2despun_gd12 = mms_dss_xdespin( sunpulse, gd12_bcs.t_gd12, 'EDI1_GUN' );
	spin2despun_gd21 = mms_dss_xdespin( sunpulse, gd21_bcs.t_gd21, 'EDI2_GUN' );

	% Transform
	fv_gd12_dsl  = mrvector_rotate( spin2despun_gd12, fv_gd12_smpa  );
	fv_gd21_dsl  = mrvector_rotate( spin2despun_gd21, fv_gd21_smpa  );


%------------------------------------%
% Spin Up Gun Positions              %
%------------------------------------%

	% Build rotation matrices
	spin2despun_gun1 = mms_dss_xdespin( sunpulse, gd12_bcs.t_gd12, 'EDI1_GUN',      'SpinUp' );
	spin2despun_det2 = mms_dss_xdespin( sunpulse, gd12_bcs.t_gd12, 'EDI1_DETECTOR', 'SpinUp' );
	
	spin2despun_gun2 = mms_dss_xdespin( sunpulse, gd21_bcs.t_gd21, 'EDI2_GUN',      'SpinUp' );
	spin2despun_det1 = mms_dss_xdespin( sunpulse, gd21_bcs.t_gd21, 'EDI2_DETECTOR', 'SpinUp' );
	
	% Transform
	gun_gd12_dsl = mrvector_rotate( spin2despun_gun1, gun_gd12_smpa );
	det_gd12_dsl = mrvector_rotate( spin2despun_det2, det_gd12_smpa );
	
	gun_gd21_dsl = mrvector_rotate( spin2despun_gun2, gun_gd21_smpa );
	det_gd21_dsl = mrvector_rotate( spin2despun_det1, det_gd21_smpa );

%------------------------------------%
% Rotate to GSE                      %
%------------------------------------%
	% GEI -> Despun
%	gei2despun_gd12 = mms_fdoa_xgei2despun(attitude, edi1_bcs.t_gd12, 'L');
%	gei2despun_gd21 = mms_fdoa_xgei2despun(attitude, edi1_bcs.t_gd21, 'L');
%	
%	% Despun -> GEI
%	despun2gei_gd12 = permute(gei2despun_gd12, [2, 1, 3]);
%	despun2gei_gd21 = permute(gei2despun_gd21, [2, 1, 3]);
%	
%	%
%	% Transform firing vectors and positions
%	%
%	
%	% GD12
%	fv_gd12_gei  = mrvector_rotate(despun2gei_gd12, fv_gd12_dsl);
%	gun_gd12_gei = mrvector_rotate(despun2gei_gd12, gun_gd12_dsl);
%	det_gd12_gei = mrvector_rotate(despun2gei_gd12, det_gd12_dsl);
%	
%	% GD21
%	fv_gd21_gei  = mrvector_rotate(despun2gei_gd21, fv_gd21_dsl);
%	gun_gd21_gei = mrvector_rotate(despun2gei_gd21, gun_gd21_dsl);
%	det_gd21_gei = mrvector_rotate(despun2gei_gd21, det_gd21_dsl);
%
%	%
%	% Transformation matrix GEI -> GSE
%	%   - Modified Julian Date (mjd).
%	%   - UTC seconds since midnight (ssm).
%	%
%	
%	% GD12
%	timevec = MrCDF_Epoch_Breakdown( gd12_bcs.t_gd12 );
%	mjd     = date2mjd(timevec(1,:), timevec(2,:), timevec(3,:));
%	ssm = timevec(4,:) * 3600.0 + ...
%		  timevec(5,:) * 60.0   + ...
%		  timevec(6,:)          + ...
%		  timevec(7,:) * 1e-3   + ...
%		  timevec(8,:) * 1e-6   + ...
%		  timevec(9,:) * 1e-9;
%	GEI2GSE_gd12 = gei2gse(mjd, ssm);
%	
%	% GD21
%	timevec = MrCDF_Epoch_Breakdown( gd21_bcs.t_gd21 );
%	mjd     = date2mjd(timevec(1,:), timevec(2,:), timevec(3,:));
%	ssm = timevec(4,:) * 3600.0 + ...
%		  timevec(5,:) * 60.0   + ...
%		  timevec(6,:)          + ...
%		  timevec(7,:) * 1e-3   + ...
%		  timevec(8,:) * 1e-6   + ...
%		  timevec(9,:) * 1e-9;
%	GEI2GSE_gd21 = gei2gse(mjd, ssm);
%
%	%
%	% Transform to GSE
%	%
%	
%	% GD12
%	fv_gd12_gse  = mrvector_rotate(GEI2GSE_gd12, fv_gd12_gei);
%	gun_gd12_gse = mrvector_rotate(GEI2GSE_gd12, gun_gd12_gei);
%	det_gd12_gse = mrvector_rotate(GEI2GSE_gd12, det_gd12_gei);

%------------------------------------%
% Return                             %
%------------------------------------%
	
	% EDI1 output structure
	%   - EDI1 contains gun1 and detector2
	gd12_dsl = struct( 't_gd12',        gd12_bcs.t_gd12,   ...
	                   'gun_gd12_dsl',  gun_gd12_dsl, ...
	                   'det_gd12_dsl',  det_gd12_dsl, ...
	                   'fv_gd12_dsl',   fv_gd12_dsl );
	
	% EDI2 output structure
	%   - EDI2 contains gun1 and detector2
	gd21_dsl = struct( 't_gd21',        gd21_bcs.t_gd21,   ...
	                   'gun_gd21_dsl',  gun_gd21_dsl, ...
	                   'det_gd21_dsl',  det_gd21_dsl, ...
	                   'fv_gd21_dsl',   fv_gd21_dsl );


end