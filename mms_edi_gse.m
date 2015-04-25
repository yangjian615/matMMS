%
% Name
%   mms_edi_read_efield
%
% Purpose
%   Read MMS EDI electric field mode data files.
%
% Calling Sequence
%   [GD12_GSE, GD21_GSE] = mms_edi_read_efield(SC, INSTR, MODE, LEVEL, TSTART, TEND)
%     Read EDI electric field mode data captured by spacecraft SC
%     (e.g. 'mms3'), instrument INSTR (e.g. 'edi'), from telemetry
%     mode MODE and data product level LEVEL between the time
%     interval of [TSTART, TEND]. Times should be provided in ISO format:
%     'yyyy-mm-ddThh:mm_ss'.
%
%   [GD12_GSE, GD21_GSE] = mms_edi_read_efield(__, 'ParamName', ParamValue)
%     Include any of the parameter name-value pairs below.
%
%   [..., GD12_DSL, GD21_DSL] = mms_edi_read_efield(__)
%     Return gun and detector information in despun coordinates.
%
%   [..., GD12_BCS, GD21_BCS] = mms_edi_read_efield(__)
%     Return gun and detector information in the spinning body coordinate system (BCS).
%
% Parameters
%   SC              in, required, type=char/cell
%   INSTR           in, required, type=char
%   MODE            in, required, type=char
%   LEVEL           in, required, type=char
%   TSTART          in, required, type=char
%   TEND            in, required, type=char
%   'AttDir'        in, optional, type=char, default='DataDir'
%                   Directory in which to find definitive attitude data.
%   'DataDir'       in, optional, type=char, default=pwd()
%                   Directory in which to find EDI data.
%   'SunPulseDir'   in, optional, type=char, default='DataDir'
%                   Directory in which to find sun pulse data. If provided,
%                     sun pulse times will be used to depsin data. Otherwise,
%                     spin phase data from the devinitive attitude files is
%                     used.
%
% Returns
%   GD12_DSL        out, required, type=structure
%                   Fields are:
%                     't_gd12'       -  TT2000 Epoch time for gun 1 and detector 2.
%                     'gun_gd12_dsl' -  Position of the gun for gd12 pair.
%                     'det_gd12_dsl' -  Position of the detector for gd12 pair.
%                     'fv_gd12_dsl'  -  Firing vector for gd12 pair.
%   GD21_DSL        out, required, type=structure
%                   Fields are:
%                     't_gd21'       -  TT2000 Epoch time for gun 2 and detector 1.
%                     'gun_gd21_dsl' -  Position of the gun for gd21 pair.
%                     'det_gd21_dsl' -  Position of the detector for gd21 pair.
%                     'fv_gd21_dsl'  -  Firing vector for gd21 pair.
%   GD12_BCS        out, optional, type=structure
%                   Fields are:
%                     't_gd12'       -  TT2000 Epoch time for gun 1 and detector 2.
%                     'gun_gd12_bcs' -  Position of the gun for gd12 pair.
%                     'det_gd12_bcs' -  Position of the detector for gd12 pair.
%                     'fv_gd12_bcs'  -  Firing vector for gd12 pair.
%   GD21_BCS        out, optional, type=structure
%                   Fields are:
%                     't_gd21'       -  TT2000 Epoch time for gun 2 and detector 1.
%                     'gun_gd21_bcs' -  Position of the gun for gd21 pair.
%                     'det_gd21_bcs' -  Position of the detector for gd21 pair.
%                     'fv_gd21_bcs'  -  Firing vector for gd21 pair.
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-19      Written by Matthew Argall
%   2015-04-20      Data is despun with attitude data unless SunPulseDir is given. - MRA
%
function [gd12_dmpa, gd21_dmpa, gd12_bcs, gd21_bcs] = mms_edi_gse(sc, instr, mode, level, tstart, tend, varargin)

	% Defaults
	edi_dir      = pwd();
	sunpulse_dir = '';
	att_dir      = '';

	% 'AttDir' and 'SunPulseDir' are used here. The rest of the optional
	% parameters are passed on to mms_edi_bcs. Find 'AttDir' and 'SunPulseDir'
	% then remove them from varargin.
	iChar = find( [ cellfun(@ischar, varargin) ] );
	[tf_dir, iDir] = ismember({'AttDir', 'SunPulseDir'}, varargin(iChar));

	% Remove AttDir
	if tf_dir(1)
		idx     = iChar(iDir(1));
		att_dir = varargin{idx+1};
		varargin(idx:idx+1) = [];
	end
	% Remove SunPulseDir
	if tf_dir(2)
		idx          = iChar(iDir(2));
		sunpulse_dir = varargin{idx+1};
		varargin(idx:idx+1) = [];
	end

	% We must have either sunpulse or attitude data.
	%   - Use attitude by default
	%   - Assume the calibration files are stored with the data.
	if isempty(att_dir) && isempty(sunpulse_dir)
		att_dir = edi_dir;
	end
	
%------------------------------------%
% Get EDI Data                       %
%------------------------------------%
	% EDI
	[gd12_bcs, gd21_bcs] = mms_edi_bcs(sc, instr, mode, level, tstart, tend, varargin{:});

%------------------------------------%
% Rotate to SMPA                     %
%------------------------------------%
	% Get to SMPA
	if isempty(att_dir)
		% Attitude data
		attitude = [];
		warning('MMS_EDI_GSE:SMPA', 'No Attitude data. Cannot transform to SMPA.');
		
		% Cannot transform to SMPA, so copy variables
		fv_gd12_smpa  = gd12_bcs.fv_gd12_bcs;
		gun_gd12_smpa = gd12_bcs.gun_gd12_bcs;
		det_gd12_smpa = gd12_bcs.det_gd12_bcs;
		gun1_smpa     = gd12_bcs.gun1_bcs;

		fv_gd21_smpa  = gd21_bcs.fv_gd21_bcs;
		gun_gd21_smpa = gd21_bcs.gun_gd21_bcs;
		det_gd21_smpa = gd21_bcs.det_gd21_bcs;
		gun2_smpa     = gd21_bcs.gun2_bcs;
	else
		% Read attitude data
		[attitude, att_hdr] = mms_fdoa_read_defatt(sc, tstart, tend, att_dir);

		% Build matrix
		bcs2smpa = mms_fdoa_xbcs2smpa( att_hdr.('zMPA')(:,1)' );

		% Transform
		fv_gd12_smpa  = mrvector_rotate(bcs2smpa, gd12_bcs.fv_gd12_bcs);
		gun_gd12_smpa = mrvector_rotate(bcs2smpa, gd12_bcs.gun_gd12_bcs);
		det_gd12_smpa = mrvector_rotate(bcs2smpa, gd12_bcs.det_gd12_bcs);
	
		fv_gd21_smpa  = mrvector_rotate(bcs2smpa, gd12_bcs.fv_gd21_bcs);
		gun_gd21_smpa = mrvector_rotate(bcs2smpa, gd12_bcs.gun_gd21_bcs);
		det_gd21_smpa = mrvector_rotate(bcs2smpa, gd12_bcs.det_gd21_bcs);
	end

%------------------------------------%
% Despin Firing Vectors              %
%------------------------------------%

	%
	% Assume the principle axis of inertia (z-MPA)
	% is the same as the angular momentum vector (L)
	%

	% Despin using definitive attitude
	if isempty(sunpulse_dir)
		% Build matrix
		smpa2dmpa_gd12 = mms_fdoa_xdespin(attitude, t, 'L', 'EDI1_GUN');
		smpa2dmpa_gd21 = mms_fdoa_xdespin(attitude, t, 'L', 'EDI2_GUN');
	
	% Despin using sun pulse times.
	else
		% Read sun pulse data
		sunpulse = mms_dss_read_sunpulse(sc, tstart, tend, ...
		                                 'UniquePulse', true, ...
		                                 'Directory', sunpulse_dir);

		% Build matrix
		smpa2dmpa_gd12 = mms_dss_xdespin( sunpulse, gd12_bcs.epoch_gd12, 'EDI1_GUN' );
		smpa2dmpa_gd21 = mms_dss_xdespin( sunpulse, gd21_bcs.epoch_gd21, 'EDI2_GUN' );
	end

	% Transform
	fv_gd12_dmpa  = mrvector_rotate( smpa2dmpa_gd12, fv_gd12_smpa  );
	fv_gd21_dmpa  = mrvector_rotate( smpa2dmpa_gd21, fv_gd21_smpa  );


%------------------------------------%
% Spin Up Gun Positions              %
%------------------------------------%

	% Spin up using definitive attitude
	if isempty(sunpulse_dir)
		% Build matrix
		smpa2dmpa_gun1 = mms_fdoa_xdespin(attitude, t, 'L', 'EDI1_GUN');
		smpa2dmpa_det2 = mms_fdoa_xdespin(attitude, t, 'L', 'EDI1_DETECTOR');

		smpa2dmpa_gun2 = mms_fdoa_xdespin(attitude, t, 'L', 'EDI2_GUN');
		smpa2dmpa_det1 = mms_fdoa_xdespin(attitude, t, 'L', 'EDI2_DETECTOR');
	
	% Spin up using sunpulse
	else
		% Build rotation matrices
		smpa2dmpa_gun1 = mms_dss_xdespin( sunpulse, gd12_bcs.epoch_gd12, 'EDI1_GUN');
		smpa2dmpa_det2 = mms_dss_xdespin( sunpulse, gd12_bcs.epoch_gd12, 'EDI1_DETECTOR');
	
		smpa2dmpa_gun2 = mms_dss_xdespin( sunpulse, gd21_bcs.epoch_gd21, 'EDI2_GUN' );
		smpa2dmpa_det1 = mms_dss_xdespin( sunpulse, gd21_bcs.epoch_gd21, 'EDI2_DETECTOR' );
	end
	
	% Transform
	gun_gd12_dmpa = mrvector_rotate( smpa2dmpa_gun1, gun_gd12_smpa );
	det_gd12_dmpa = mrvector_rotate( smpa2dmpa_det2, det_gd12_smpa );
	gun1_dmpa     = mrvector_rotate( smpa2dmpa_gun1, gun1_smpa     );

	gun_gd21_dmpa = mrvector_rotate( smpa2dmpa_gun2, gun_gd21_smpa );
	det_gd21_dmpa = mrvector_rotate( smpa2dmpa_det1, det_gd21_smpa );
	gun2_dmpa     = mrvector_rotate( smpa2dmpa_gun2, gun2_smpa     );

%------------------------------------%
% Rotate to GSE                      %
%------------------------------------%
	% Can only rotate to GSE if we have attitude data.
	if ~isempty(attitude)
		% GEI -> Despun
		gei2despun_gd12 = mms_fdoa_xgei2despun(attitude, edi1_bcs.t_gd12, 'L');
		gei2despun_gd21 = mms_fdoa_xgei2despun(attitude, edi1_bcs.t_gd21, 'L');
	
		% Despun -> GEI
		despun2gei_gd12 = permute(gei2despun_gd12, [2, 1, 3]);
		despun2gei_gd21 = permute(gei2despun_gd21, [2, 1, 3]);
	
		%
		% Transform firing vectors and positions
		%
	
		% GD12
		fv_gd12_gei  = mrvector_rotate(despun2gei_gd12, fv_gd12_dsl);
		gun_gd12_gei = mrvector_rotate(despun2gei_gd12, gun_gd12_dsl);
		det_gd12_gei = mrvector_rotate(despun2gei_gd12, det_gd12_dsl);
	
		% GD21
		fv_gd21_gei  = mrvector_rotate(despun2gei_gd21, fv_gd21_dsl);
		gun_gd21_gei = mrvector_rotate(despun2gei_gd21, gun_gd21_dsl);
		det_gd21_gei = mrvector_rotate(despun2gei_gd21, det_gd21_dsl);

		%
		% Transformation matrix GEI -> GSE
		%   - Modified Julian Date (mjd).
		%   - UTC seconds since midnight (ssm).
		%
	
		% GD12
		timevec = MrCDF_Epoch_Breakdown( gd12_bcs.t_gd12 );
		mjd     = date2mjd(timevec(1,:), timevec(2,:), timevec(3,:));
		ssm = timevec(4,:) * 3600.0 + ...
			  timevec(5,:) * 60.0   + ...
			  timevec(6,:)          + ...
			  timevec(7,:) * 1e-3   + ...
			  timevec(8,:) * 1e-6   + ...
			  timevec(9,:) * 1e-9;
		GEI2GSE_gd12 = gei2gse(mjd, ssm);
	
		% GD21
		timevec = MrCDF_Epoch_Breakdown( gd21_bcs.t_gd21 );
		mjd     = date2mjd(timevec(1,:), timevec(2,:), timevec(3,:));
		ssm = timevec(4,:) * 3600.0 + ...
			  timevec(5,:) * 60.0   + ...
			  timevec(6,:)          + ...
			  timevec(7,:) * 1e-3   + ...
			  timevec(8,:) * 1e-6   + ...
			  timevec(9,:) * 1e-9;
		GEI2GSE_gd21 = gei2gse(mjd, ssm);

		%
		% Transform to GSE
		%
	
		% GD12
		fv_gd12_gse  = mrvector_rotate(GEI2GSE_gd12, fv_gd12_gei);
		gun_gd12_gse = mrvector_rotate(GEI2GSE_gd12, gun_gd12_gei);
		det_gd12_gse = mrvector_rotate(GEI2GSE_gd12, det_gd12_gei);
	
		% GD21
		fv_gd21_gse  = mrvector_rotate(GEI2GSE_gd21, fv_gd21_gei);
		gun_gd21_gse = mrvector_rotate(GEI2GSE_gd21, gun_gd21_gei);
		det_gd21_gse = mrvector_rotate(GEI2GSE_gd21, det_gd21_gei);
	end

%------------------------------------%
% Return                             %
%------------------------------------%
	
	% EDI1 output structure
	%   - EDI1 contains gun1 and detector2
	gd12_dmpa = struct( 'epoch_gd12',     gd12_bcs.epoch_gd12,   ...
	                    'gun_gd12_dmpa',  gun_gd12_dmpa, ...
	                    'det_gd12_dmpa',  det_gd12_dmpa, ...
	                    'gun1_dmpa',      gun1_dmpa,     ...
	                    'fv_gd12_dmpa',   fv_gd12_dmpa,  ...
	                    'q_gd12',         gd12_bcs.q_gd12,  ...
	                    'tof_gd12',       gd12_bcs.tof_gd12 );
	
	% EDI2 output structure
	%   - EDI2 contains gun1 and detector2
	gd21_dmpa = struct( 'epoch_gd21',     gd21_bcs.epoch_gd21,   ...
	                    'gun_gd21_dmpa',  gun_gd21_dmpa, ...
	                    'det_gd21_dmpa',  det_gd21_dmpa, ...
	                    'gun2_dmpa',      gun2_dmpa,     ...
	                    'fv_gd21_dmpa',   fv_gd21_dmpa,  ...
	                    'q_gd21',         gd21_bcs.q_gd21,  ...
	                    'tof_gd21',       gd21_bcs.tof_gd21 );
end