%
% Name
%   mms_edi_create_l2pre_efield
%
% Purpose
%   Turn level 1A EDI electric-field mode data into level 2-pre quality data
%   for use with BESTARG or the beam convergence method for calculating
%   the electric field. L2 implies despun data in the spacecraft reference
%   frame (no VxB removal) and in GSE coordinates. The "pre" indicates that
%   we have not used the beam information to calculate the electric field yet.
%
% Calling Sequence
%   EDI = mms_edi_create_l2pre_efield(FILENAMES, TSTART, TEND)
%     Read EDI electric field mode data from files named FILENAMES
%     between the time interval of [TSTART, TEND]. Return a data
%     structure EDI with fields given below. Times should be provided
%     in ISO format: 'yyyy-mm-ddThh:mm_ss'.
%
%   [GD12_GSE, GD21_GSE] = mms_edi_create_l2pre_efield(__, 'ParamName', ParamValue)
%     Include any of the parameter name-value pairs below.
%
% Parameters
%   FILENAMES       in, required, type=char/cell
%   TSTART          in, required, type=char
%   TEND            in, required, type=char
%   'Attitude'      in, optional, type = struct, default = []
%                   Structure of definitive attitude data returned by mms_fdoa_read_defatt.m
%   'CS_123'        in, optional, type=boolean, default=false
%                   If true, data in the EDI sensor 123 coordinate system is included
%                     in the returned data structure.
%   'CS_BCS'        in, optional, type=boolean, default=false
%                   If true, data in the BCS is included in the returned data structure.
%   'CS_SMPA'       in, optional, type=boolean, default=false
%                   If true, data in the SMPA coordinate system is included in the
%                     returned data structure.
%   'CS_DMPA'       in, optional, type=boolean, default=false
%                   If true, data in the DMPA coordinate system is included in the
%                     returned data structure.
%   'CS_GSE'        in, optional, type=boolean, default=false
%                   If true, data in the GSE coordinate system is included in the
%                     returned data structure.
%   'Quality'       in, optional, type=integer/array, default=[]
%                   The quality of beams to process. Options are 0, 1, 2, 3, or any
%                     combination of the four. By default, beams of all qualities are
%                     processed.
%   'SunPulse'      in, optional, type=struct, default=[]
%                   Structure of HK 101 sunpulse data returned by mms_dss_read_sunpulse.m
%   'zMPA'          in, optional, type=boolean, default=false
%                   The z-MPA axis used for rotating from BCS to SMPA coordinates.
%                     Available in the returned header structure from mms_fdoa_read_defatt.m
%
% Returns
%   EDI             out, required, type=structure
%                   In addition to fields returned by mms_edi_bcs, we have:
%                     'gun_gd12_*'     -  Gun1 position.
%                     'det_gd12_*'     -  Detector2 position.
%                     'virtual_gun1_*' -  Gun1 position on virtual spacecraft.
%                     'fv_gd12_*'      -  Firing vectors from gun1.
%
%                     'gun_gd21_*'     -  Gun2 position.
%                     'det_gd21_*'     -  Detector1 position.
%                     'virtual_gun2_*' -  Gun2 position on virtual spacecraft.
%                     'fv_gd21_*'      -  Firing vectors from gun1.
%                   Where "*" indicates the coordinate system:
%                     '123', 'bcs', 'smpa', 'dmpa', or 'gse'.
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-19  Written by Matthew Argall
%   2015-04-20  Data is despun with attitude data unless SunPulseDir is given. - MRA
%   2015-05-14  Added detail to header. rlm
%   2015-06-20  Renamed from mms_edi_bcs to mms_edi_create_l2. - MRA
%   2015-07-20  Renamed from mms_edi_create_l2 to mms_edi_create_l2pre_efield. - MRA
%
function edi = mms_edi_create_l2(filenames, tstart, tend, varargin)

	% Defaults
	sunpulse = [];
	attitude = [];
	z_mpa    = [];
	cs_123   = false;
	cs_bcs   = false;
	cs_smpa  = false;
	cs_dmpa  = false;
	cs_gse   = true;
	quality  = [];
	
	% Optional parameters
	nOptArg = length(varargin);
	for ii = 1 : 2 : nOptArg
		switch varargin{ii}
			case 'Attitude'
				attitude = varargin{ii+1};
			case 'Sunpulse'
				sunpulse = varargin{ii+1};
			case 'CS_123'
				cs_123 = varargin{ii+1};
			case 'CS_BCS'
				cs_bcs = varargin{ii+1};
			case 'CS_SMPA'
				cs_smpa = varargin{ii+1};
			case 'CS_DMPA'
				cs_dmpa = varargin{ii+1};
			case 'CS_GSE'
				cs_gse = varargin{ii+1};
			case 'Quality'
				quality = varargin{ii+1};
			case 'zMPA'
				z_mpa = varargin{ii+1};
			otherwise
				error( ['Optional parameter not recognized: "' varargin{ii+1} '".'] );
		end
	end
	
%------------------------------------%
% Get EDI Data                       %
%------------------------------------%
	% EDI
	edi = mms_edi_create_l1b(filenames, tstart, tend, ...
	                         'CS_123',  cs_123, ...
	                         'CS_BCS',  true, ...
	                         'Quality', quality);

%------------------------------------%
% Rotate to SMPA                     %
%------------------------------------%
	% Get to SMPA
	if isempty(z_mpa)
		% Attitude data
		attitude = [];
		warning('MMS_EDI_GSE:SMPA', 'zMPA not given. Cannot transform to SMPA.');
		cs_smpa = false;
		
		% Cannot transform to SMPA, so copy variables
		fv_gd12_smpa      = edi.fv_gd12_bcs;
		gun_gd12_smpa     = edi.gun_gd12_bcs;
		det_gd12_smpa     = edi.det_gd12_bcs;
		virtual_gun1_smpa = edi.virtual_gun1_bcs;

		fv_gd21_smpa      = edi.fv_gd21_bcs;
		gun_gd21_smpa     = edi.gun_gd21_bcs;
		det_gd21_smpa     = edi.det_gd21_bcs;
		virtual_gun2_smpa = edi.virtual_gun2_bcs;
	else
		% Build matrix
		bcs2smpa = mms_fdoa_xbcs2smpa( z_mpa );

		% Transform
		fv_gd12_smpa      = mrvector_rotate(bcs2smpa, edi.fv_gd12_bcs);
		gun_gd12_smpa     = mrvector_rotate(bcs2smpa, edi.gun_gd12_bcs);
		det_gd12_smpa     = mrvector_rotate(bcs2smpa, edi.det_gd12_bcs);
		virtual_gun1_smpa = mrvector_rotate(bcs2smpa, edi.virtual_gun1_bcs);
	
		fv_gd21_smpa      = mrvector_rotate(bcs2smpa, edi.fv_gd21_bcs);
		gun_gd21_smpa     = mrvector_rotate(bcs2smpa, edi.gun_gd21_bcs);
		det_gd21_smpa     = mrvector_rotate(bcs2smpa, edi.det_gd21_bcs);
		virtual_gun2_smpa = mrvector_rotate(bcs2smpa, edi.virtual_gun2_bcs);
	end

%------------------------------------%
% Despin Firing Vectors              %
%------------------------------------%

	%
	% Assume the principle axis of inertia (z-MPA)
	% is the same as the angular momentum vector (L)
	%

	% Despin using definitive attitude
	if isempty(sunpulse)
		smpa2dmpa_gd12 = mms_fdoa_xdespin( attitude, edi.tt2000_gd12, 'L' );
		smpa2dmpa_gd21 = mms_fdoa_xdespin( attitude, edi.tt2000_gd21, 'L' );
	
	% Despin using sun pulse times.
	else
		smpa2dmpa_gd12 = mms_dss_xdespin( sunpulse, edi.tt2000_gd12 );
		smpa2dmpa_gd21 = mms_dss_xdespin( sunpulse, edi.tt2000_gd21 );
	end

	% Despin
	fv_gd12_dmpa      = mrvector_rotate( smpa2dmpa_gd12, fv_gd12_smpa );
	gun_gd12_dmpa     = mrvector_rotate( smpa2dmpa_gd12, gun_gd12_smpa );
	det_gd12_dmpa     = mrvector_rotate( smpa2dmpa_gd12, det_gd12_smpa );
	virtual_gun1_dmpa = mrvector_rotate( smpa2dmpa_gd12, virtual_gun1_smpa );
	
	fv_gd21_dmpa      = mrvector_rotate( smpa2dmpa_gd21, fv_gd21_smpa );
	gun_gd21_dmpa     = mrvector_rotate( smpa2dmpa_gd21, gun_gd21_smpa );
	det_gd21_dmpa     = mrvector_rotate( smpa2dmpa_gd21, det_gd21_smpa );
	virtual_gun2_dmpa = mrvector_rotate( smpa2dmpa_gd21, virtual_gun2_smpa );

%------------------------------------%
% Rotate to GSE                      %
%------------------------------------%
	% Can only rotate to GSE if we have attitude data.
	if ~isempty(attitude) && cs_gse
		% GEI -> Despun
		gei2despun_gd12 = mms_fdoa_xgei2despun(attitude, edi.tt2000_gd12, 'L');
		gei2despun_gd21 = mms_fdoa_xgei2despun(attitude, edi.tt2000_gd21, 'L');
	
		% Despun -> GEI
		despun2gei_gd12 = permute(gei2despun_gd12, [2, 1, 3]);
		despun2gei_gd21 = permute(gei2despun_gd21, [2, 1, 3]);
	
		%
		% Transform firing vectors and positions
		%
	
		% GD12
		fv_gd12_gei      = mrvector_rotate(despun2gei_gd12, fv_gd12_dmpa);
		gun_gd12_gei     = mrvector_rotate(despun2gei_gd12, gun_gd12_dmpa);
		det_gd12_gei     = mrvector_rotate(despun2gei_gd12, det_gd12_dmpa);
		virtual_gun1_gei = mrvector_rotate(despun2gei_gd12, virtual_gun1_dmpa);
	
		% GD21
		fv_gd21_gei      = mrvector_rotate(despun2gei_gd21, fv_gd21_dmpa);
		gun_gd21_gei     = mrvector_rotate(despun2gei_gd21, gun_gd21_dmpa);
		det_gd21_gei     = mrvector_rotate(despun2gei_gd21, det_gd21_dmpa);
		virtual_gun2_gei = mrvector_rotate(despun2gei_gd21, virtual_gun2_dmpa);

		%
		% Transformation matrix GEI -> GSE
		%   - Modified Julian Date (mjd).
		%   - UTC seconds since midnight (ssm).
		%
	
		% GD12
		timevec = MrCDF_Epoch_Breakdown( gd12_bcs.tt2000_gd12 );
		mjd     = date2mjd(timevec(1,:), timevec(2,:), timevec(3,:));
		ssm = timevec(4,:) * 3600.0 + ...
			  timevec(5,:) * 60.0   + ...
			  timevec(6,:)          + ...
			  timevec(7,:) * 1e-3   + ...
			  timevec(8,:) * 1e-6   + ...
			  timevec(9,:) * 1e-9;
		GEI2GSE_gd12 = gei2gse(mjd, ssm);
	
		% GD21
		timevec = MrCDF_Epoch_Breakdown( gd21_bcs.tt2000_gd21 );
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
		fv_gd12_gse      = mrvector_rotate(GEI2GSE_gd12, fv_gd12_gei);
		gun_gd12_gse     = mrvector_rotate(GEI2GSE_gd12, gun_gd12_gei);
		det_gd12_gse     = mrvector_rotate(GEI2GSE_gd12, det_gd12_gei);
		virtual_gun1_gse = mrvector_rotate(GEI2GSE_gd12, virtual_gun1_gei);
	
		% GD21
		fv_gd21_gse      = mrvector_rotate(GEI2GSE_gd21, fv_gd21_gei);
		gun_gd21_gse     = mrvector_rotate(GEI2GSE_gd21, gun_gd21_gei);
		det_gd21_gse     = mrvector_rotate(GEI2GSE_gd21, det_gd21_gei);
		virtual_gun2_gse = mrvector_rotate(GEI2GSE_gd21, virtual_gun2_gei);
	else
		warning('MMS_EDI:GSE', 'No attitude data. Cannot rotate to GSE.')
		cs_gse = false;
	end

%------------------------------------%
% Return                             %
%------------------------------------%
	
	%
	% Add data to structure
	% 
	
	if ~cs_bcs
		edi = rmfield( edi, { 'gun_gd12_bcs', ...
		                      'det_gd12_bcs', ...
		                      'virtual_gun1_bcs', ...
		                      'fv_gd12_bcs', ...
		                      'gun_gd21_bcs', ...
		                      'det_gd21_bcs', ...
		                      'virtual_gun2_bcs', ...
		                      'fv_gd21_bcs' } );
	end
	
	if cs_smpa
		edi.( 'gun_gd12_smpa' )     = gun_gd12_smpa;
		edi.( 'det_gd12_smpa' )     = det_gd12_smpa;
		edi.( 'virtual_gun1_smpa' ) = virtual_gun1_smpa;
		edi.( 'fv_gd12_smpa' )      = fv_gd12_smpa;
	
		edi.( 'gun_gd21_smpa' )     = gun_gd21_smpa;
		edi.( 'det_gd21_smpa' )     = det_gd21_smpa;
		edi.( 'virtual_gun2_smpa' ) = virtual_gun2_smpa;
		edi.( 'fv_gd21_smpa' )      = fv_gd21_smpa;
	end
	
	if cs_dmpa
		edi.( 'gun_gd12_dmpa' )     = gun_gd12_dmpa;
		edi.( 'det_gd12_dmpa' )     = det_gd12_dmpa;
		edi.( 'virtual_gun1_dmpa' ) = virtual_gun1_dmpa;
		edi.( 'fv_gd12_dmpa' )      = fv_gd12_dmpa;
	
		edi.( 'gun_gd21_dmpa' )     = gun_gd21_dmpa;
		edi.( 'det_gd21_dmpa' )     = det_gd21_dmpa;
		edi.( 'virtual_gun2_dmpa' ) = virtual_gun2_dmpa;
		edi.( 'fv_gd21_dmpa' )      = fv_gd21_dmpa;
	end
	
	if cs_gse
		edi.( 'gun_gd12_gse' )     = gun_gd12_gse;
		edi.( 'det_gd12_gse' )     = det_gd12_gse;
		edi.( 'virtual_gun1_gse' ) = virtual_gun1_gse;
		edi.( 'fv_gd12_gse' )      = fv_gd12_gse;
	
		edi.( 'gun_gd21_gse' )     = gun_gd21_gse;
		edi.( 'det_gd21_gse' )     = det_gd21_gse;
		edi.( 'virtual_gun2_gse' ) = virtual_gun2_gse;
		edi.( 'fv_gd21_gse' )      = fv_gd21_gse;
	end
end