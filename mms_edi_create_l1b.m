%
% Name
%   mms_edi_create_l1b
%
% Purpose
%   Turn level 1A EDI electric-field mode data into level 1B quality data.
%   L1B implies calibrated data in the spinning spacecraft body coordinate
%   system.
%
% Calling Sequence
%   EDI = mms_edi_create_l1b(FILENAME, TSTART, TEND)
%     Read EDI electric field mode data from files named FILENAMES
%     between the time interval of [TSTART, TEND]. Return a data
%     structure EDI with fields given below. Times should be provided
%     in ISO format: 'yyyy-mm-ddThh:mm_ss'.
%
%   [GD12_GSE, GD21_GSE] = mms_edi_create_l2(__, 'ParamName', ParamValue)
%     Include any of the parameter name-value pairs below.
%
% Parameters
%   FILENAMES       in, required, type=char/cell
%   TSTART          in, required, type=char
%   TEND            in, required, type=char
%   'CS_123'        in, optional, type=boolean, default=false
%                   If true, data in the EDI sensor 123 coordinate system is included
%                     in the returned data structure.
%   'CS_BCS'        in, optional, type=boolean, default=false
%                   If true, data in the BCS is included in the returned data structure.
%   'Quality'       in, optional, type=integer/array, default=[]
%                   The quality of beams to process. Options are 0, 1, 2, 3, or any
%                     combination of the four. By default, beams of all qualities are
%                     processed.
%
% Returns
%   EDI             out, required, type=structure
%                   In addition to the fields return by mms_edi_read_l1a_efield, we have:
%                     'gun_gd12_*'     -  Gun1 position in BCS.
%                     'det_gd12_*'     -  Detector2 position in BCS.
%                     'virtual_gun1_*' -  Position of gun1 on virtual spacecraft in BCS.
%                     'fv_gd12_*'      -  Firing vectors from gun1 in BCS.
%
%                     'gun_gd21_*'     -  Gun2 position in BCS.
%                     'det_gd21_*'     -  Detector1 position in BCS.
%                     'virtual_gun2_*' -  Position of gun2 on virtual spacecraft in BCS.
%                     'fv_gd21_*'      -  Firing vectors from gun1 in BCS.
%                   Where "*" indicates the coordinate system:
%                     '123', 'bcs'.
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-18      Written by Matthew Argall
%   2015-04-23      Assigned time from gd12 to output gd21 structure. Fixed. - MRA
%   2015-04-24      Renamed EDI_DIR to DATADIR and made a name-value pair.
%                     Added "Quality" parameter. - MRA
%   2015-04-29      Throw error if quality filter returns zero beams. - MRA
%   2015-04-29      Remove quality selection and firing vector calculations.
%                     Accept filenames as inputs. Select coordinate system. - MRA
%   2015-05-21      Renamed from mms_edi_bcs to mms_edi_create_l1b. - MRA
%   2015-08-11      Store gun and detector positions as 3xN instead of Nx3 arrays. - MRA
%
function edi = mms_edi_create_l1b(filenames, tstart, tend, varargin)

	% Defaults
	quality = [];
	cs_123  = false;
	cs_bcs  = true;
	
	% Optional parameters
	nOptArg = length(varargin);
	for ii = 1 : 2 : nOptArg
		switch varargin{ii}
			case 'CS_123'
				cs_123 = varargin{ii+1};
			case 'CS_BCS'
				cs_bcs = varargin{ii+1};
			case 'Quality'
				quality = varargin{ii+1};
			otherwise
				error( ['Optional parameter not recognized: "' varargin{ii+1} '".'] );
		end
	end
	
%------------------------------------%
% Get EDI Data                       %
%------------------------------------%
	% Guns and detector positions in BCS
	gun_gd12_bcs = mms_instr_origins_ocs('EDI1_GUN')';
	det_gd12_bcs = mms_instr_origins_ocs('EDI1_DETECTOR')';
	gun_gd21_bcs = mms_instr_origins_ocs('EDI2_GUN')';
	det_gd21_bcs = mms_instr_origins_ocs('EDI2_DETECTOR')';
	
	% Read EDI efield data
	edi = mms_edi_read_l1a_efield(filenames, tstart, tend, 'Quality', quality);

%------------------------------------%
% Transform to BCS                   %
%------------------------------------%
	
	% Transformations from GUN1 and GUN2 to BCS
	edi12bcs = mms_instr_xxyz2ocs('EDI1_GUN');
	edi22bcs = mms_instr_xxyz2ocs('EDI2_GUN');

	% Transform firing vectors
	fv_gd12_bcs = mrvector_rotate( edi12bcs, edi.fv_gd12_123 );
	fv_gd21_bcs = mrvector_rotate( edi22bcs, edi.fv_gd21_123 );
	
%------------------------------------%
% Output                             %
%------------------------------------%

	% Remove data in 123 system?
	if ~cs_123
		edi = rmfield( edi, {'fv_gd12_123', 'fv_gd21_123'} );
	end
	
	% Add data in BCS
	if cs_bcs
		edi.('gun_gd12_bcs')     = gun_gd12_bcs;
		edi.('det_gd12_bcs')     = det_gd12_bcs;
		edi.('virtual_gun1_bcs') = gun_gd12_bcs - det_gd21_bcs;
		edi.('fv_gd12_bcs')      = fv_gd12_bcs;
	
		edi.('gun_gd21_bcs')     = gun_gd21_bcs;
		edi.('det_gd21_bcs')     = det_gd21_bcs;
		edi.('virtual_gun2_bcs') = gun_gd21_bcs - det_gd12_bcs;
		edi.('fv_gd21_bcs')      = fv_gd21_bcs;
	end
end