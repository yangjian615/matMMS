%
% Name
%   mms_edi_read_efield
%
% Purpose
%   Read MMS EDI electric field mode data files.
%
% Calling Sequence
%   [EDI1_BCS, EDI2_BCS] = mms_edi_bcs(SC, INSTR, MODE, LEVEL, TSTART, TEND, EDI_DIR)
%     Read EDI electric field mode data captured by spacecraft SC
%     (e.g. 'mms3'), instrument INSTR (e.g. 'edi'), from telemetry
%     mode MODE and data product level LEVEL between the time
%     interval of [TSTART, TEND]. Data can be found in directory
%     EDI_DIR. Times should be provided in ISO format:
%     'yyyy-mm-ddThh:mm_ss'. Firing angles are converted to firing
%     vectors in BCS and returned in the output structures.
%
% Parameters
%   SC              in, required, type=char/cell
%   INSTR           in, required, type=char
%   MODE            in, required, type=char
%   LEVEL           in, required, type=char
%   TSTART          in, required, type=char
%   TEND            in, required, type=char
%   'DataDir'       in, optional, type=char, default=pwd()
%                   Directory in which to find EDI data.
%   'Quality'       in, optional, type=integer, default=[]
%                   Quality of beams to select. Can be a scalar or array with
%                     values of 0, 1, 2, or 3. The default is to select all
%                     beams.
%
% Returns
%   EDI1_BCS        out, required, type=structure
%                   Fields are:
%                     't_gd12'        -  TT2000 Epoch time for gun 1 and detector 2.
%                     'gun_gd12_bcs'  -  Gun1 position in BCS.
%                     'det_gd12_bcs'  -  Detector2 position in BCS.
%                     'fv_gd12_bcs'   -  Firing vectors from gun1 in BCS.
%                     'q_gd12'        -  Quality flag.
%                     'tof_gd12'      -  Time of flight.
%   EDI2_BCS        out, required, type=structure
%                   Fields are:
%                     't_gd21'        -  TT2000 Epoch time for gun 2 and detector 1.
%                     'gun_gd21_bcs'  -  Gun2 position in BCS.
%                     'det_gd21_bcs'  -  Detector1 position in BCS.
%                     'fv_gd21_bcs'   -  Firing vectors from gun2 in BCS.
%                     'q_gd21'        -  Quality flag.
%                     'tof_gd21'      -  Time of flight.
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-18      Written by Matthew Argall
%   2015-04-23      Assigned time from gd12 to output gd21 structure. Fixed. - MRA
%   2015-04-24      Renamed EDI_DIR to DATADIR and made a name-value pair.
%                     Added "Quality" parameter. - MRA
%
function [gd12, gd21] = mms_edi_bcs(sc, instr, mode, level, tstart, tend, varargin)

	quality = [];
	edi_dir = pwd();
	
	nOptArg = length(varargin);
	for ii = 1 : 2 : nOptArg
		switch varargin{ii}
			case 'DataDir'
				edi_dir = varargin{ii+1};
			case 'Quality'
				quality = varargin{ii+1};
			otherwise
				error( ['Optional parameter not recognized: "' varargin{ii+1} '".'] );
		end
	end
	
%------------------------------------%
% Get EDI Data                       %
%------------------------------------%
	% Gun positions and detectors in BCS
	gun_gd12_bcs = mms_instr_origins_ocs('EDI1_GUN');
	det_gd12_bcs = mms_instr_origins_ocs('EDI1_DETECTOR');
	gun_gd21_bcs = mms_instr_origins_ocs('EDI2_GUN');
	det_gd21_bcs = mms_instr_origins_ocs('EDI2_DETECTOR');
	
	% Read EDI efield data
	[gd12, gd21] = mms_edi_read_efield(sc, instr, mode, level, tstart, tend, edi_dir);

	% Convert to radians
	deg2rad      = pi / 180.0;
	polar_gd12   = gd12.polar_gd12   * deg2rad;
	azimuth_gd12 = gd12.azimuth_gd12 * deg2rad;
	polar_gd21   = gd21.polar_gd21   * deg2rad;
	azimuth_gd21 = gd21.azimuth_gd21 * deg2rad;

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
% Filter by Quality                  %
%------------------------------------%
	if ~isempty(quality)
		if length(quality) > 4 || max(quality) > 3 || min(quality) < 0
			error('Quality can have only these unique values: 0, 1, 2, 3.');
		end
		
		% Find quality
		iq_gd12 = ismember(gd12.q_gd12, quality);
		iq_gd21 = ismember(gd21.q_gd21, quality);
		
		% Select data
		fv_gd12         = fv_gd12(:, iq_gd12);
		gd12.epoch_gd12 = gd12.epoch_gd12(iq_gd12);
		gd12.q_gd12     = gd12.q_gd12(iq_gd12);
		gd12.tof_gd12   = gd12.tof_gd12(iq_gd12);
		
		fv_gd21         = fv_gd21(:, iq_gd21);
		gd21.epoch_gd21 = gd21.epoch_gd21(iq_gd21);
		gd21.q_gd21     = gd21.q_gd21(iq_gd21);
		gd21.tof_gd21   = gd21.tof_gd21(iq_gd21);
	end
keyboard
%------------------------------------%
% Transform to BCS                   %
%------------------------------------%
	
	% Transformations from GUN1 and GUN2 to BCS
	edi12bcs = mms_instr_xxyz2ocs('EDI1_GUN');
	edi22bcs = mms_instr_xxyz2ocs('EDI2_GUN');

	% Transform firing vectors
	fv_gd12_bcs = mrvector_rotate( edi12bcs, fv_gd12 );
	fv_gd21_bcs = mrvector_rotate( edi22bcs, fv_gd21 );
	
%------------------------------------%
% Output                             %
%------------------------------------%
	%
	% I would like to create a new structure GD12_BCS and copy fields from
	% GD12 into it. I do not see a simple way of doing this in MATLAB, so
	% resort to appending new fields to the old structure. This means having
	% to do away with old fields.
	%

	% Remove fields from the structures
	gd12 = rmfield(gd12, {'polar_gd12', 'azimuth_gd12'});
	gd21 = rmfield(gd21, {'polar_gd21', 'azimuth_gd21'});
	
	% EDI1 output structure
	gd12.('gun_gd12_bcs') = gun_gd12_bcs;
	gd12.('det_gd12_bcs') = det_gd12_bcs;
	gd12.('gun1_bcs')     = gun_gd12_bcs - det_gd21_bcs;
	gd12.('fv_gd12_bcs')  = fv_gd12_bcs;
	
	% EDI2 output structure
	gd21.('gun_gd21_bcs') = gun_gd21_bcs;
	gd21.('det_gd21_bcs') = det_gd21_bcs;
	gd21.('gun2_bcs')     = gun_gd21_bcs - det_gd12_bcs;
	gd21.('fv_gd21_bcs')  = fv_gd21_bcs;
end