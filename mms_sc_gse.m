%
% Name
%   mms_sc_gse
%
% Purpose
%   Take l1a search coil data, calibrate, despin, and rotate into
%   GSE.
%
%   Process:
%     1. Find, read, calibrate l1a data
%     2. Despin
%     3. Transform to GSE
%
% Calling Sequence
%   [T, B_GSE] = mms_sc_bcs(SC, INSTR, MODE, TSTART, TEND)
%     Find and read search coil magnetometer and transfer function
%     data from MMS scacecraft SC (e.g. 'mms1') from instrument
%     INSTR (i.e. 'scm') in telemetry mode MODE (e.g. 'comm') during
%     the time interval [TSTART, TEND]. Calibrate the data and
%     transform into the spacecraft body frame (BCS), despin, and
%     transform into GSE. Return the data and its time tags B_GSE
%     and T.
%
%   [..., B_DSL] = mms_sc_bcs(SC, INSTR, MODE, TSTART, TEND)
%     Also return the magnetic field in DSL.
%
%   [..., B_SMPA] = mms_sc_bcs(SC, INSTR, MODE, TSTART, TEND)
%     Also return the magnetic field in SMPA.
%
%   [..., B_BCS] = mms_sc_bcs(SC, INSTR, MODE, TSTART, TEND)
%     Also return the magnetic field in BCS.
%
%   [__] = mms_sc_bcs(__, 'ParamName', ParamValue)
%     Any parameter name-value pair found below.
%
% Parameters
%   SC              in, required, type = char
%   INSTR           in, required, type = char
%   MODE            in, required, type = char
%   TSTART          in, optional, type = char
%   TEND            in, optional, type = char
%   'AttDir'        in, optional, type = char default = 'DataDir'
%                   Directory in which to find attitude data.
%   'CalDir'        in, optional, type = char default = 'DataDir'
%                   Directory in which to find calibration data.
%   'DataDir'       in, optional, type = char default = 'DataDir'
%                   Directory in which to find search coil data.
%   'Duration'      in, optional, type = char default = 'DataDir'
%                   Time duration over which to calibrate data.
%
% Returns
%   T               out, required, type=1xN int64
%   B_GSE           out, required, type=3xN double
%   B_DSL           out, optional, type=3xN double
%   B_SMPA          out, optional, type=3xN double
%   B_BCS           out, optional, type=3xN double
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-15      Written by Matthew Argall
%
function [t, b_gse, b_dsl, b_smpa, b_bcs] = mms_sc_gse(sc, instr, mode, tstart, tend, varargin)

%------------------------------------%
% Inputs                             %
%------------------------------------%
	% Check if the attitude directory was given.
	%   - All others are passed to mms_sc_bcs.m
	%   - ismember requires all strings
	iChar = find( [ cellfun(@ischar, varargin) ] );
	[tf_att, iAtt] = ismember('AttDir', varargin(iChar));
	if tf_att
		idx                 = iChar(iAtt);
		att_dir             = varargin{idx+1};
		varargin(idx:idx+1) = [];
	else
		att_dir = pwd();
	end

%------------------------------------%
% Calibrated Mag in BCS              %
%------------------------------------%
	[b_bcs, t] = mms_sc_bcs(sc, instr, mode, tstart, tend, varargin{:});

%------------------------------------%
% Rotate to SMPA                     %
%------------------------------------%
	% Read attitude data
	[attitude, att_hdr] = mms_fdoa_read_defatt(sc, tstart, tend, att_dir);

	% Build matrix
	bcs2smpa = mms_fdoa_xbcs2smpa( att_hdr.('zMPA')(:,1)' );

	% Transform
	b_smpa = mrvector_rotate(bcs2smpa, b_bcs);

%------------------------------------%
% Despin                             %
%------------------------------------%

	%
	% Assume the principle axis of inertia (z-MPA)
	% is the same as the angular momentum vector (L)
	%
	smpa2dsl = mms_fdoa_xdespin(attitude, t, 'L', [ instr '_123' ]);

	% Rotate to DSL
	b_dsl = mrvector_rotate(smpa2dsl, b_smpa);

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
	timevec = MrCDF_Epoch_Breakdown(t)';
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