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
%   [..., B_DSL] = mms_fg_bcs(SC, INSTR, MODE, TSTART, TEND)
%     Also return the magnetic field in DSL.
%
%   [..., B_SMPA] = mms_fg_bcs(SC, INSTR, MODE, TSTART, TEND)
%     Also return the magnetic field in SMPA.
%
%   [..., B_BCS] = mms_fg_bcs(SC, INSTR, MODE, TSTART, TEND)
%     Also return the magnetic field in BCS.
%
%   [__] = mms_fg_bcs(__, 'ParamName', ParamValue)
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
%                   Directory in which to find fluxgate data.
%   'SunPulseDir'   in, optional, type=char, default='DataDir'
%                   Directory in which to find sun pulse data. If provided,
%                     sun pulse times will be used to depsin data. Otherwise,
%                     spin phase data from the devinitive attitude files is
%                     used.
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
%   2015-04-20      Added 'SunPulseDir' parameter.
%
function [t, b_gse, b_dmpa, b_smpa, b_bcs, b_omb] = mms_fg_gse(sc, instr, mode, tstart, tend, varargin)

%------------------------------------%
% Inputs                             %
%------------------------------------%
	% Check if the attitude directory was given.
	%   - All others are passed to mms_fg_bcs.m
	[tf_dir, iDir] = ismember({'SunPulseDir', 'AttDir'}, varargin);
	
	% Sun Pulse directory
	if tf_dir(1)
		sunpulse_dir = varargin{iDir(1) + 1};
		varargin(iDir(1):iDir(1)+1) = [];
	else
		sunpulse_dir = '';
	end
	
	% Attitude directory
	if tf_dir(2)
		att_dir = varargin{iDir(2) + 1};
		varargin(iDir(2):iDir(2)+1) = [];
	else
		att_dir = '';
	end
	
	% Must have attitude or sunpulse to at least despin
	if isempty(sunpulse_dir) && isempty(att_dir)
		attitude = pwd();
	end

%------------------------------------%
% Calibrated Mag in BCS              %
%------------------------------------%
	[b_bcs, t, b_smpa, b_omb] = mms_fg_bcs(sc, instr, mode, tstart, tend, varargin{:});

%------------------------------------%
% Despin                             %
%------------------------------------%
	%
	% Assume the principle axis of inertia (z-MPA)
	% is the same as the angular momentum vector (L)
	%
	
	% Read attitude data?
	%   - We need it to despin and to rotate into GSE.
	%   - If attitude data not available, cannot rotate into GSE
	if ~isempty(att_dir)
		attitude = mms_fdoa_read_defatt(sc, tstart, tend, att_dir);
	else
		attitude = [];
	end
	
	% Despin using attitude data
	if isempty(sunpulse_dir)
		% Build matrix
		smpa2dmpa = mms_fdoa_xdespin( attitude, t, 'L', [upper(instr) '_123'] );
	
	% Despin using sunpulse directory
	else
		% Read sun pulse data
		sunpulse = mms_dss_read_sunpulse(sc, tstart, tend, ...
		                                 'UniquePulse', true, ...
		                                 'Directory', sunpulse_dir);

		% Build matrix
		smpa2dmpa = mms_dss_xdespin( sunpulse, t );
	end

	% Rotate to DSL
	b_dmpa = mrvector_rotate(smpa2dmpa, b_smpa);

%------------------------------------%
% Rotate to GSE                      %
%------------------------------------%
	if ~isempty(attitude)
		% Transform into GEI
		gei2despun = mms_fdoa_xgei2despun(attitude, t, 'L');
		despun2gei = permute(gei2despun, [2, 1, 3]);
		b_gei      = mrvector_rotate(despun2gei, b_dmpa);

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
	else
		warning('FG:GSE', 'No attitude data. Cannot transform into GSE.');
		b_gse = zeros(size(b_smpa));
	end
end