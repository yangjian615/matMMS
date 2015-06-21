%
% Name
%   mms_sc_gse
%
% Purpose
%   Turn level 1A AFG and DFG data into level 2 quality data. L2 implies
%   despun data in the spacecraft reference frame (no VxB removal) and in
%   GSE coordinates.
%
%   Process:
%     1. Find, read, calibrate l1a data
%     2. Despin
%     3. Transform to GSE
%
% Calling Sequence
%   [T, B_GSE] = mms_fg_create_l1a(FILES, HICAL_FILE, LOCAL_FILE)
%     Read, calibrate and rotate into GSE level 1A fluxgate data
%     from files with names FILES, using hi- and lo-range
%     calibration data found in files HICAL_FILE and LOCAL_FILE.
%
%   [T, B_GSE] = mms_fg_create_l1a(__, TSTART, TEND)
%     Return data within the time interval beginning at TSTART and
%     ending at TEND.
%
%   [__, B_DMPA, B_SMPA, B_OMB, B_123] = mms_fg_create_l1a(__)
%     Also return data in the DMPA, SMPA, OMB, and 123 coordinate systems.
%
%   [__] = mms_fg_create_l1a(..., 'ParamName', ParamValue)
%     Any parameter name-value pair found below.
%
% Parameters
%   FILES           in, required, type = char/cell
%   HICAL_FILE      in, required, type = char
%   LOCAL_FILE      in, required, type = char
%   TSTART          in, optional, type = char, default = ''
%   TEND            in, optional, type = char, default = ''
%   'Attitude'      in, optional, type = struct, default = []
%                   Structure of definitive attitude data returned by mms_fdoa_read_defatt.m
%   'SunPulse'      in, optional, type=struct, default=[]
%                   Structure of HK 101 sunpulse data returned by mms_dss_read_sunpulse.m
%
% Returns
%   T               out, required, type=1xN int64
%   B_GSE           out, required, type=3xN double
%   B_DMPA          out, optional, type=3xN double
%   B_BCS           out, optional, type=3xN double
%   B_SMPA          out, optional, type=3xN double
%   B_OMB           out, optional, type=3xN double
%   B_123           out, optional, type=3xN double
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-15      Written by Matthew Argall
%   2015-04-20      Added 'SunPulseDir' parameter.
%   2015-05-21      Take filenames as input, not pieces of filenames.
%
function [t, b_gse, b_dmpa, b_smpa, b_bcs, b_omb, b_123] = mms_fg_create_l1a(files, hiCal_file, loCal_file, tstart, tend, varargin)

%------------------------------------%
% Inputs                             %
%------------------------------------%
	attitude = [];
	sunpulse = [];
	nOptArgs = length(varargin);
	
	for ii = 1 : 2 : nOptArgs
		switch varargin{ii}
			case 'Attitude'
				attitude = varargin{ii+1};
			case 'SunPulse'
				sunpulse = varargin{ii+1};
			otherwise
				error(['Parameter name not recognized: "' varargin{ii+1} '".']);
		end
	end
	
	% Must have attitude or sunpulse to at least despin
	if isempty(sunpulse) && isempty(attitude)
		error( 'Attitude and/or Sunpulse data must be given.' );
	end

%------------------------------------%
% Calibrated Mag in BCS              %
%------------------------------------%
	switch nargout()
		case 7
			[t, b_bcs, b_smpa, b_omb, b_123] = mms_fg_create_l1a(files, hiCal_file, loCal_file, tstart, tend);
		case 6
			[t, b_bcs, b_smpa, b_omb] = mms_fg_create_l1a(files, hiCal_file, loCal_file, tstart, tend);
		case 5
			[t, b_bcs, b_smpa] = mms_fg_create_l1a(files, hiCal_file, loCal_file, tstart, tend);
		otherwise
			[t, ~, b_smpa] = mms_fg_create_l1a(files, hiCal_file, loCal_file, tstart, tend);
	end

%------------------------------------%
% Despin                             %
%------------------------------------%
	%
	% Assume the principle axis of inertia (z-MPA)
	% is the same as the angular momentum vector (L)
	%
	
	% Despin using attitude data
	if isempty(sunpulse)
		smpa2dmpa = mms_fdoa_xdespin( attitude, t, 'L', [upper(instr) '_123'] );
	
	% Despin using sunpulse directory
	else
		smpa2dmpa = mms_dss_xdespin( sunpulse, t );
	end

	% Rotate to DMPA
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