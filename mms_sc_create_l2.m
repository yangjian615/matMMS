%
% Name
%   mms_sc_gse
%
% Purpose
%   Turn level 1A SCM data into level 2 quality data. L2 implies
%   despun data in the spacecraft reference frame (no VxB removal) and in
%   GSE coordinates.
%
%   Process:
%     1. Find, read, calibrate l1a data
%     2. Despin
%     3. Transform to GSE
%
% Calling Sequence
%   [T, B_GSE] = mms_sc_create_l2(SC_FILES, CAL_FILE, TSTART, TEND)
%     Find and read search coil magnetometer and transfer function
%     data from files SC_FILES and CAL_FILE during the time interval
%     [TSTART, TEND]. Return the data and its time tags B_GSE and T.
%
%   [..., B_DMPA] = mms_sc_create_l2(__)
%     Also return the magnetic field in DMPA.
%
%   [..., B_SMPA] = mms_sc_create_l2(__)
%     Also return the magnetic field in SMPA.
%
%   [..., B_BCS] = mms_sc_create_l2(__)
%     Also return the magnetic field in BCS.
%
%   [__] = mms_sc_create_l2(__, 'ParamName', ParamValue)
%     Any parameter name-value pair found below.
%
% Parameters
%   SC_FILES        in, required, type = char/cell
%   CAL_FILE        in, required, type = char
%   TSTART          in, optional, type = char
%   TEND            in, optional, type = char
%   'Attitude'      in, optional, type = struct, default = []
%                   Structure of definitive attitude data returned by mms_fdoa_read_defatt.m
%   'Duration'      in, optional, type = char default = 'DataDir'
%                   Time duration over which to calibrate data.
%   'SunPulse'      in, optional, type=struct, default=[]
%                   Structure of HK 101 sunpulse data returned by mms_dss_read_sunpulse.m
%   'zMPA'          in, optional, type=boolean, default=false
%                   The z-MPA axis used for rotating from BCS to SMPA coordinates.
%                     Available in the returned header structure from mms_fdoa_read_defatt.m
%
% Returns
%   T               out, required, type=1xN int64
%   B_GSE           out, required, type=3xN double
%   B_DMPA          out, optional, type=3xN double
%   B_SMPA          out, optional, type=3xN double
%   B_BCS           out, optional, type=3xN double
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-15      Written by Matthew Argall
%   2015-05-07      Added 'SunPulseDir'. - MRA
%   2015-06-21      Renamed from mms_sc_gse to mms_sc_create_l2. - MRA
%
function [t, b_gse, b_dmpa, b_smpa, b_bcs, b_cal, b_123] = mms_sc_create_l2(sc_files, cal_file, tstart, tend, varargin)

%------------------------------------%
% Inputs                             %
%------------------------------------%
	% Defaults
	duration = 600;
	attitude = [];
	sunpulse = [];
	zMPA     = [];

	nOptArg = length(varargin)
	for ii = 1 : 2 : nOptArg
		switch varargin{ii}
			case 'Attitude'
				attitude = varargin{ii+1};
			case 'SunPulse'
				sunpulse = varargin{ii+1};
			case 'Duration'
				duration = varargin{ii+1};
			case 'zMPA'
				zMPA     = varargin{ii+1};
			otherwise
				error( ['Unknown input "' varargin{ii+1} '".'] );
		end
	end
	
	% Must have attitude or sunpulse to at least despin
	assert(~isempty(sunpulse) && ~isempty(attitude), 'Sunpulse and/or attitude data must be given.')

%------------------------------------%
% Calibrated & Rotate to BCS         %
%------------------------------------%
	[t, b_bcs, b_cal, b_123] = mms_sc_create_l1b(sc_files, cal_file, tstart, tend, duration);

%------------------------------------%
% BCS --> SMPA                       %
%------------------------------------%
	if isempty(zMPA)
		% Attitude data
		warning('MMS_SC_GSE:SMPA', 'MPA axis not given. Cannot transform BCS -> SMPA.');
		
		% Cannot transform to SMPA, so copy variables
		b_smpa = b_bcs;
	else
		% Build matrix
		bcs2smpa = mms_fdoa_xbcs2smpa( zMPA );
		smpa2bcs = permute(bcs2smpa, [2, 1, 3]);

		% Transform
		b_smpa  = mrvector_rotate(smpa2bcs, b_bcs);
	end

%------------------------------------%
% Despin                             %
%------------------------------------%

	%
	% Assume the principle axis of inertia (z-MPA)
	% is the same as the angular momentum vector (L)
	%
	
	% Despin sunpulse or attitude
	if isempty(attitude)
		smpa2dmpa = mms_dss_xdespin( sunpulse, t );
	else
		smpa2dmpa = mms_fdoa_xdespin( attitude, t, 'L' );
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