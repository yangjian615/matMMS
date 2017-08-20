%
% Name
%   mms_scm_l2_create
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
%   [T, B_GSE] = mms_scm_l2_create(SC_FILES, CAL_FILE, TSTART, TEND)
%     Find and read search coil magnetometer and transfer function
%     data from files SC_FILES and CAL_FILE during the time interval
%     [TSTART, TEND]. Return the data and its time tags B_GSE and T.
%
%   [..., B_DMPA] = mms_scm_l2_create(__)
%     Also return the magnetic field in DMPA.
%
%   [..., B_SMPA] = mms_scm_l2_create(__)
%     Also return the magnetic field in SMPA.
%
%   [..., B_BCS] = mms_scm_l2_create(__)
%     Also return the magnetic field in BCS.
%
%   [__] = mms_scm_l2_create(__, 'ParamName', ParamValue)
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
%   2015-08-09      Direct warnings to 'logwarn' with mrfprintf. - MRA
%
function [t, b_gsm, b_gse, b_gei, b_dmpa, b_bcs, b_smpa, b_omb, b_123] ...
	= mms_scm_l2_create(sc_files, cal_file, tstart, tend, varargin)

%------------------------------------%
% Inputs                             %
%------------------------------------%
	% Defaults
	duration = 64;
	attitude = [];
	sunpulse = [];
	zMPA     = [];

	nOptArg = length(varargin);
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
				error( ['Unknown input "' varargin{ii} '".'] );
		end
	end
	
	% Must have attitude or sunpulse to at least despin
	assert( ~( isempty(sunpulse) && isempty(attitude) ), 'Sunpulse and/or attitude data must be given.')

%------------------------------------%
% Calibrated & Rotate to BCS         %
%------------------------------------%
	[t, b_bcs, b_smpa, b_omb, b_123] = mms_scm_l1b_create(sc_files, cal_file, tstart, tend, duration);

%------------------------------------%
% BCS --> SMPA                       %
%------------------------------------%
	if isempty(zMPA)
		% Attitude data
		mrfprintf( 'logwarn', 'MMS_SC_GSE:SMPA', ...
		           'MPA axis not given. Cannot transform SMPA -> BCS.');
		
		% Cannot transform to SMPA, so copy variables
		b_bcs = b_smpa;
	else
		% Build matrix
		bcs2smpa = mms_fdoa_xbcs2smpa( zMPA );
		smpa2bcs = permute(bcs2smpa, [2, 1, 3]);

		% Transform
		b_bcs = mrvector_rotate(smpa2bcs, b_smpa);
	end

%------------------------------------%
% Despin                             %
%------------------------------------%

	%
	% Assume the principle axis of inertia (z-MPA)
	% is the same as the angular momentum vector (L)
	%
	
	% Despin sunpulse or attitude
	%   - Use sunpulse if it was given, otherwise use attitude
	if ~isempty(sunpulse)
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
		[b_gsm, b_gse, b_gei] = mms_rot_despun2gsm(t, b_dmpa, attitude, 'P');
	else
		mrfprintf( 'logwarn', 'FG:GSE', 'No attitude data. Cannot transform into GSE.');
		b_gse = NaN(size(b_dmpa));
		b_gsm = NaN(size(b_dmpa));
	end
end