%
% Name
%   mms_fdoa_read_defatt
%
% Purpose
%   Read MMS Attitude data from FDOA ASCII files.
%
% Calling Sequence
%   ATTITUDE = mms_fdoa_read_defatt(SC, TSTART, TEND, ATT_DIR)
%     Return structure of attitude data ATTITUDE, from MMS spacecraft SC
%     (e.g., 'mms2') between the times TSTART and TEND. Attitude data
%     files will be searched for in ATT_DIR.
%
%   [ATTITUDE, ATT_HDR] = mms_fdoa_read_defatt(__)
%     Also return a structure of header information from all files read.
%
% Parameters
%   SC              in, required, type=char
%   TSTART          in, required, type=char
%   TEND            in, required, type=char
%   ATT_DIR         in, required, type=char
%
% Returns
%   ATTITUDE         out, required, type=struct
%                    Fields are:
%                      'UTC'  -  Time tags in UTC
%                      'TAI'  -  Time tags in TAI as seconds since MMS reference epoch
%                      'q'    -  Quaternions
%                      'w'    -  Spin rate and phase.
%                      'z'    -  Body z-axis components and phase.
%                      'L'    -  Angular momentum axis and phase.
%                      'P'    -  Principle axis and phase.
%                      'Nut'  -  Nutation information.
%                      'QF'   -  Quality flag
%   ATT_HDR          out, required, type=struct
%
% Examples
%   Given the files
%     MMS2_DEFATT_2015078_2015079.V00
%     MMS2_DEFATT_2015079_2015080.V00
%
%   Read all the data from
%     >> sc      = 'mms2';
%     >> tstart  = '2015-03-20T00:00:00Z';
%     >> tend    = '2015-03-21T00:00:00Z';
%     >> att_dir = '/Users/argall/Documents/Work/Data/MMS/Ephemeris/'
%     >> data    = mms_fdoa_read_defatt(sc, tstart, tend, att_dir)
%        data    = UTC: {278887x1 cell}
%                  TAI: [278887x1 double]
%                    q: [278887x4 double]
%                    w: [278887x4 double]
%                    z: [278887x3 double]
%                    L: [278887x3 double]
%                    P: [278887x3 double]
%                  Nut: [278887x1 double]
%                   QF: {278887x1 cell}
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-06      Written by Matthew Argall
%
function [attitude, att_hdr] = mms_fdoa_read_defatt(sc, tstart, tend, att_dir)

%------------------------------------%
% Find Definitive Attitude File      %
%------------------------------------%
	% MMS#_DEFATT_%Y%D_%Y%D.V00
	fname = MrFile_Search( fullfile(att_dir, [upper(sc) '_DEFATT_%Y%D_%Y%D.V*']), ...
	                      'VersionRegex', 'V([0-9]{2})', ...
	                      'TStart',       tstart, ...
	                      'TEnd',         tend, ...
	                      'TimeOrder',    '%Y%D' );
	nFiles = length(fname);
	
	% Make sure the file exists
	assert( ~isempty( fname ), ...
	        ['Definitive attitude file not found or does not exist: "' fname '".']);

%------------------------------------%
% Read Header from Each File         %
%------------------------------------%
	% Read the header
	att_hdr = mms_fdoa_read_defatt_header(fname{1});
	
	% Step through the rest of the files
	for ii = 2 : nFiles
		temp = mms_fdoa_read_defatt_header(fname{ii});
		
		% Append the header information
		att_hdr.('att_err')    = [ att_hdr.('att_err'),    temp.('att_err')    ];
		att_hdr.('DataStart')  = [ att_hdr.('DataStart'),  temp.('DataStart')  ];
		att_hdr.('InertiaCal') = [ att_hdr.('InertiaCal'), temp.('InertiaCal') ];
		att_hdr.('StartTime')  = [ att_hdr.('StartTime'),  temp.('StartTime')  ];
		att_hdr.('StopTime')   = [ att_hdr.('StopTime'),   temp.('StopTime')  ];
		att_hdr.('rate_err')   = [ att_hdr.('rate_err'),   temp.('rate_err')   ];
		att_hdr.('zMPA')       = [ att_hdr.('zMPA'),       temp.('zMPA')       ];
		
		% Make sure the MPA axis has not changed
		if ~min( att_hdr.('zMPA')(:,ii) == temp.('zMPA') )
			warning('MMS_FDOA_Read_DefAtt:MPA', 'zMPA has changed.')
		end
	end

%------------------------------------%
% Read the Data                      %
%------------------------------------%
	% Describe the file
	column_names = {'UTC', 'TAI', 'q', 'q', 'q', 'q', 'w', 'w', 'w', 'w', ...
	              'z', 'z', 'z', 'L', 'L', 'L', 'P', 'P', 'P', 'Nut', 'QF'};
	data_start = att_hdr.('DataStart')(1);

	% Read the files
	attitude = MrFile_Read_nAscii(fname, ...
	                              'ColumnNames', column_names, ...
	                              'DataStart',   data_start, ...
	                              'nFooter',     1);

%------------------------------------%
% Time Interval                      %
%------------------------------------%
	% Convert input times to TT2000
	trange = MrTimeParser({tstart, tend}, '%Y-%M-%dT%H:%m:%S', '%Y-%M-%dT%H:%m:%S.%1%2%3');
	trange = spdfparsett2000(trange);
	
	% Convert TAI to TT2000
	tt2000 = mms_fdoa_epoch2tt2000(attitude.('TAI'), 'AttTAI', true);
	
	% Indices to keep
	irange    = zeros(1, 2);
	irange(1) = find( tt2000 >= trange(1), 1, 'first' );
	irange(2) = find( tt2000 <= trange(2), 1, 'last' );
	
	% Number of data returned
	names   = fieldnames(attitude);
	nFields = length(names);
	
	% Get rid of unwanted data
	for ii = 1 : nFields
		attitude.( names{ii} ) = attitude.( names{ii} )(irange(1):irange(2), :);
	end
end


%
% Name
%   mms_fdoa_read_defatt_header
%
% Purpose
%   Parse information for the header of MMS definitive attitude files.
%
% Calling Sequence
%   ATT_HDR = mms_fdoa_read_defatt_header(FILENAME)
%     Parse information from the header in file FILENAME and return a
%     structure ATT_HDR.
%
% Parameters
%   FILENAME        in, required, type=char
%
% Returns
%   ATT_HDR         out, required, type=structure
%                   A structure with the following fields:
%                     'att_err'        EKF mean attitude error
%                     'DataStart'      Line at which data begins
%                     'InertialCal'    Calibration file used with attitude data.
%                     'rate_err'       EKF mean rate error
%                     'StartTime'      First time in the file (converted to TT2000)
%                     'StopTime'       Last time in the file (converted to TT2000)
%                     'zMPA'           The z-MPA axis in BCS
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-09      Written by Matthew Argall
%
function att_hdr = mms_fdoa_read_defatt_header(filename)
	%
	% Fields to look for
	%   COMMENT  Latest inertia tensor calibration file = MMS2_INERTIACAL_2015064.V03
	%
	%   COMMENT  Major principal axis of inertia in BCS = -0.00004885  -0.00009667   0.99999999
	%
	%   COMMENT  Estimated EKF mean attitude and rate errors (1-sigma) about the body X,Y,Z axes are:
	%   COMMENT  0.0034  0.0034  0.0087 (deg)    0.003  0.003  0.004 (deg/s)
	%
	%   START_TIME      = 2015-079T00:59:04.000
	%
	%   STOP_TIME       = 2015-080T00:39:04.000
	%
	%   DATA_START
	%

	% Open the file
	fID = fopen(filename);

	% Where does data start?
	data_start  = 0;
	nHeader     = 0;
	errors_next = false;
	
	% Go until the end of the header
	while data_start == 0
		% Get a line of data
		line = fgetl(fID);
		
		%
		% Check for flags of interest
		%
		
		% INERTIALCAL
		if ~isempty(regexp(line, 'INERTIACAL', 'once'))
			parts      = regexp(line, '[ ]+', 'split');
			inertiacal = parts{end};
		
		% MPA
		elseif ~isempty(regexp(line, 'Major principal axis', 'once'))
			parts = regexp(line, '[ ]+', 'split');
			zMPA  = cellfun(@str2num, parts(end-2:end) )';
		
		% EKF Errors
		elseif ~isempty(regexp(line, 'Estimated EKF.*errors', 'once'))
			errors_next = true;
			
		% START_TIME
		elseif ~isempty(regexp(line, 'START_TIME', 'once'))
			parts      = regexp(line, '[ ]+', 'split');
			start_time = mms_fdoa_epoch2tt2000( parts{end} );
		
		% END_TIME
		elseif ~isempty(regexp(line, 'STOP_TIME', 'once'))
			parts     = regexp(line, '[ ]+', 'split');
			stop_time = mms_fdoa_epoch2tt2000( parts{end} );
		
		% DATA_START
		elseif ~isempty(regexp(line, 'DATA_START', 'once'))
			data_start = nHeader + 3;  % Must account for the nHeader at end of loop
			nHeader    = nHeader + 1;
		
		% ERRORS
		elseif errors_next
			parts    = regexp(line, '[ ]*', 'split');
			att_err  = cellfun(@str2num, parts(end-3:end-1) )';
			rate_err = cellfun(@str2num, parts(end-7:end-5) )';
			
			% Do not read errors any more.
			errors_next = false;
		end
		
		% Next line
		nHeader = nHeader + 1;
	end

	% Close the file
	fclose(fID);
	
	% Create a structure
	att_hdr = struct( 'att_err',     att_err,    ...
	                  'DataStart',   data_start, ...
	                  'InertiaCal',  inertiacal, ...
	                  'StartTime',   start_time, ...
	                  'StopTime',    stop_time,  ...
	                  'rate_err',    rate_err,   ...
	                  'zMPA',        zMPA );
end

% CCSDS_AEM_VERS = 1.0
% CREATION_DATE  = 2015-080T16:21:19
% ORIGINATOR     = GSFC
%  
% META_START
% COMMENT  This file includes definitive attitude data displayed as:
% COMMENT    time (UTC), elapsed seconds since reference epoch, quaternion (ECI-to-BCS),
% COMMENT    X,Y,Z rotation rate components (deg/s) (instantaneous spin axis in body frame),
% COMMENT    w-phase (Sun-to-body-X dihedral angle about rotation rate vector) (deg),
% COMMENT    right ascension (deg) and declination (deg) of body Z-axis,
% COMMENT    Z-phase (Sun-to-body-X dihedral angle about body Z-axis) (deg),
% COMMENT    right ascension (deg) and declination (deg) of angular momentum (L),
% COMMENT    L-phase (Sun-to-body-X dihedral angle about angular momentum vector L) (deg),
% COMMENT    right ascension (deg) and declination (deg) of major principal axis (P),
% COMMENT    P-phase (Sun-to-body-X dihedral angle about major principal axis P) (deg),
% COMMENT    nutation angle (deg),
% COMMENT    and quality flag: EKF=good extended Kalman filter solution with star tracker data
% COMMENT                      CNV=filter not yet converged
% COMMENT                      SUN=no tracker data; spin phase obtained from Sun sensor data
% COMMENT                      INT=no tracker or Sun sensor data; phase interpolated from neighboring tracker data
% COMMENT                      BAD=no tracker or Sun sensor data for a time span too large for interpolation
% COMMENT
% COMMENT  REPORT SOURCE = AGS OPS
% COMMENT
% COMMENT  Latest inertia tensor calibration file = MMS2_INERTIACAL_2015064.V03
% COMMENT
% COMMENT  Major principal axis of inertia in BCS = -0.00004885  -0.00009667   0.99999999
% COMMENT
% COMMENT  Estimated EKF mean attitude and rate errors (1-sigma) about the body X,Y,Z axes are:
% COMMENT  0.0034  0.0034  0.0087 (deg)    0.003  0.003  0.004 (deg/s)
% COMMENT
% OBJECT_NAME     = MMS2
% OBJECT_ID       = 2015-011B
% REF_FRAME_A     = J2000
% REF_FRAME_B     = BCS
% ATTITUDE_DIR    = A2B
% ATTITUDE_TYPE   = QUATERNION/RATE/SPIN/NUTATION
% QUATERNION_TYPE = LAST
% TIME_SYSTEM     = UTC
% COMMENT           The time tags are given twice. First as calendar UTC time,
% COMMENT           and then as TAI, expressed as elapsed SI seconds since
% COMMENT           the mission reference epoch, 1958-001T00:00:00 UTC.
% COMMENT
% START_TIME      = 2015-079T00:59:04.000
% STOP_TIME       = 2015-080T00:39:04.000
% META_STOP
%  
% DATA_START
% COMMENT   Time (UTC)    Elapsed Sec      q1       q2       q3       qc     wX     wY     wZ   w-Phase   Z-RA    Z-Dec Z-Phase   L-RA    L-Dec L-Phase   P-RA    P-Dec P-Phase    Nut    QF
% 2015-079T00:59:04.312 1805504379.312  0.26383  0.21707 -0.58189  0.73802  0.000  0.000 18.000 284.283 271.192  50.044 284.283 271.188  50.045 284.283 271.183  50.046 284.283   0.003  CNV