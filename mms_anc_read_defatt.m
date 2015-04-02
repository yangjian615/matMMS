%
% Name
%   mms_anc_read_defatt
%
% Purpose
%   Read definitive attitude file.
%
% Calling Sequence
%   [DEFATT] = mms_anc_read_deffat(FILENAME)
%     Read definitive attitude data from file FILENAME.
%
%   [DEFATT, META] = mms_anc_read_deffat(FILENAME)
%     Also return the file's metadata. Each line of the file is stored as
%     an element of a cell array. Includes all lines between META_START and
%     META_STOP.
%
%   [DEFATT, META, MPA] = mms_anc_read_deffat(FILENAME)
%     Also return the averate major principle axis, as viewed from the body
%     coordinate system (BCS).
%
% Parameters
%   FILENAME        in, required, type=char
%
% Returns
%   ATT             out, required, type=1x21 cell
%                   Cells are:
%                      1  Time           (UTC) 'yyyy-dddTHH:MM:SS.sss'
%                      2  Elapsed Sec    Elapsed seconds since reference epoch
%                      3  q1             Quaternion (ECI-to-BCS)
%                      4  q2             Quaternion (ECI-to-BCS)
%                      5  q3             Quaternion (ECI-to-BCS)
%                      6  qc             Quaternion (ECI-to-BCS)
%                      7  wX             Rotation rate components (rad/s) -- Instantaneous spin axis in BCS
%                      8  wY             Rotation rate components (rad/s) -- Instantaneous spin axis in BCS
%                      9  wZ             Rotation rate components (rad/s) -- Instantaneous spin axis in BCS
%                     10  w-Phase        Sun-to-body-X dihedral angle about rotation rate vector (deg)
%                     11  Z-RA           Right ascension of body Z-axis
%                     12  Z-Dec          Declination of body Z-axis
%                     13  Z-Phase        Sun-to-body-X dihedral angle about body Z-axis (deg)
%                     14  L-RA           Right ascension of body L-axis
%                     15  L-Dec          Declination of body L-axis
%                     16  L-Phase        Sun-to-body-X dihedral angle about angular momentum vector L (deg)
%                     17  P-RA           Right ascension of body P-axis
%                     18  P-Dec          Declination of body P-axis
%                     19  P-Phase        Sun-to-body-X dihedral angle about major principal axis P (deg)
%                     20  Nut            Nutation angle (deg)
%                     21  QF             Quality flag
%                                          EKF = good extended Kalman filter solution with star tracker data
%                                          CNV = filter not yet converged
%                                          SUN = no tracker data; spin phase obtained from Sun sensor data
%                                          INT = no tracker or Sun sensor data; phase interpolated from neighboring
%                                          BAD = no tracker or Sun sensor data for a time span too large for interpolation
%   META            out, optional, type=cell
%   MPA             out, optional, type=1x3 double
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-03-25      Written by Matthew Argall
%
function [att, meta, mpa] = mms_anc_read_defatt(filename)

	% File must exist.
	assert( exist(filename, 'file') == 2, ...
		      ['Attitude file does not exist: "' filename '".']);

	% Open the file
	fID = fopen(filename);

%------------------------------------%
% Skip Over Header                   %
%------------------------------------%
	
	% Find META_START
	line = '';
	while isempty(regexp(line, 'META_START', 'once'))
		line = fgetl(fID);
	end
	
	% Find META_STOP
	count   = 1;
	meta    = cell(1, 100);
	meta(:) = {''};
	while isempty(regexp(line, 'META_STOP', 'once'))
		line        = fgetl(fID);
		meta{count} = line;
		count       = count + 1;
	end
	meta = meta(1:count-1);
	
	% Find DATA_START
	while isempty(regexp(line, 'DATA_START', 'once'))
		line = fgetl(fID);
	end
	
	% Read the column headers
	line = fgetl(fID);

%------------------------------------%
% Read the Data                      %
%------------------------------------%
	nCols  = 21;
	format = ['%s ' repmat('%f ', 1, nCols-2) '%f'];
	att = textscan(fID, format, 'MultipleDelimsAsOne', true);
	
	% Close the file
	fclose(fID);

%------------------------------------%
% Major Principle Axis               %
%------------------------------------%
	if nargout() == 3
		% Test each line of meta data to see if it contains the MPA
		mpa_regex    = cell(1, count-1);
		mpa_regex(:) = {'Major principal axis'};
		tf_mpa       = cellfun(@regexp, meta, mpa_regex, 'UniformOutput', false);
		
		% Extract the line with the MPA in it.
		impa     = find(cellfun(@isempty, tf_mpa) == 0, 1, 'first');
		mpa_line = meta{impa};
		
		% Extract the major principle axis
		%   - COMMENT  Major principal axis of inertia in BCS = -0.00008026  -0.00006849   0.99999999
		mpa_line = textscan(mpa_line, '%s %s %s %s %s %s %s %s %s %f %f %f', 'MultipleDelimsAsOne', true);
		mpa      = vertcat(mpa_line{end-2:end})';
	end
end



% CCSDS_AEM_VERS = 1.0
% CREATION_DATE  = 2015-075T02:39:50
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
% COMMENT  Latest inertia tensor calibration file = MMS2_INERTIACAL_2015075.V00
% COMMENT
% COMMENT  Major principal axis of inertia in BCS = -0.00008026  -0.00006849   0.99999999
% COMMENT
% COMMENT  Estimated EKF mean attitude and rate errors (1-sigma) about the body X,Y,Z axes are:
% COMMENT  0.0049  0.0049  0.0126 (deg)    0.003  0.003  0.004 (deg/s)
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
% START_TIME      = 2015-073T03:34:17.000
% STOP_TIME       = 2015-074T03:13:30.000
% META_STOP
%  
% DATA_START
% COMMENT   Time (UTC)    Elapsed Sec      q1       q2       q3       qc     wX     wY     wZ   w-Phase   Z-RA    Z-Dec Z-Phase   L-RA    L-Dec L-Phase   P-RA    P-Dec P-Phase    Nut    QF
% 2015-073T03:34:20.062 1804995295.062  0.23509  0.10334 -0.83693  0.48332  0.000 -0.000 18.000 251.081 233.736  60.239 251.081 233.732  60.241 251.080 233.727  60.244 251.079   0.004  CNV


