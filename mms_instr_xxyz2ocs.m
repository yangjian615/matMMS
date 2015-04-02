%
% Name
%   mms_instr_xxyz2ocs
%
% Purpose
%   Create a rotation matrix to rotate an instruments coordinate system
%   into the Observatory Coordinate System (OCS).
%
% Calling Sequence
%   XYZ2OCS = mms_instr_xxyz2ocs(INSTRUMENT);
%     Create a coordinate system transformation XYZ2OCS that rotates the
%     instrument coordinate system into OCS.
%
%   [XYZ2OCS, ANGLES, SEQUENCE] = mms_instr_xxyz2ocs(INSTRUMENT);
%     Return the euler rotation angles and the sequence in which they were
%     applied.
%
% Parameters
%   EUL             in, required, type=double
%   SEQUENCE        in, required, type=char, defualt='ZYX'
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-03-22      Written by Matthew Argall
%
function [xyz2ocs, angles, sequence] = mms_instr_xxyz2ocs(instrument)

	% If an instrument was given
	if nargin > 0
		% Ensure only one instrument name was given
		assert(ischar(instrument) && isrow(instrument), ...
		       'Input instrument must be a scalar string.');
		
		% Case insensitive -- all uppercase
		instr = upper(instrument);
	end

%-------------------------------------------------------
% Define Rotations to OCS //////////////////////////////
%-------------------------------------------------------
	%
	% NOTE!!
	%   - It is important to other programs that the rotation about z take place last --
	%       after the instrument z-axis has been aligned with the OCS z-axis. This allows
	%       the EUL(1) to provide information about relative angular offsets between
	%       instrument coordinate systems.
	%
	euler_angles = containers.Map('KeyType', 'char', 'ValueType', 'any');
	%
	% ORDER = 'ZXY' =>
	%   Az = eul(1)
	%   Ax = eul(2)
	%   Ay = eul(3)
	%
	%   x' = Az * Ay * Ax * x
	%             INSTRUMENT            EULER ANGLES           ORDER
	euler_angles('ADP1')          = { [   0.0,  0.0,    0.0 ], 'ZXY' };
	euler_angles('ADP2')          = { [   0.0,  0.0,  180.0 ], 'ZXY' };
	euler_angles('AFG_XYZ')       = { [ -90.0,  0.0,  -90.0 ], 'ZXY' };     % AFG_XYZ To AFG BOOM
	euler_angles('AFG_123')       = { [   0.0,  0.0,    0.0 ], 'ZXY' };     % AFG_123 To AFG_BOOM
	euler_angles('AFG_BOOM')      = { [ 135.0,  0.0,    0.0 ], 'ZXY' };
	euler_angles('BCS')           = { [   0.0,  0.0,    0.0 ], 'ZXY' };
	euler_angles('DFG_XYZ')       = { [ -90.0,  0.0,  -90.0 ], 'ZXY' };     % DFG_XYZ To DFG_BOOM
	euler_angles('DFG_123')       = { [ 180.0,  0.0,    0.0 ], 'ZXY' };     % DFG_123 To DFG_BOOM
	euler_angles('DFG_BOOM')      = { [ -45.0,  0.0,    0.0 ], 'ZXY' };
	euler_angles('DSS')           = { [ -76.0,  0.0,    0.0 ], 'ZXY' };
	euler_angles('EDI1')          = { [ 221.0,  0.0,  -90.0 ], 'ZXY' };
	euler_angles('EDI1_GUN')      = { [ 221.0,  0.0,  -90.0 ], 'ZXY' };
	euler_angles('EDI1_DETECTOR') = { [ 221.0,  0.0,  -90.0 ], 'ZXY' };
	euler_angles('EDI2')          = { [  41.0,  0.0,  -90.0 ], 'ZXY' };
	euler_angles('EDI2_GUN')      = { [  41.0,  0.0,  -90.0 ], 'ZXY' };
	euler_angles('EDI2_DETECTOR') = { [  41.0,  0.0,  -90.0 ], 'ZXY' };
	euler_angles('OCS')           = { [   0.0,  0.0,    0.0 ], 'ZXY' };
	euler_angles('OMB')           = { [ 225.0,  0.0,    0.0 ], 'ZXY' };
	euler_angles('SC')            = { [   0.0,  0.0,    0.0 ], 'ZXY' };
	euler_angles('SCM_XYZ')       = { [ 180.0, 90.0,    0.0 ], 'ZXY' };     % SCM_XYZ To SCM_BOOM
	euler_angles('SCM_123')       = { [ 180.0, 90.0,    0.0 ], 'ZXY' };     % SCM_123 To SCM_BOOM
	euler_angles('SCM_BOOM')      = { [ 135.0,  0.0,    0.0 ], 'ZXY' };
	euler_angles('SDP1')          = { [ -60.0,  0.0,  180.0 ], 'ZXY' };
	euler_angles('SDP2')          = { [ 120.0,  0.0,  180.0 ], 'ZXY' };
	euler_angles('SDP3')          = { [  30.0,  0.0,  180.0 ], 'ZXY' };
	euler_angles('SDP4')          = { [ 210.0,  0.0,  180.0 ], 'ZXY' };

	% This is how the rotations are given in the instrument manual.
% 	euler_angles('AFG123')     = {  0.0,  -90.0,  -90.0, [3,2,1] };     % AFG-to-AFG123
% 	euler_angles('DFG123')     = {  0.0,  -90.0,   90.0, [3,2,1] };     % DFG-to-DFG123
% 	euler_angles('SCM123')     = { 90.0,    0.0,  180.0, [3,1,2] };     % SCM-to-SCM123

%-------------------------------------------------------
% Names ////////////////////////////////////////////////
%-------------------------------------------------------
	% Return only the instrument names?
	if nargin == 0
		xyz2ocs = euler_angles.keys();
		return
	end

	% Case insensitive version
	tf_has = euler_angles.isKey(instr);
	if ~tf_has
		error( ['Invalid instruments given: "' instrument] );
	end

%-------------------------------------------------------
% Create Rotation Matrix ///////////////////////////////
%-------------------------------------------------------
	% Collect the data
	eulCell  = euler_angles(instr);
	angles   = eulCell{1};
	sequence = eulCell{2};
	
	% Get the rotation matrix that rotates XYZ into OCS.
	xyz2ocs = mreul2rotm(angles, sequence, 'Degrees', true);
	
%-------------------------------------------------------
% Finish Rotation to OCS ///////////////////////////////
%-------------------------------------------------------

	% Magnetometer mechanical and sensor frames must be rotated from MAG_BOOM to OCS.
	%    - e.g. 'AFG_123' and 'AFG_XYZ' must be rotated from 'AFG_BOOM' to 'OCS'
	%    - The process is the same for all three.
	if ~isempty(regexp(instr, '(AFG|DFG|SCM)_(XYZ|123)', 'once'))
		% Extract the instrument
		mag_name = regexp(instr, '^(AFG|DFG|SCM)', 'tokens');
		mag_name = mag_name{1}{1};
		
		% What we really have is the transformation to the MAG_BOOM system
		mag2mag_boom = xyz2ocs;
				
		% Finish rotation to OCS
		eulCell   = euler_angles( [mag_name '_BOOM'] );
		angles2   = eulCell{1};
		sequence2 = eulCell{2};
		
		mag_boom2ocs = mreul2rotm(angles2, sequence2, 'Degrees', true);
		xyz2ocs      = mag_boom2ocs * mag2mag_boom;
		
		% The booms lie in the OCS XY-plane, and this extra rotation is about
		% only the z-axis. Futhermore, since the rotation about Z is always
		% last, the final set of rotations is the orignal (from mag2mag_boom)
		% plus the rotation about z from mag_boom2ocs.
		angles(1) = angles(1) + angles2(1);
	end
end