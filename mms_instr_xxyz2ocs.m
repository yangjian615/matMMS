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
%   INSTRUMENT      in, required, type=char
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-03-22      Written by Matthew Argall
%   2015-04-16      Magnetometer 123 and XYZ systems are now directly
%                     related to BCS. Removed OMB, as it is with respect
%                     to SMPA. - MRA
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
	euler_angles('AFG_XYZ')       = { [  45.0,  0.0,  -90.0 ], 'ZXY' };
	euler_angles('AFG_123')       = { [ 135.0,  0.0,    0.0 ], 'ZXY' };
	euler_angles('AFG_BOOM')      = { [ 135.0,  0.0,    0.0 ], 'ZXY' };
	euler_angles('BCS')           = { [   0.0,  0.0,    0.0 ], 'ZXY' };
	euler_angles('DFG_XYZ')       = { [-135.0,  0.0,  -90.0 ], 'ZXY' };
	euler_angles('DFG_123')       = { [ 135.0,  0.0,    0.0 ], 'ZXY' };
	euler_angles('DFG_BOOM')      = { [ -45.0,  0.0,    0.0 ], 'ZXY' };
	euler_angles('DSS')           = { [ -76.0,  0.0,    0.0 ], 'ZXY' };
	euler_angles('EDI1')          = { [ 221.0,  0.0,  -90.0 ], 'ZXY' };
	euler_angles('EDI1_GUN')      = { [ 221.0,  0.0,  -90.0 ], 'ZXY' };
	euler_angles('EDI1_DETECTOR') = { [ 221.0,  0.0,  -90.0 ], 'ZXY' };
	euler_angles('EDI2')          = { [  41.0,  0.0,  -90.0 ], 'ZXY' };
	euler_angles('EDI2_GUN')      = { [  41.0,  0.0,  -90.0 ], 'ZXY' };
	euler_angles('EDI2_DETECTOR') = { [  41.0,  0.0,  -90.0 ], 'ZXY' };
	euler_angles('OCS')           = { [   0.0,  0.0,    0.0 ], 'ZXY' };
	euler_angles('SC')            = { [   0.0,  0.0,    0.0 ], 'ZXY' };
	euler_angles('SCM_XYZ')       = { [ -45.0, 90.0,    0.0 ], 'ZXY' };
	euler_angles('SCM_123')       = { [ 135.0,  0.0,    0.0 ], 'ZXY' };
	euler_angles('SCM_BOOM')      = { [ 135.0,  0.0,    0.0 ], 'ZXY' };
	euler_angles('SDP1')          = { [ -60.0,  0.0,  180.0 ], 'ZXY' };
	euler_angles('SDP2')          = { [ 120.0,  0.0,  180.0 ], 'ZXY' };
	euler_angles('SDP3')          = { [  30.0,  0.0,  180.0 ], 'ZXY' };
	euler_angles('SDP4')          = { [ 210.0,  0.0,  180.0 ], 'ZXY' };

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
end