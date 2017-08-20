%
% Name
%   mms_instr_xxyz2ocs
%
% Purpose
%   Create a rotation matrix to rotate one instrument's coordinate system
%   (CS) into another instrument's CS
%
% Calling Sequence
%   XYZ2INSTR = mms_instr_xxyz2instr(INSTR1, INSTR2);
%     Create a coordinate system transformation XYZ2INSTR that rotates
%     INSTR1's CS into that of INSTR2. INSTR1 and INSTR2 are names of
%     instruments on MMS; see mms_instr_xxyz2ocs.m for details.
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
function xyz2instr = mms_instr_xxyz2instr(instr1, instr2)

	% Get the rotation angles
	instr12ocs = mms_instr_xxyz2ocs(instr1);
	instr22ocs = mms_instr_xxyz2ocs(instr2);
	
	% For orthogonal matrices, A inverse = A transpose
	xyz2instr = instr22ocs' * instr12ocs;
end