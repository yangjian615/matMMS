%
% Name
%   mms_instr_xocs2mpa
%
% Purpose
%   Create a rotation matrix from the observator coordinate system (OCS) to
%   the major principle axis (MPA) system.
%
%   Note that OCS differs from the spacecraft body coordinate system (BCS)
%   by a translation along z-OCS, and so does not affect a coordinate
%   transformation.
%
% Calling Sequence
%   OCS2MPA = mms_instr_xocs2mpa(SC, TSTART, ATT_DIR)
%     Create a rotation matrix from the OCS to MPA. Using the MMS
%     spacecraft SC and start time TSTART, extract the MPA z-axis from the
%     definitive attitude files found in directory ATT_DIR.
%
%   [OCS2MPA, zMPA] = mms_instr_xocs2dmpa(__)
%     Also return the z-MPA axis as it is viewed from OCS.
%
% Parameters
%   SC              in, required, type=char
%   TSTART          in, optional, type=char
%   ATT_DIR         in, optional, type=char
%
% Returns
%   DESPUN          out, required, type=3xN float
%
% See Also
%   MrTimeParser.m
%   mms_anc_read_defatt.m
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-06      Written by Matthew Argall
%
function [ocs2mpa, zMPA] = mms_instr_xocs2mpa(sc, tstart, att_dir)

%------------------------------------%
% Find Definitive Attitude File      %
%------------------------------------%
	% MMS#_DEFATT_%Y%D_%Y%D.V00

	% Transform the start time into a file time.
	att_start    = MrTimeParser(tstart, '%Y-%M-%dT%H:%m:%S', '%Y%D');
	
	% Build the file name and search for it.
	fname_test = fullfile(att_dir, [ upper(sc) '_DEFATT_' att_start '*' ]);
	defatt     = dir(fname_test);
	fname      = fullfile(att_dir, defatt.name);
	
	% Make sure the file exists
	assert( ~isempty( fname ), ...
            ['Definitive attitude file not found or does not exist: "' fname_test '".']);
	
%------------------------------------%
% Create Rotation to MPA             %
%------------------------------------%
	% Read the attitude data
	[~, ~, zMPA]  = mms_anc_read_defatt( fname );
	
	% Create the matrix to transform from OCS into DMPA.
	%   - z_hat = zMPA
	%   - x_hat = [0 1 0] x z_hat    -- [0 1 0] is yOCS
	%   - y_hat = z_hat   x x_hat
	ocs2mpa       = zeros(3, 3);
	ocs2mpa(3, :) = zMPA;
	ocs2mpa(1, :) = cross([0, 1, 0], zMPA);
	ocs2mpa(2, :) = cross( ocs2mpa(3, :), ocs2mpa(1, :) );
end