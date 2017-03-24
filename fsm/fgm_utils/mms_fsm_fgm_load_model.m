%
% Name
%   mms_fsm_fgm_load_model
%
% Purpose
%   Find a MAT file containing a filter model for FGM and load the model.
%
% Calling Sequence
%   MODEL = mms_fsm_fgm_load_model(ROOT, SC, INSTR, RANGE, MODE)
%     Load one of David Fischer's filter models that lives in directory ROOT.
%     The model applies to MMS spacecraft SC, instrument INSTR, and operating
%     mode MODE, and field range RANGE. Values for MODE and RANGE are:
%       RANGE                  MODE
%         0 = Lo-Range          0 = ADCA (afg)   DEC64 (dfg)
%         1 = Hi-Range          1 = ADCB (afg)   DEC32 (dfg)
%
% Parameters
%   ROOT            in, required, type = char
%   SC              in, required, type = char
%   INSTR           in, required, type = char
%   RANGE           in, required, type = int8
%   MODE            in, required, type = int8
%
% Returns
%   MODEL           out, required, type = struct
%
% MATLAB release(s) MATLAB 9.0.0.341360 (R2016a)
% Required Products None
%
% History:
%   2017-01-31      Adapted from David Fischer's merge_snippet.m - MRA
%
function model = mms_fsm_fgm_model(root, sc, instr, range, mode)
	
	% Create the compensation model file name
	%   - File name strings have FS=128 always
	%   - File path has FSNEW=X, but we will remove the directory anyway
	isc    = str2num( sc(4) );
	fmodel = createmodelfilestring( instr, isc, range, mode, 128, 0 );
	
	% Pre-pend the file path
	%   - Remove a leading subdirectory
	[path, name, ext] = fileparts( fmodel );
	if nargin() == 5
		fmodel = fullfile( root, [name ext] );
	else
		fmodel = [name ext];
	end
	
	% Check if model exists
	if exist(fmodel, 'file') == 2
		% Print information about the model
		mrfprintf('logtext', 'Loading model:' );
		mrfprintf('logtext', ['    ' fmodel] );
	else
		mrfprintf('logtext', 'Model file does not exist:' );
		mrfprintf('logtext', ['    ' fmodel] );
		error( 'Cannot load FGM model.' );
	end
	
	% Load instrument model file for AFG/DFG
	%   - File name strings have FS=128 always
	load( fmodel );
	model = xfgmodel;
	clear xfgmodel
end