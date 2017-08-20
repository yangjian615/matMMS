%
% NAME
%   mms_dissect_filename
%
% PURPOSE
%   Dissect an MMS file name. The file name format is:
%     scId_instrumentId_mode_dataLevel_optionalDataProductDescriptor_startTime_vX.Y.Z.cdf
%
% Calling Sequence:
%   [SC, INSTR, MODE, LEVEL, DESC, TSTART, VERSION] = mms_dissect_filename(FILENAME);
%     Dissect an MMS data file name FILENAME and return the spacecraft ID,
%     SC, instrument ID, INSTR, data level, LEVEL, optional descriptor,
%     DESC, start time of the file TSTART, and file version number,
%     VERSION. Note that if DESC is not present in the file name, an empty
%     string will be returned.
%
%   [__] = mms_dissect_filename(FILENAME);
%     Dissect a cell array of file names. In this case all outputs are cell
%     arrays the same length as FILENAME.
%
% :Examples:
%   Dissect a DFG file name:
%     >> varname = 'mms2_dfg_lo_l2pre_G';
%     >> [sc, instr, name, desc] = mms_dissect_varname(varname);
%     >> fprintf('SC: %s\nINSTR: %s\nNAME: %s\nDESC: %s\n', sc, instr, name, desc);
% 				SC: mms2
% 				INSTR: dfg
% 				NAME: lo
% 				DESC: l2pre
%
%   Dissect two file names:
%     >> varnames = {'mms2_dfg_lo_l2pre_G' ...
%                    'mms2_dfg_lo_l2pre_dTheta' ...
%                    'mms2_dfg_lo_l2pre_dPhi' ...
%                    'mms2_dfg_lo_l2pre_U3' ...
%                    'mms2_dfg_lo_l2pre_O'}
%     >> [sc, instr, mode, level, desc, tstart, version] = mms_dissect_filename(filenames);
%     
%
% :Params:
%   FILENAME:       in, required, type=char/cell
%
% :Returns:
%   SC:             out, required, type=same as FILENAME
%   INSTR:          out, optional, type=same as FILENAME
%   MODE:           out, optional, type=same as FILENAME
%   LEVEL:          out, optional, type=same as FILENAME
%   DESC:           out, optional, type=same as FILENAME
%   TSTART:         out, optional, type=same as FILENAME
%   VERSION:        out, optional, type=same as FILENAME
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-03-22      Written by Matthew Argall
%
function [sc, instr, name, desc] = mms_dissect_varname(varname)
    
    % Form a recular expression to take apart the file name.
    parts = regexp(varname, ['(mms[1-4])_', ...              % Spacecraft ID
                              '([a-z-]+)_', ...              % Instrument ID
                              '([a-z0-4]+)_?', ...           % Param Name
															'([a-zA-Z0-9_]*)'], ...        % Optional Descriptor
                              'tokens');
    
    % Make sure the file name is dissectable.
    assert(isempty(parts) == 0, ['Filename not recognized: "' varname '".']);
		
		% Simplify getting retults
		parts = vertcat(parts{:});
		if iscell(varname)
			parts = vertcat(parts{:});
			
			% Extract the parts
			sc    = parts(:,1);
			instr = parts(:,2);
			name  = parts(:,3);
			desc  = parts(:,4);
		
		% A single file name was given.
		else
			sc    = parts{1};
			instr = parts{2};
			name  = parts{3};
			desc  = parts{4};
		end
    
end
