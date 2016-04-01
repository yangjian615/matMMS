%
% NAME
%   mms_dissect_filename
%
% PURPOSE
%   Dissect an MMS file name. The file name format is:
%     scId_instrumentId_mode_dataLevel_optionalDataProductDescriptor_startTime_vX.Y.Z.cdf
%
% Calling Sequence:
%   [SC, INSTR, MODE, LEVEL, TSTART, VERSION] = mms_dissect_filename(FILENAME);
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
%     >> filename = 'mms3_dfg_hr_l2p_duration-1h1m_20010704_v1.0.0.cdf'
%     >> [sc, instr, mode, level, desc, tstart, version] = mms_dissect_filename(filename);
%     >> fprintf('SC: %s\nINSTR: %s\nMODE: %s\nLEVEL: %s\nDESC: %s\nTSTART: %s\nVERSION: %s\n', sc, instr, mode, level, desc, tstart, version);
%         SC: mms3
%         INSTR: dfg
%         MODE: hr
%         LEVEL: l2p
%         DESC: duration-1h1m
%         TSTART: 20010704
%         VERSION: v1.0.0
%
%   Dissect two file names:
%     >> filenames = {'mms3_dfg_hr_l2p_duration-1h1m_20010704_v1.0.0.cdf' ...
%                     'mms3_dfg_hr_l2p_duration-1h1m_20010704_v1.0.0.cdf'}
%     >> [sc, instr, mode, level, desc, tstart, version] = mms_dissect_filename(filenames);
%     >> fprintf('SC: %s\nINSTR: %s\nMODE: %s\nLEVEL: %s\nDESC: %s\nTSTART: %s\nVERSION: %s\n', sc, instr, mode, level, desc, tstart, version);
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
function [sc, instr, mode, level, tstart, version, desc] = mms_dissect_filename(filename)
	
	% Extract the directory
	if iscell(filename)
		[path, fname, ext] = cellfun(@fileparts, filename, 'UniformOutput', false);
		fname              = strcat(fname, ext);
	else
		[path, fname, ext] = fileparts(filename);
		fname              = [fname ext];
	end
	
	% Form a recular expression to take apart the file name.
	parts = regexp(fname, ['(mms[1-4])_', ...                % Spacecraft ID
	                       '([a-z-]+)_', ...                 % Instrument ID
	                       '([a-z0-9]+)_', ...               % Instrument Mode
	                       '([a-z0-4]+)_', ...               % Data Level
	                       '([a-zA-Z0-9-]*)_?', ...          % Optional Descriptor
	                       '(20[0-9]{2}[0-9]{2}[0-9]{2}[0-9]*)_', ...  % Start Time
	                       '(v[0-9]+\.[0-9]+\.[0-9]+)', ...  % Version
	                       '\.cdf'], ...                     % Extension
	                       'tokens');

	% Make sure the file name is dissectable.
	assert(isempty(parts) == 0, ['Filename not recognized: "' filename '".']);
	
	% Un-nest the cells.
	parts = vertcat(parts{:});
	if iscell(filename)
		parts = vertcat(parts{:});
	
		% Extract the parts as row vectors
		sc      = parts(:,1)';
		instr   = parts(:,2)';
		mode    = parts(:,3)';
		level   = parts(:,4)';
		desc    = parts(:,5)';
		tstart  = parts(:,6)';
		version = parts(:,7)';
	
	% A single file name was given.
	else
		sc      = parts{1};
		instr   = parts{2};
		mode    = parts{3};
		level   = parts{4};
		desc    = parts{5};
		tstart  = parts{6};
		version = parts{7};
	end
end
