%
% Name
%   mms_fsm_l2plus_read
%
% Purpose
%   Read Level 2 Plus (QL) fluxgate-searchcoil-merged (FSM) data from MMS.
%
% Calling Sequence
%   FSM_QL = mms_fsm_read_ql(FILES)
%     Read FSM quick-look data from files with file names FILES.
%
%   FSM_QL = mms_fsm_read_ql(FILES, TSTART, TEND)
%     Read fluxgate quick-look data from files with file names FILES
%     between the time interval beginning at TSTART and ending at
%     TEND. TSTART and TEND should be ISO-8601 strings: 
%     yyyy-mm-ddThh:mm:ss
%
% Parameters
%   FILES           in, required, type = char/cell
%   TSTART          in, required, type = char, default = ''
%   TEND            in, required, type = char, default = ''
%
% Returns
%   FSM_QL          out, required, type=struct
%                   Fields are:
%                       'tt2000'      - FSM TT2000 epoch time tags
%                       'tt2000_fgm'  - FGM TT2000 epoch time tags
%                       'B_DMPA'      - Merged magnetic field in DMPA coordinates
%                       'B_FGM_DMPA'  - FGM magnetic field in DMPA coordinates
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-07-23      Written by Matthew Argall
%
function fsm_l2plus = mms_fsm_l2plus(arg1, arg2, arg3, arg4)

	% Defaults
	if nargin() < 2
		tstart = '';
	end
	if nargin() < 3
		tend = '';
	end

%------------------------------------%
% Check Inputs                       %
%------------------------------------%
	if nargin == 3
		files  = arg1;
		tstart = arg2;
		tend   = arg3;
	elseif nargin == 4
		sc     = arg1;
		mode   = arg2;
		tstart = arg3;
		tend   = arg4;
		
		[files, ~, str] = mms_file_search( sc, 'dfg-scm', mode, 'l2plus', ...
		                       'Directory', '/nfs/fsm/temp',  ...
		                       'OptDesc', 'fsm-split',        ...
		                       'TStart',  tstart,             ...
		                       'TEnd',    tend);
		assert(~isempty(files), ['Unable to find FSM L2PLUS files: "' str '".']);
	end


%------------------------------------%
% Read FSM Data                      %
%------------------------------------%
	% Number of files given
	if iscell(files)
		nFiles = length(files);
	else
		nFiles = 1;
	end

	% Dissect the file name
	[sc, instr, mode, level, ~, ~, optdesc] = mms_dissect_filename(files);

	% Check file names
	assert( min( strcmp(instr, 'dfg-scm') ) == 1, 'Only FSM files are allowed.' );
	if nFiles > 1
		assert( min( strcmp(mode, mode{1}) ) == 1, 'All files must have the same telemetry mode.' );
	else
	assert( min( strcmp(level,   'l2plus'   ) ) == 1, 'Only L2PLUS files are allowed.' );
	assert( min( strcmp(optdesc, 'fsm-split') ) == 1, 'Only FSM files are allowed.' );

	% We now know all the files match, so keep on the the first value.
	if nFiles > 1
		sc      = sc{1};
		instr   = instr{1};
		mode    = mode{1};
		level   = level{1};
		optdesc = optdesc{1};
	end

%------------------------------------%
% Read Mag Data                      %
%------------------------------------%
	% Variables names use underscores instead of hyphens
	vinstr = strrep(instr, '-', '_');

	% Create variable names
	b_dmpa_name = mms_construct_varname(sc, vinstr, 'b', 'dmpa');
	b_bcs_name  = mms_construct_varname(sc, vinstr, 'b', 'bcs');
	b_gse_name  = mms_construct_varname(sc, vinstr, 'b', 'gse');
	b_gsm_name  = mms_construct_varname(sc, vinstr, 'b', 'gsm');

	% Read the magnetometer data
	[b_dmpa, t_fsm] = MrCDF_nRead(files, b_dmpa_name, 'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	b_bcs           = MrCDF_nRead(files, b_bcs_name,  'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	b_gse           = MrCDF_nRead(files, b_gse_name,  'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	b_gsm           = MrCDF_nRead(files, b_gsm_name,  'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);

	% Transpose the data to be row vectors.
	fsm_l2plus = struct( 'tt2000', t_fsm,  ...
	                     'b_dmpa', b_dmpa, ...
	                     'b_bcs',  b_bcs,  ...
	                     'b_gse',  b_gse,  ...
	                     'b_gsm',  b_gsm   ...
	                    );
end