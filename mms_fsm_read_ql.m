%
% Name
%   mms_fsm_read_ql
%
% Purpose
%   Read Quick-Look (QL) fluxgate-searchcoil-merged (FSM) data from MMS.
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
function fg_ql = mms_fsm_read_ql(files, tstart, tend)

	% Defaults
	if nargin() < 2
		tstart = '';
	end
	if nargin() < 3
		tend = '';
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
	assert( min( strcmp(level, 'ql') ) == 1, 'Only QL files are allowed.' );

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

	% Create variable names
	b_dmpa_name = mms_construct_varname(sc, instr, 'b', 'dmpa');
	b_fgm_name  = mms_construct_varname(sc, instr, 'b', 'dfg_dmpa');

	% Read the magnetometer data
	[b_fsm_dmpa, t_fsm] = MrCDF_nRead(files, b_dmpa_name, 'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	[b_fgm_dmpa, t_fgm] = MrCDF_nRead(files, b_fgm_name,  'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);

	% Transpose the data to be row vectors.
	fg_ql = struct( 'tt2000',      t_fsm,      ...
	                'tt2000_fgm',  t_fgm,      ...
	                'b_dmpa',      b_fsm_dmpa, ...
	                'b_fgm_dmpa',  b_fgm_dmpa );
end