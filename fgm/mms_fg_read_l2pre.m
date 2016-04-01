%
% Name
%   mms_fg_read_ql
%
% Purpose
%   Read Quick-Look (QL) fluxgate magnetometer data from MMS.
%
% Calling Sequence
%   FG_QL = mms_fg_read_ql(FILES)
%     Read fluxgate level 1B data from files with file names FILES.
%
%   FG_QL = mms_fg_read_ql(FILES, TSTART, TEND)
%     Read fluxgate level 1B data from files with file names FILES
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
%   FG_QL           out, required, type=struct
%                   Fields are:
%                       'tt2000'   - TT2000 epoch times
%                       'B_BCS'    - Magnetic field in BCS coordinates
%                       'B_DMPA'   - Magnetic field in DMPA coordinates
%                       'B_GSE'    - Magnetic field in GSE coordinates
%                       'B_GSM'    - Magnetic field in GSM coordinates
%                       'POS_GSE'  - Spacecraft position in GSE coordinates
%                       'POS_GSM'  - Spacecraft position in GSM coordinates
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-11-19      Written by Matthew Argall
%
function fg_l2pre = mms_fg_read_l2pre(files, tstart, tend)

	% Defaults
	if nargin() < 2
		tstart = '';
	end
	if nargin() < 3
		tend = '';
	end

%------------------------------------%
% Read Mag Data                      %
%------------------------------------%
	% Number of files given
	if iscell(files)
		nFiles = length(files);
	else
		nFiles = 1;
	end

	% Dissect the file name
	[sc, instr, mode, level, ~, ~, optdesc] = mms_dissect_filename(files);

	% AFG or DFG
	assert( min( strcmp(instr, 'afg') | strcmp(instr, 'dfg') ) == 1, ...
	        'Only AFG and DFG files are allowed.' );

	% Instr, Level, Mode
	if nFiles > 1
		assert( min( strcmp(instr, instr{1}) ) == 1, 'All files must be from the same instrument.' );
		assert( min( strcmp(mode, mode{1}) )   == 1, 'All files must have the same telemetry mode.' );
	end
	assert( min( strcmp(level, 'l2pre') ) == 1, 'Only QL files are allowed.' );

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
	b_bcs_name  = mms_construct_varname(sc, instr, mode, [level '_bcs']);
	b_dmpa_name = mms_construct_varname(sc, instr, mode, [level '_dmpa']);
	b_gse_name  = mms_construct_varname(sc, instr, mode, [level '_gse']);
	b_gsm_name  = mms_construct_varname(sc, instr, mode, [level '_gsm']);
	x_gse_name  = mms_construct_varname(sc, 'pos', 'gse');
	x_gsm_name  = mms_construct_varname(sc, 'pos', 'gsm');

	% Read the magnetometer data
	b_bcs       = MrCDF_nRead(files, b_bcs_name,  'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	[b_dmpa, t] = MrCDF_nRead(files, b_dmpa_name, 'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	b_gse       = MrCDF_nRead(files, b_gse_name,  'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	b_gsm       = MrCDF_nRead(files, b_gsm_name,  'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	pos_gse     = MrCDF_nRead(files, x_gse_name,  'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	pos_gsm     = MrCDF_nRead(files, x_gsm_name,  'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);

	% Transpose the data to be row vectors.
	fg_l2pre = struct( 'tt2000',      t,       ...
	                   'b_bcs',       b_bcs,   ...
	                   'b_dmpa',      b_dmpa,  ...
	                   'b_gse',       b_gse,   ...
	                   'b_gsm',       b_gsm,   ...
	                   'pos_gse',     pos_gse, ...
	                   'pos_gsm',     pos_gsm );
end