%
% Name
%   mms_edp_read_ql
%
% Purpose
%   Read Quick-Look (QL) electric field double probe (EDP) data from MMS.
%
% Calling Sequence
%   FG_QL = mms_fg_read_ql(FILES)
%     Read electric field double probe Quick-Look (QL) data from
%     files with file names FILES.
%
%   EDP_QL = mms_edp_read_ql(FILES, TSTART, TEND)
%     Read electric field double probe Quick-Look (QL) data from
%     files with file names FILES between the time interval beginning at
%     TSTART and ending at TEND. TSTART and TEND should be ISO-8601
%     strings: yyyy-mm-ddThh:mm:ss
%
% Parameters
%   FILES           in, required, type = char/cell
%   TSTART          in, required, type = char, default = ''
%   TEND            in, required, type = char, default = ''
%
% Returns
%   EDP_QL          out, required, type=struct
%                   Fields are:
%                       'tt2000'      - TT2000 epoch times
%                       'E_DSL'       - Electric field in DSL coordinates
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-06-07      Written by Matthew Argall
%
function edp_ql = mms_fg_read_ql(files, tstart, tend)

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

	% Instr, Level, Mode
	assert( min( strcmp(instr, 'edp') == 1 ), 'Only EDP files are allowed.' );
	assert( min( strcmp(level, 'ql') ) == 1, 'Only QL files are allowed.' );
	if nFiles > 1
		assert( min( strcmp(instr, instr{1}) ) == 1, 'All files must be from the same instrument.' );
		assert( min( strcmp(mode, mode{1}) )   == 1, 'All files must have the same telemetry mode.' );
	else

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
	e_dsl_name = mms_construct_varname(sc, instr, 'dce', 'xyz_dsl');

	% Read the magnetometer data
	[e_dsl, t] = MrCDF_nRead(files, e_dsl_name, 'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);

	% Transpose the data to be row vectors.
	edp_ql = struct( 'tt2000', t,     ...
	                 'e_dsl',  e_dsl );
end