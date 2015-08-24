%
% Name
%   mms_hk_read_0x10e
%
% Purpose
%   Read housekeeping 0x10e files that contain sensor temperatures.
%
% Calling Sequence
%   HK_10E = mms_hk_read_0x10e(FILES)
%     Read housekeeping L1B 10E data from files with file names FILES.
%
%   HK_10E = mms_hk_read_0x10e(FILES, TSTART, TEND)
%     Read housekeeping L1B 10E data from files with file names FILES
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
%   HK_10E          out, required, type=struct
%                   Fields are:
%                       'tt2000'       - TT2000 epoch times
%                       'afg_stemp'    - Sensor temperature for AFG
%                       'dfg_stemp'    - Sensor temperature for DFG
%                       'afg_etemp'    - Electronics temperature for AFG
%                       'dfg_etemp'    - Electronics temperature for DFG
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-07-27      Written by Matthew Argall
%
function hk_10e = mms_hk_read_0x10e(files, tstart, tend)

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

	% Make sure correct files were given.
	assert( min( strcmp(instr,   'fields') ) == 1, 'Only "fields" files are allowed.' );
	assert( min( strcmp(mode,    'hk') )     == 1, 'Only housekeeping (HK) files are allowed' );
	assert( min( strcmp(level,   'l1b') )    == 1, 'Only L1B files are allowed.' );
	assert( min( strcmp(optdesc, '10e') )    == 1, 'Only 10E files are allowed.' );

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
	afg_stemp_name = mms_construct_varname(sc, optdesc, 'afgstemp');
	dfg_stemp_name = mms_construct_varname(sc, optdesc, 'dfgstemp');
	afg_etemp_name = mms_construct_varname(sc, optdesc, 'afgetemp');
	dfg_etemp_name = mms_construct_varname(sc, optdesc, 'dfgetemp');

	% Read the magnetometer data
	[afg_stemp, tt2000] = MrCDF_nRead(files, afg_stemp_name, 'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	dfg_stemp           = MrCDF_nRead(files, dfg_stemp_name, 'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	afg_etemp           = MrCDF_nRead(files, afg_etemp_name, 'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	dfg_etemp           = MrCDF_nRead(files, dfg_etemp_name, 'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);

	% Save data structure.
	hk_10e = struct( 'tt2000',      tt2000,    ...
	                 'afg_stemp',   afg_stemp, ...
	                 'afg_etemp',   afg_etemp, ...
	                 'dfg_stemp',   dfg_stemp, ...
	                 'dfg_etemp',   dfg_etemp   );
end