%
% Name
%   mms_fg_read_l1b
%
% Purpose
%   Read L1B fluxgate magnetometer data from MMS.
%
% Calling Sequence
%   FG_L1B = mms_fg_read_l1a(FILES)
%     Read fluxgate level 1B data from files with file names FILES.
%
%   FG_L1B = mms_fg_read_l1a(FILES, TSTART, TEND)
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
%   FG_L1B          out, required, type=struct
%                   Fields are:
%                       'tt2000' - TT2000 epoch times
%                       'B_BCS'  - Magnetic field in BCS coordinates
%                       'B_OMB'  - Magnetic field in OMB coordinates
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-05-21      Written by Matthew Argall
%
function fg_l1b = mms_fg_read_l1b(files, tstart, tend)

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
	else
	assert( min( strcmp(level, 'l1b') ) == 1, 'Only L1B files are allowed.' );

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
	b_omb_name = mms_construct_varname(sc, instr, mode, 'omb');
	b_bcs_name = mms_construct_varname(sc, instr, mode, 'bcs');

	% Read the magnetometer data
	[b_omb, t] = MrCDF_nRead(files, b_omb_name, 'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	b_bcs      = MrCDF_nRead(files, b_bcs_name, 'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);

	% Transpose the data to be row vectors.
	fg_l1b = struct( 'tt2000', t,     ...
	                 'b_bcs',  b_bcs, ...
	                 'b_omb',  b_omb );
end