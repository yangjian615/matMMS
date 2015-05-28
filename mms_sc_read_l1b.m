%
% Name
%   mms_sc_read_l1b
%
% Purpose
%   Read L1B search coil magnetometer data from MMS.
%
% Calling Sequence
%   SC_L1B = mms_sc_read_l1b(FILES)
%     Read search coil level 1B data from files with file names FILES.
%
%   SC_L1B = mms_sc_read_l1b(FILES, TSTART, TEND)
%     Read search coil level 1A data from files with file names FILES
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
%   SC_L1B          out, required, type=struct
%                   Fields are:
%                       'tt2000'        - TT2000 epoch times for B_OMB
%                       'b_omb'         - Magnetic field in OMB frame
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-05-21      Written by Matthew Argall
%
function sc_l1b = mms_sc_read_l1b(files, tstart, tend)

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
	assert( min( strcmp(instr, 'scm') ) == 1, 'Only SCM files are allowed.' );

	% Instr, Level, Mode
	if nFiles > 1
		assert( min( strcmp(mode, mode{1}) ) == 1, 'All files must have the same telemetry mode.' );
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
	b_name = mms_construct_varname(sc, instr, optdesc, 'scm123');

	% Read the magnetometer data
	[b_omb, t] = MrCDF_nRead(files, b_name,  'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);

	% Transpose the data to be row vectors.
	sc_l1b = struct( 'tt2000',        t, ...
	                 'b_omb',         b_omb );
end