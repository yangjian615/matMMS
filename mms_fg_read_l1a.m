%
% Name
%   mms_fg_read_l1a
%
% Purpose
%   Read L1A fluxgate magnetometer data from MMS.
%
% Calling Sequence
%   FG_L1A = mms_fg_read_l1a(FILES)
%     Read fluxgate level 1A data from files with file names FILES.
%
%   FG_L1A = mms_fg_read_l1a(FILES, TSTART, TEND)
%     Read fluxgate level 1A data from files with file names FILES
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
%   FG_L1A          out, required, type=struct
%                   Fields are:
%                       'tt2000'       - TT2000 epoch times for B_123
%                       'tt2000_ts'    - TT2000 epoch times for RANGE and SAMPLE_RATE
%                       'b_123'        - Magnetic field in sensor frame
%                       'range'        - Instrument range flag.
%                       'sample_rate'  - Sample rate
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-05-21      Written by Matthew Argall
%
function fg_l1a = mms_fg_read_l1a(files, tstart, tend)

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
	assert( min( strcmp(level, 'l1a') ) == 1, 'Only L1A files are allowed.' );

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
	b_name  = mms_construct_varname(sc, instr, '123');
	sr_name = mms_construct_varname(sc, instr, 'samplerate');
	if strcmp(instr, 'afg')
		range_name = mms_construct_varname(sc, instr, 'hirange');
	else
		range_name = mms_construct_varname(sc, instr, 'range');
	end

	% Read the magnetometer data
	[b_123, t]       = MrCDF_nRead(files, b_name,     'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	[range, t_range] = MrCDF_nRead(files, range_name, 'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);
	sample_rate      = MrCDF_nRead(files, sr_name,    'sTime', tstart, 'eTime', tend, 'ColumnMajor', true);

	% Transpose the data to be row vectors.
	fg_l1a = struct( 'tt2000',      t, ...
	                 'tt2000_ts',   t_range, ...
	                 'b_123',       b_123, ...
	                 'range',       range, ...
	                 'sample_rate', 1.0 ./ sample_rate );
end