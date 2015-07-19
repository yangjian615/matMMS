%
% Name
%   mms_fdoa_read_defeph
%
% Purpose
%   Read MMS ephemeris data from FDOA ASCII files.
%
% Calling Sequence
%   EPHEMERIS = mms_fdoa_read_defatt(EPH_FILES, TSTART, TEND)
%     Return structure of ephemeris data EPHEMERIS, from MMS spacecraft SC
%     (e.g., 'mms2') between the times TSTART and TEND from ephemeris data
%     named EPH_FILES.
%
%   [EPHEMERIS, EPH_HDR] = mms_fdoa_read_defatt(__)
%     Also return a structure of header information from all files read.
%
% Parameters
%   EPH_FILES       in, required, type=char/cell
%   TSTART          in, required, type=char
%   TEND            in, required, type=char
%
% Returns
%   EPHEMERIS       out, required, type=struct
%                   Fields are:
%                     'UTC'       -  Time tags in UTC.
%                     'TAI'       -  Time tags in TAI as days since MMS reference epoch.
%                     'Position'  -  (x,y,z) components of spacecraft position.
%                     'Velocity'  -  (x,y,z) components of spacecraft velocity.
%                     'Mass'      -  Spacecraft mass.
%   EPH_HDR         out, optional, type=struct
%
% Examples
%   Given the files
%     MMS2_DEFEPH_2015078_2015079.V00
%     MMS2_DEFEPH_2015079_2015080.V00
%
%   Read all the data from
%     >> sc      = 'mms2';
%     >> tstart  = '2015-03-20T00:00:00Z';
%     >> tend    = '2015-03-21T00:00:00Z';
%     >> att_dir = '/Users/argall/Documents/Work/Data/MMS/Ephemeris/'
%     >> data    = mms_fdoa_read_defatt(sc, tstart, tend, att_dir)
%        data    =      UTC: {4019x1 cell}
%                       TAI: [4019x1 double]
%                  Position: [4019x3 double]
%                  Velocity: [4019x3 double]
%                      Mass: [4019x1 double]
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-12      Written by Matthew Argall
%   2015-07-16      Take file names as inputs. - MRA
%
function [ephemeris, eph_hdr] = mms_fdoa_read_defeph(eph_files, tstart, tend)

%------------------------------------%
% Check Files                        %
%------------------------------------%
	% Number of files given
	if iscell(eph_files)
		nFiles = length(eph_files);
		assert(nFiles > 0, 'At least one file name must be given.');
	else
		assert(ischar(eph_files) && isrow(eph_files), 'EPH_FILES must be a file name or cell array of file names.')
		nFiles = 1;
	end
	
	if nargin() < 3
		tend = '';
	end
	if nargin() < 2
		tstart = '';
	end

%------------------------------------%
% Read Header from Each File         %
%------------------------------------%
	% Read the header
	eph_hdr = mms_fdoa_read_ephem_header(eph_files{1});

	% Step through the rest of the files
	for ii = 2 : nFiles
		temp = mms_fdoa_read_ephem_header(eph_files{ii});
		
		% Append the header information
		eph_hdr.('DataStart')  = [ eph_hdr.('DataStart'),  temp.('DataStart')  ];
		eph_hdr.('Header')     = [ eph_hdr.('Header'),     temp.('Header')     ];
		eph_hdr.('StartTime')  = [ eph_hdr.('StartTime'),  temp.('StartTime')  ];
		eph_hdr.('StopTime')   = [ eph_hdr.('StopTime'),   temp.('StopTime')   ];
	end

%------------------------------------%
% Read the Data                      %
%------------------------------------%
	% Describe the file
	column_names = {'UTC', 'TAI', 'Position', 'Position', 'Position', ...
	                'Velocity', 'Velocity', 'Velocity', 'Mass'};
	column_types = {'char', 'double', 'double', 'double', 'double', ...
	                'double', 'double', 'double', 'double'};
	data_start = eph_hdr.('DataStart')(1);

	% Read the files
	ephemeris = MrFile_Read_nAscii(eph_files,                   ...
	                               'ColumnNames', column_names, ...
	                               'ColumnTypes', column_types, ...
	                               'DataStart',   data_start);
	
	% Convert TAI to TT2000
	tt2000 = mms_fdoa_epoch2tt2000(ephemeris.('TAI'), 'EphemTAI', true);

%------------------------------------%
% Time Interval                      %
%------------------------------------%
	trange = zeros(1, 2, 'int64');

	% Start of range
	if isempty(tstart)
		trange(1) = tt2000(1);
	else
		temp      = MrTimeParser(tstart, '%Y-%M-%dT%H:%m:%S', '%Y-%M-%dT%H:%m:%S.%1%2%3');
		trange(1) = spdfparsett2000(temp);
	end
	
	% End of range
	if isempty(tend)
		trange(2) = tt2000(end);
	else
		temp      = MrTimeParser(tend, '%Y-%M-%dT%H:%m:%S', '%Y-%M-%dT%H:%m:%S.%1%2%3');
		trange(2) = spdfparsett2000(temp);
	end
	
	% Indices to keep
	irange    = zeros(1, 2);
	irange(1) = find( tt2000 >= trange(1), 1, 'first' );
	irange(2) = find( tt2000 <= trange(2), 1, 'last' );
	
	% Number of data returned
	names   = fieldnames(ephemeris);
	nFields = length(names);

	% Get rid of unwanted data
	for ii = 1 : nFields
		ephemeris.( names{ii} ) = ephemeris.( names{ii} )(:, irange(1):irange(2));
	end
	
	% Add tt2000 times to the structure
	ephemeris.( 'tt2000' ) = tt2000(irange(1):irange(2));

%------------------------------------%
% Remove Duplicates                  %
%------------------------------------%
	%
	% Definitive ephemeris files overlap, meaning if more than
	% one file is found, it could have repeated data. Since
	% attitude data will often be interpolated to data sample
	% times, and MATLAB's interp1() function does not consider
	% repeated values to be strictly monotonically increasing,
	% we have to remove repeats.
	%
	if nFiles > 1
		[~, iuniq]         = unique(ephemeris.TAI);
		ephemeris.UTC      = ephemeris.UTC(      1, iuniq );
		ephemeris.TAI      = ephemeris.TAI(      1, iuniq );
		ephemeris.tt2000   = ephemeris.tt2000(   1, iuniq );
		ephemeris.position = ephemeris.Position( :, iuniq );
		ephemeris.velocity = ephemeris.Velocity( :, iuniq );
		ephemeris.mass     = ephemeris.Mass(     1, iuniq );
	end
end


%
% Name
%   mms_fdoa_read_ephem_header
%
% Purpose
%   Parse information for the header of MMS definitive ephemeris files.
%
% Calling Sequence
%   EPH_HDR = mms_fdoa_read_defatt_header(FILENAME)
%     Parse information from the header in file FILENAME and return a
%     structure EPH_HDR.
%
% Parameters
%   FILENAME        in, required, type=char
%
% Returns
%   EPH_HDR         out, required, type=structure
%                   A structure with the following fields:
%                     'DataStart'      Line at which data begins
%                     'Header'         Each line of the header as a cell array of strings.
%                     'StartTime'      First time in the file (converted to TT2000)
%                     'StopTime'       Last time in the file (converted to TT2000)
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-12      Written by Matthew Argall
%
function eph_hdr = mms_fdoa_read_ephem_header(filename)

	% File must exist.
	assert( exist(filename, 'file') == 2, ...
	        ['Ephemeris file does not exist: "' filename '".']);

	% Open the file
	fileID = fopen(filename);

%------------------------------------%
% Skip Over Header                   %
%------------------------------------%
	
	% Find the first empty line
	line       = 'temp';
	nCheck     = 100;
	count      = 0;
	data_start = -1;
	meta       = repmat({''}, 1, nCheck);
	
	% Step through each line until we get to the data
	while data_start < 0 && count <= nCheck
		% Read a line of data
		line        = fgetl(fileID);
		
		if ~isempty(regexp(line, 'StartTime', 'once'))
			parts      = regexp(line, '[ =()]+', 'split');
			start_time = mms_fdoa_epoch2tt2000( parts{end-1}, 'EphemTAI', true );
		
		elseif ~isempty(regexp(line, 'StopTime', 'once'))
			parts     = regexp(line, '[ =()]+', 'split');
			stop_time = mms_fdoa_epoch2tt2000( parts{end-1}, 'EphemTAI', true );
		
		elseif ~isempty(regexp(line, 'Spacecraft', 'once'))
			parts = regexp(line, '[ =()]+', 'split');
			sc    = parts{end};
		
		elseif ~isempty(regexp(line, '^[ \t]*(+|-)?[0-9]+', 'once'))
			data_start = count + 1;
		end
		
		% Count the line and save it.
		count       = count + 1;
		meta{count} = line;
	end
	
	% Trim off the last line of data.
	meta = meta(1:count-1);
	
	% Close the file
	fclose(fileID);
	
	% Create a structure
	eph_hdr = struct( 'DataStart', data_start, ...
	                  'Header',    { meta },   ...
	                  'StartTime', start_time, ...
	                  'StopTime',  stop_time );
end


% Definitive Orbit Ephemeris
% Spacecraft = MMS2
% StartTime = 2015-073/03:00:25.000 UTC  (MMS TAI: 20891.125694445)
% StopTime = 2015-074/05:28:55.000 UTC  (MMS TAI: 20892.228819445)
% MMS TAI Reference Epoch = 1958-001/00:00:00 UTC
% CentralBody = Earth
% ReferenceFrame = Mean of J2000
% PrincipalPlane = Equatorial
% Project = PPFA_ProductGeneration_DefinitiveEphem.m
% Source = FDGSS OPS
% FileCreationDate = Mar 16 2015 04:15:12.290 UTC
%  
% Epoch (UTC)             Epoch MMS TAI         X                       Y                       Z                       VX                      VY                      VZ                      Mass             
%                                               Km                      Km                      Km                      Km/Sec                  Km/Sec                  Km/Sec                  Kg               
% 2015-073/03:00:25.000   20891.125694445       8947.173812221          -8931.236682260         -6773.106446976         0.7986656112000         6.1909841475000         2.6588927471000         1353.6114694391  