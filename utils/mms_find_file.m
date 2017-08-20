%
% Construct an MMS file name. The file name format is:
%
%   scId_instrumentId_mode_dataLevel_optionalDataProductDescriptor_startTime_vX.Y.Z.cdf
%
% :Examples:
%   Dissect and construct a filename
%     filename = 'mms3_dfg_hr_l2p_duration-1h1m_20010704_v1.0.0.cdf';
%     [sc, instr, mode, level, desc, tstart, version] = mms_dissect_filename(filenames);
%     fname = mms_construct_filename(sc, instr, mode, level, desc, tstart, version)
%     fname =
%              mms3_dfg_hr_l2p_duration-1h1m_20010704_v1.0.0.cdf
%
% Calling Sequence:
%   FNAME = mms_construct_filename(SC, INSTR, MODE, LEVEL)
%     Construct an MMS file name FNAME using the spacecraft ID, SC,
%     instrumnet identifier INSTR, telemetry mode MODE and data level
%     LEVEL. Other mandatory components of the file name will be the
%     wildcard character '*' and optional components will be the empty
%     string ''.
%
%   FNAME = mms_construct_filename(..., 'ParamName', ParamValue)
%     Any parameter name and value described below.
%
% Parameters
%   SC:                 in, required, type=char
%                       mms1, mms2, mms3, mms4
%   INSTR:              in, required, type=char
%                       Instrument or investigation identifier
%                           hpca
%                           aspoc
%                           epd
%                           epd-eis
%                           epd-feeps
%                           fpi
%                           des
%                           dis
%                           des-dis
%                           fields
%                           edi
%                           adp
%                           sdp
%                           adp-sdp
%                           afg
%                           dfg
%                           dsp
%                           afg-dfg
%                           scm
%   MODE:               in, required, type=string
%                       Instrument telemetry mode:
%                           fast
%                           slow
%                           brst
%                           srvy
%   LEVEL:              in, required, type=string
%                       Level of the data product:
%                           l1a
%                           l1b
%                           l2
%                           ql
%                           l2pre
%                           l2plus
%   'Directory':        in, required, type=char
%                       A directory to be appended to the file name.
%   'OptDesc':          in, optional, type=char, default=''
%                       Optional data product descriptor. Should be short
%                           (3-8 characters). Hyphens used to separate
%                           multiple components.
%   'RelaxedTStart':    in, required, type=boolean, default=false
%                       Since MMS filenames contain only 'TStart', matching
%                           files to the time interval provided is not
%                           perfect. Normally, the first file found will be
%                           the closest matching file that starts /before/
%                           'TStart'. If that is the case, the default is
%                           to exclude if from the search results. This
%                           behavior is not ideal, however, for DFG calibration
%                           files, which start on 2010-10-01. Set this
%                           parameter to true to keep the first file, regardless
%                           of its start time.
%   'SDCroot':          in, required, type=char, default='/nfs/'
%                       If your directory structure models that of the SDC,
%                           then provide its root directory.
%   'SubDirs':          in, required, type=cell, default={'%Y', '%M'}
%                       Date-time subdirectories, specified as MrTokens tokens
%                           (e.g. %Y = year, %M = month, %d = day). Normally,
%                           the SDC directory chain ends with /%Y/%M (/year/month/).
%                           If that is not the case, use 'SubDirs' to specify the
%                           level of subdirectories. The empty string, '',
%                           indicates no date-time subdirectories.
%   'TStart':           in, required, type=char
%                       Start time of the data product, formatted as:
%                           'yyyy-mm-ddTHH:MM:SS'. Least significant fields can
%                           be dropped when files start on regular hourly
%                           or minute boundaries.
%   'TEnd':             in, required, type=char
%                       Start time of the data product, formatted as:
%                           'yyyy-mm-ddTHH:MM:SS'. Least significant fields can
%                           be dropped when files start on regular hourly
%                           or minute boundaries.
%   'TimeOrder':        in, required, type=char, default='%Y%M%d'
%                       In order to compare 'TStart' and 'TEnd' with the start
%                           time of each file, the file's start time must be
%                           in a predictable order. Specify the date-time format
%                           of the file names using MrTokens tokens. Most files
%                           have [year month day], which is the default.
%   'Version':          in, required, type=char
%                       Version number in the form: "vX.Y.Z"
%                           X - Interface number. Increments represent
%                               significant changes that will break code or
%                               require code changes in analysis software.
%                           Y - Quality number. Represents change in
%                               quality of the, such as calibration or
%                               fidelity. Should not impact software.
%                           Z - Bug fix/Revision number. Minor changes to
%                               the contents of the file due to
%                               reprocessing of missing data. Dependent
%                               data products should be reprocessed.
%
function [files, nFiles, searchstr] = mms_find_file(sc, instr, mode, level, varargin)

%------------------------------------%
% Inputs                             %
%------------------------------------%
	tstart         = '';
	tend           = '';
	optdesc        = '';
	version        = '';
	subdirs        = '';
	timeorder      = '';
	directory      = '';
	dropbox        = '';
	sdc_root       = '/nfs/';
	relaxed_tstart = false;
	nOptArgs       = length(varargin);

	% Optional parameters
	for ii = 1 : 2 : nOptArgs
		switch varargin{ii}
			case 'Directory'
				directory = varargin{ii+1};
			case 'Dropbox'
				dropbox = varargin{ii+1};
			case 'OptDesc'
				optdesc = varargin{ii+1};
			case 'RelaxedTStart'
				relaxed_tstart = varargin{ii+1};
			case 'SDCroot'
				sdc_root = varargin{ii+1};
			case 'SubDirs'
				subdirs = varargin{ii+1};
			case 'TStart'
				tstart = varargin{ii+1};
			case 'TEnd'
				tend = varargin{ii+1};
			case 'TimeOrder'
				timeorder = varargin{ii+1};
			case 'Version'
				version = varargin{ii+1};
			otherwise
				error([ 'Parameter name not recognized: "' varargin{ii} '"' ]);
		end
	end
	
	% Time order
	if strcmp(timeorder, '')
		if strcmp(mode, 'brst')
			timeorder = '%Y%M%d%H%m%S';
		else
			timeorder = '%Y%M%d';
		end
	end
	
	% Check time formats
	if ~isempty(tstart)
		assert( MrTokens_IsMatch(tstart, '%Y-%M-%dT%H:%m:%S'), ...
		        'TStart must be an ISO-8601 string: "yyyy-mm-ddThh:mm:ss".' );
	end
	if ~isempty(tend)
		assert( MrTokens_IsMatch(tend, '%Y-%M-%dT%H:%m:%S'), ...
		        'TEnd must be an ISO-8601 string: "yyyy-mm-ddThh:mm:ss".' );
	end
	
	% Subdirectories
	if isempty(subdirs)
		if strcmp(mode, 'brst')
			subdirs = {'%Y', '%M', '%d'};
		else
			subdirs = {'%Y', '%M'};
		end
	else
		% Turn char to cell for ease of use.
		assert( ischar(subdirs) || iscell(subdirs), 'SUBDIRS must be a string or cell array of strings.' );
		if ischar(subdirs)
			subdirs = { subdirs };
		end
	end

%------------------------------------%
% Directory                          %
%------------------------------------%
	if ~isempty(directory)
		sdc_dir = directory;
	else
		% First half of directory
		sdc_dir  = fullfile(sdc_root, sc, instr, mode, level, optdesc);
		
		% Test the directory
%		assert( exist(sdc_dir, 'dir') == 7, ['SDC directory does not exist: "' sdc_dir '".'] );
		
		% Date subdirectories
		sdc_dir = fullfile( sdc_dir, subdirs{:} );
	end

%------------------------------------%
% Find File in Database              %
%------------------------------------%
	% Create the file name
	%   - Use the SDC-like structure
	searchstr = mms_construct_filename(sc, instr, mode, level, ...
	                                   'Directory', sdc_dir, ...
	                                   'OptDesc',   optdesc, ...
	                                   'Version',   version, ...
	                                   'TStart',    timeorder);

	% Search for the file.
	[files, nFiles] = MrFile_Search(searchstr,               ...
	                                'TStart',    tstart,    ...
	                                'TEnd',      tend,      ...
	                                'TimeOrder', timeorder, ...
	                                'Closest',   true);

%------------------------------------%
% Find File in Dropbox               %
%------------------------------------%
	if ~isempty(dropbox)
		% Reformat the file name
		[~, base, ext] = fileparts(searchstr);
		searchstr      = fullfile(dropbox, [base ext]);

		% Search in dropbox
		[files_db, nDropbox] = MrFile_Search(searchstr,              ...
		                                     'TStart',    tstart,    ...
		                                     'TEnd',      tend,      ...
		                                     'TimeOrder', timeorder, ...
		                                     'Closest',   true);

		% Append to initial results
		files = [ files, files_db ];
	end

%------------------------------------%
% Refine Results                     %
%------------------------------------%

	%
	% Because MMS file names contain a start time, but no end time,
	% the first file returned by MrFile_Search may not lie between
	% TSTART and TEND. Here, we compare the date within the file name
	% and ensure that it is within the same day as what was given
	%
	% This is only done if more than one files are found. One reason
	% is that the FGM Calibration files start on 2010-10-01 and there
	% is a single file for the entire mission. Therefore, I assume
	% that if one file is found, it contains data that extends into
	% your time period of interest.
	%
	if nFiles > 0 && ~isempty(tstart) && ~relaxed_tstart
		% Get the start time of the first file.
		if nFiles == 1
			[~, ~, ~, ~, fstart] = mms_dissect_filename( files );
		else
			[~, ~, ~, ~, fstart] = mms_dissect_filename( files{1} );
		end

		% If the files's start time does not match the requested tstart time
		if ~strcmp( fstart(1:4), tstart(1:4) ) || ...
		   ~strcmp( fstart(5:6), tstart(6:7) ) || ...
		   ~strcmp( fstart(7:8), tstart(9:10) )
		   
			% Remove the first file.
			nFiles = nFiles - 1;
			if nFiles > 1
				files = files(2:end);
			else
				files = '';
			end
		end
	end
end