%
% Name
%   mms_file_search
%
% Purpose
%   Find MMS data files.
%
% Calling Sequence
%   FILES = mms_file_search(SC, INSTR, MODE, LEVEL)
%     Find all files in the current directory with spacecraft id SC (e.g.,
%     'mms1') instrument id INSTR (e.g. 'scm'), data rate mode MODE (e.g.
%     'sc128'), and data product level LEVEL (e.g. 'l2pre').
%
%   FILES = mms_file_search(__, 'ParamName', ParamValue)
%     Filter results using any of the parameter name-value pairs below.
%
% Parameters
%   SC              in, required, type = char
%   INSTR           in, required, type = char
%   MODE            in, required, type = char
%   LEVEL           in, required, type = char
%   'Closest'       in, optional, type = boolean, default = false
%                   Find the file that begins closest to TStart. This
%                     option is ignored unless TStart is specified.
%   'Directory'     in, optional, type = char, default = pwd()
%                   Look in this directory instead of the present working
%                     directory. The directory path may include tokens
%                     recognized by MrTokens.m
%   'Newest'        in, optional, type = boolean
%                   Return only the newest version of each file found.
%                     Unless 'Version' is given, this is the default
%                     behavior.
%   'OptDesc'       in, optional, type = char, default = '*'
%                   Optional descriptor included in the file names.
%   'Version'       in, optional, type = char, default = ''
%                   Return a specific version of a files.
%   'TStart'        in, optional, type = char, default = ''
%                   An ISO-8601 string. All files that end before this date
%                     are excluded from the results. All files that start
%                     on or after this date are included.
%   'TEnd'          in, optional, type = char, default = ''
%                   An ISO-8601 string. All files that start on or after
%                     this date are excluded from the results. All files
%                     that start before this date are included.
%
% Returns:
%   FILES           out, required, type = cell
%
% Examples
%   Given the directory structure
%     >> directory = '/Users/argall/Documents/Work/Data/MMS/SCM';
%
%   with contents
%     >> ls(directory)
%       mms2_scm_comm_l1a_sc128_20150317_v0.3.0.cdf
%       mms2_scm_comm_l1a_sc128_20150318_v0.0.0.cdf
%       mms2_scm_comm_l1a_sc128_20150318_v0.2.0.cdf
%       mms2_scm_comm_l1a_sc128_20150318_v0.7.0.cdf
%       mms2_scm_comm_l1a_sc128_20150319_v0.3.0.cdf
%       mms2_scm_comm_l1a_sc128_20150320_v0.3.0.cdf
%
%   1. Find the newest verion of each file
%     >> files = mms_file_search('mms2', 'scm', 'comm', 'l1a', ...
%                                'Directory', directory,       ...
% 															 'OptDesc',   'sc128');
%     >> vertcat(files{:})
%       mms2_scm_comm_l1a_sc128_20150317_v0.3.0.cdf
%       mms2_scm_comm_l1a_sc128_20150318_v0.7.0.cdf
%       mms2_scm_comm_l1a_sc128_20150319_v0.3.0.cdf
%       mms2_scm_comm_l1a_sc128_20150320_v0.3.0.cdf
%
%   2. Find files that begin at or after 2015-03-19T00:00:00Z
%     >> files = mms_file_search('mms2', 'scm', 'comm', 'l1a', ...
%                                'Directory', directory,       ...
% 															 'OptDesc',   'sc128',         ...
%                                'TStart',    '2015-03-19T00:00:00Z');
%     >> vertcat(files{:})
%       mms2_scm_comm_l1a_sc128_20150319_v0.3.0.cdf
%       mms2_scm_comm_l1a_sc128_20150320_v0.3.0.cdf
%
%   3. Also exclude files that start at or after 2015-03-20T00:00:00Z
%     >> files = mms_file_search('mms2', 'scm', 'comm', 'l1a',        ...
%                                'Directory', directory,              ...
% 															 'OptDesc',   'sc128',                ...
%                                'TStart',    '2015-03-19T00:00:00Z', ...
%                                'TEnd',      '2015-03-20T00:00:00Z');
%     >> files{:}
%       mms2_scm_comm_l1a_sc128_20150319_v0.3.0.cdf
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-03      Written by Matthew Argall
%
function [files] = mms_file_search(sc, instr, mode, level, varargin)

%------------------------------------%
% Inputs                             %
%------------------------------------%
	
	% Defaults
	closest   = false;
	newest    = false;
	directory = '';
	tstart    = '';
	tend      = '';
	optDesc   = '';
	version   = '';

	% Check for optional arguments
	nOptArgs = length(varargin);
	for ii = 1 : 2: nOptArgs
		switch varargin{ii}
			case 'Closest'
				closest = varargin{ii+1};
			case 'Newest'
				newest = varargin{ii+1};
			case 'Directory'
				directory = varargin{ii+1};
			case 'OptDesc'
				optDesc = varargin{ii+1};
			case 'Version'
				version = varargin{ii+1};
			case 'TStart'
				tstart = varargin{ii+1};
			case 'TEnd'
				tend = varargin{ii+1};
			otherwise
				error( ['Unknown parameter "' varargin{ii} '".'] );
		end
	end
	
	% Look in the current directory
	if isempty(directory)
		directory = pwd();
	end
	
	% If a specific version is not specified, choose the newest
	if isempty(version)
		newest = true;
	end
	
	% Newest and Version are mutually exclusive.
	%   - One or both must not be set.
	assert( isempty(version) || ~newest,  'Version and Newest are mutually exclusive.' );
	assert( isempty(tend)    || ~closest, 'Closest and TEnd are mutually exclusive.' );

%------------------------------------%
% Find Files                         %
%------------------------------------%
	
	% Create a file name
	%   - Use tokens so that all times are found.
	%   - Do not specify version so all versions are found.
	fname = mms_construct_filename(sc, instr, mode, level, ...
		                             'Tokens',    true, ...
																 'OptDesc',   optDesc, ...
																 'Directory', directory);
	
	% Search for the files.
	allFiles = MrFile_Search(fname);

%------------------------------------%
% Find Copies                        %
%------------------------------------%

	%
	% Copies differ only in their version numbers.
	%
	
	% Breakdown files into parts
	[allSC, allInstr, allMode, allLevel, allTStart, allVersion, allDesc] = mms_dissect_filename(allFiles);
	
	% Separate the optional description from TStart
	%   - But only if it is not empty
	iDesc          = ~cellfun(@isempty, allDesc);
	allDesc(iDesc) = strcat( allDesc(iDesc), {'_'} );

	% Build up the file names again, but without version numbers.
	%   - TODO: Remove extra '_' when desc = ''
	base = strcat( allSC, {'_'}, allInstr, {'_'}, allMode, {'_'}, allLevel, {'_'}, allDesc, allTStart );
	
	% Sort the results
	%   - BASE_UNIQUE(IUNIQ) is the same as BASE, so IUNIQ maps the
	%     elements of BASE onto the unique elements in BASE_UNIQ. So, to look
	%     for copies in BASE, we look for repeated indices in IUNIQ.
	%   - IBASE is the unique elements of BASE, so contains the maximum
	%     number of elements we must search through.
	[files, ibase, iuniq] = unique(base);
	nFiles                = length(files);
	
	% Now append the version and file extension onto the base
	base = strcat(base, {'_'}, allVersion, {'.cdf'});

%------------------------------------%
% Return Newest                      %
%------------------------------------%
	if newest
		% Search through all of the uniq elements
		for ii = 1 : length(ibase)
			
			% Pick out copies of the current file
			iCopies = find(ibase(iuniq) == ibase(ii));
			nCopies = length(iCopies);
			
			% We only need to compare if nCopies > 1
			newestVersion = base{ iCopies(1) };
			
			% Step through all copies
			for jj = 2 : nCopies

				% Compare copies and keep the newest version.
				if mms_comp_version( base{ iCopies(jj) }, newestVersion ) == 1
					newestVersion = base{ iCopies(jj) };
				end
			end
			
			% Keep the newest
			files{ii} = newestVersion;
		end
	
%------------------------------------%
% Specific Version                   %
%------------------------------------%
	elseif ~isempty(version)
		% Search through all of the uniq elements
		for ii = 1 : length(ibase)
			
			% Pick out copies of the current file
			iCopies = find(ibase(iuniq) == ibase(ii));
			nCopies = length(iCopies);
			
			% We only need to compare if nCopies > 1
			thisVersion = [];
			
			% Step through copies until we find a match.
			jj = 2;
			while isempty(thisVersion) && jj <= nCopies

				% Compare copies and keep a matching version.
				if mms_comp_version( base{ iCopies(jj) }, version ) == 0
					thisVersion = base{ iCopies(jj) };
				end
			end
			
			% Keep the matching version
			files{ii} = thisVersion;
		end
	end

%------------------------------------%
% Specific Time                      %
%------------------------------------%
	if ~isempty(tstart) || ~isempty(tend)
		% Extract times from file names
		[~, ~, ~, ~, fileStart] = mms_dissect_filename(files);
		
		% Make sure FILESTART has all 14 characters: YYYYMMddHHmmSS
		%   - Append '0' until it does.
		for ii = 1 : nFiles
			nAddZeros     = 14 - length( fileStart{ii} );
			zeroStr       = repmat( '0', 1, nAddZeros );
			fileStart(ii) = strcat( fileStart(ii), zeroStr );
		end
			
		% Find files that start at or after TSTART
		%   - Times are in descending order, so can be compared directly
		fileStart = strcat( {'int64('}, fileStart, {')'} );
		fileStart = cellfun(@str2num, fileStart);
		
		% Create a regular expression to dissect TSTART and TEND
		tpattern = '%Y-%M-%dT%H:%m:%S';
		regex    = MrTokens_ToRegex(tpattern);
				
	%------------------------------------%
	% Start Time                         %
	%------------------------------------%
		if ~isempty(tstart)
			% Dissect TSTART and assemble it without delimiters
			tparts = regexp(tstart, regex, 'tokens');
			tparts = strcat(tparts{:});
			tstart = strcat(tparts{:});
			
			%Convert file times to int64
			tstart = str2num( ['int64(' tstart ')'] );
			
			% Find all files such that FILESTART >= TSTART
			files     = files( fileStart >= tstart );
			fileStart = fileStart( fileStart >= tstart );
				
		%------------------------------------%
		% Closest Time                       %
		%------------------------------------%
			if closest && length(files) > 1
				% Find the file that starts closest to TSTART
				[~, iClosest] = min( fileStart - tstart );
				
				% Select only that one file
				files     = files( iClosest );
				fileStart = fileStart( iClosest );
			end
		end
		
	%------------------------------------%
	% End Time                           %
	%------------------------------------%
		if ~isempty(tend)
			% Dissect TEND and assemble it without delimiters
			tparts = regexp(tend, regex, 'tokens');
			tparts = strcat(tparts{:});
			tend   = strcat(tparts{:});
			
			%Convert file times to int64
			tend = str2num( ['int64(' tend ')'] );
			
			% Find all files such that FILESTART < TEND
			files     = files( fileStart < tend );
			fileStart = fileStart( fileStart < tstart );
		end
	end
		
	%------------------------------------%
	% Directory                          %
	%------------------------------------%
	nFiles      = length(files);
	dir_cell    = cell(1, nFiles);
	dir_cell(:) = { directory };
	files       = cellfun(@fullfile, dir_cell, files, 'UniformOutput', false);
end