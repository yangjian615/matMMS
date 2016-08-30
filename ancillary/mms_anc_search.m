%
% Name
%   mms_anc_search
%
% Purpose
%   Search for MMS data files.
%
%   At the SDC, data files are initially saved in a dropbox location before
%   being mored to their final resting place in the directory structure. This
%   program provide a mechanism for searching in one folder or the other, or
%   both.
%
% Calling Sequence
%   FILE = mms_anc_search(FOLDER, SC, INSTR, TSTART);
%     Search for a file in directory DIR with spacecraft identifier SC,
%     instrument identifier INSTR, telemetry mode MODE, data level LEVEL,
%     and data start time TSTART. File names of all versions are returned
%     in FILE.
%
%   [__] = mms_anc_search(..., 'ParamName', ParamValue);
%     Also supply any of the parameter name-value pairs listed below.
%
% Parameters
%   FOLDER:         in, required, type=char
%   SC:             in, required, type=char
%   INSTR:          in, required, type=char
%   TSTART:         in, required, type=char
%   'Version':      in, optional, type=char, default='*'
%                   File version number, formatted as 'X.Y.Z', where X, Y, and
%                     Z are integers.
%   'RootDir':      in, optional, type=logical/char, default='*'
%                   If boolean and true, then `FOLDER` is taken to be the root of
%                     an SDC-like directory structure. If char, then "RootDir"
%                     is itself taken to be the root of an SDC-like directory
%                     structure. In the latter case, results from both ROOTDIR
%                     and FOLDER are returned.
%
% MATLAB release(s) MATLAB 8.2.0.701 (R2013b)
% Required Products None
%
% History:
%   2016-04-01      Written by Matthew Argall
%
function files = mms_anc_search(folder, sc, instr, tstart, varargin)

%------------------------------------%
% Inputs                             %
%------------------------------------%

	% Defaults
	version  = '*';
	root_dir = false;
	files    = {};

	% Check optional parameters
	nOptArgs = length(varargin);
	for ii = 1 : 2 : nOptArgs
		switch varargin{ii}
			case 'Version'
				version = varargin{ii+1};
			case 'RootDir'
				root_dir = varargin{ii+1};
			otherwise
				error('MMS:badOptArg', 'Optional argument not recognized: "%s".', varargin{ii})
		end
	end
	
	% Root directory?
	if islogical(root_dir) && root_dir == true
		tf_root   = true;
		data_path = folder;
	elseif ~isempty(root_dir) && ischar(root_dir)
		tf_root   = true;
		data_path = root_dir;
	else
		tf_root   = false;
		data_path = '';
	end

%------------------------------------%
% Create File Name                   %
%------------------------------------%
	
	% Convert TSTART to day-of-year
	tvec = mms_parse_time(tstart);
	tvec = cellfun( @str2num, tvec );
	doy  = datenum(tvec(1:3)) - datenum([tvec(1), 1, 1]);
	
	% Create the file name
	filename = sprintf('%s_%s_*%04i%03i*.V%s', upper(sc), upper(instr), tvec(1), doy, version);

%------------------------------------%
% Search in FOLDER                   %
%------------------------------------%

	% Search in FOLDER
	%   - TF_ROOT == False     ==> Search only in FOLDER
	%   - FOLER   ~= DATA_PATH ==> Search in both FOLDER and ROOTDIR
	if ~tf_root || ~strcmp(data_path, folder)
		% Search for files
		files = dir( fullfile( folder, filename ) );
		
		% Append directory to file name
		%   - FullFile accepts only single file and dir names ==> CellFun
		%   - Make a copy of FOLDER for each file found
		%   - Concatenate
		if isempty(files)
			files = '';
		else
			files = strcat( folder, filesep, { files.name } );
%			files = cellfun(@fullfile, repmat( {folder}, 1, length(files) ), { files.name }, 'UniformOutput', false);
		end
	end
	
%------------------------------------%
% Search in ROOT                     %
%------------------------------------%

	% Search in RootDir
	if tf_root
		% Create the path
		data_path = fullfile(data_path, 'ancillary', sc, instr);
		
		% Search for files
		files2 = dir( fullfile( data_path, filename ) );
		if ~isempty(files2)
			files = [ files, strcat( data_path, filesep, { files2.name } ) ] ;
		end
	end

%------------------------------------%
% Find Latest Version                %
%------------------------------------%

	if ~isempty(files)
		% Remove the directory and version number (the extension)
		[~, fnames, ext] = cellfun(@fileparts, files, 'UniformOutput', false );
		
		% Extract the version numbers
		v = regexp(ext, '\.V([0-9]+)', 'tokens');
		v = vertcat( v{:} );
		v = cellfun( @str2num, vertcat( v{:} ) );
	
		% Find uniq files
		[funiq, iuniq] = unique(fnames);
		nuniq          = length(iuniq);
		f_out          = cell( 1, nuniq );
		
		% Step through each unique file name
		for ii = 1 : nuniq
			% Files of interest
			tf_files    = strcmp( fnames, fnames(iuniq(ii)) );
			these_files = files( tf_files );
			these_v     = v( tf_files );
		
			% Find the largest version number
			[~, iv]   = max(these_v);
			f_out(ii) = these_files(iv);
		end
		
		% Transfer results to output
		files = f_out;
	end

	% Return single string
	if isempty(files)
		files = '';
	elseif length(files) == 1
		files = files{1};
	end
end