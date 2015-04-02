%
% Name
%   c_find_filename
%
% Purpose
%   Create a valid Cluster file name.
%
% Calling Sequence
%   FILENAME = c_find_file(SC, TEAM, EXPERIMENT, DATE, STIME, ETIME)
%       Create a FILENAME out of the Cluster TEAM and EXPERIMENT name,
%       spacecraft number, SC, DATE of recorded data, with start and end
%       times, STIME and ETIME. The search is performed in the current
%       working directory.
%
%   FILENAME = c_find_file(..., DIRECTORY)
%       Search in a specific directory instead of pwd().
%
%   [__, STATUS] = c_find_file(__)
%       Returns the status of the search.
%
% Parameters
%   EXPERIMENT      in, required, type = char
%   SC              in, required, type = char/numeric
%   DATE            in, required, type = char
%   STIME           in, required, type = char
%   ETIME           in, required, type = char
%
% Returns:
%   FILENAME        out, required, type = char
%   STATUS          out, optional, type = integer
%                   Indicates output status:
%                      0-10: Success
%                        0 - An exact match was found.
%                        1 - A match wsa found, but with the start
%                            time > sTime and/or the end time >= eTime
%                     11-20: Fail
%                       11 - No files found for the given date.
%                       12 - Files for the same date were found, but none
%                            contain the requested time interval.
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2014-10-26      Written by Matthew Argall
%   2015-03-07      Added the TEAM parameter. Arranged as SC, TEAM, EXPERIMENT. - MRA
%
function [filename, status] = mms_find_file(sc, instr, mode, level, tstart, version, desc, directory)
  
	% Check number of arguments
	if nargin < 8
		directory = pwd();
	end
	
	if nargin < 7
		desc = '';
	end
	
	if nargin < 6
		version = '*';
	end

	% Create a file name
	test_file = mms_construct_filename(sc, instr, mode, level, tstart, version, desc);

%----------------------------------%
% Exact Match \\\\\\\\\\\\\\\\\\\\ %
%----------------------------------%
	if exist(fullfile(directory, test_file), 'file') == 2
		filename = test_file;
		status   = 1;

%----------------------------------%
% Find Partial Match \\\\\\\\\\\\\ %
%----------------------------------%
	else
		% Create a file without time information
		test_file = mms_construct_filename(sc, instr, mode, level);
		
		% Search for the file
		result_file = dir(fullfile(directory, test_file));
		nResults    = length(result_file);
		
	%----------------------------------%
	% No Partial Matches \\\\\\\\\\\\\ %
	%----------------------------------%
		if nResults == 0
			status   = 11;
			filename = '';
			
	%----------------------------------%
	% Select Partial Match \\\\\\\\\\\ %
	%----------------------------------%
		else
			% Take only the file name field
			if nResults == 1
				result_file = result_file.name;
			else
				result_file = { result_file(:).name };
			end
			
			% Get the start and end times of the file found.
			[~, ~, ~, ~, file_tstart] = mms_dissect_filename(result_file);
			
			% Concatenate into character arrays
			if nResults > 1
				file_tstart = vertcat( file_tstart{:} );
			end
			
			% Compare to original time
			file_tstart = sort( int64( str2num( file_tstart ) ) );
			iMatch = find( file_tstart <= int64( str2double( tstart ) ), 1, 'last');
			
			% No partial matches for this time interval
			if isempty(iMatch)
				status   = 12;
				filename = '';
				
			% Partial match found
			else
				status   = 1;
				filename = result_file{iMatch};
			end
		end
		
		% Append the directory
		if status < 11 && ~isempty(directory)
			filename = fullfile(directory, filename);
		else
			disp(['No match found "' fullfile(directory, test_file ) '".']);
		end
	end
end