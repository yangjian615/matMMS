%
% Name
%   mms_sdc_proc_edi_ql_efield
%
% Purpose
%   Find all EDI L1A electric field-mode data files and calculate the electric field.
%
% Calling Sequence
%   FILES = mms_sdc_proc_edi_ql_efield(SC, TSTART, TEND)
%     Read and process all data from MMS spacecraft SC between the times
%     TSTART and TEND. Data is saved to CDF files named FILES.
%
%   FILES = mms_sdc_proc_edi_ql_efield()
%     Process all data from all spacecraft from the beginning of the mission
%     until three days before the current date.
%
% Parameters
%   SC:             in, optional, type=char/cell, default={ 'mms1', 'mms2', 'mms3', 'mms4' }
%   TSTART:         in, optional, type=char, default='2015-03-11T00:00:00'
%   TEND:           in, optional, type=char, default=three days before today's date
%
% Returns
%   FILES           out, optional, type=cell
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-06-03      Written by Matthew Argall
%   2015-08-07      Generate dates instead of looking for files then extracting dates. - MRA
%
function files = mms_sdc_proc_edi_ql_efield(sc, tstart, tend)

%------------------------------------%
% Find Files                         %
%------------------------------------%
	% EDI L1A E-Field Data Files
	if nargin < 3
		tend = datestr( now() - 3.0, 'yyyy-mm-dd' );
	end
	if nargin < 2
		tstart = '2015-03-17T00:00:00';
	end
	if nargin < 1
		sc = { 'mms1' 'mms2' 'mms3' 'mms4' };
	end

%------------------------------------%
% Process Data for Each Spacecraft   %
%------------------------------------%
	% Make sure we have a cell.
	if ischar(sc)
		sc = { sc };
	end
	
	% Convert dates to day-of-year
	dates  = MrDateGen(tstart(1:10), tend(1:10));
	fstart = strcat( dates, 'T00:00:00' );
	fend   = strcat( dates, 'T24:00:00' );

	nSC    = length(sc);
	nDates = length(dates);
	files  = cell( nSC * nDates );
	count  = 0;

	% Step through each spacecraft
	for ii = 1 : length(sc)
		% Loop through each date
		for jj = 1 : length(dates)
			% Create the data
			try
				files{count+1} = mms_edi_create_ql_efield(sc{ii}, fstart{jj}, fend{jj});
				count          = count + 1;
			catch ME
				% Print error
				fprintf('Unable to create file: %s %s %s\n', sc{ii}, fstart{jj}, fend{jj});
				fprintf('  Error using %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
				fprintf('  %s\n', ME.message);
				
				% Print stack
				for ii = 2: length(ME.stack)
					fprintf('      %s (line %d)\n', ME.stack(ii).name, ME.stack(ii).line);
				end
				fprintf('\n');
			end
		end
	end
	
	% Trim files
	files = files(1:count);
end