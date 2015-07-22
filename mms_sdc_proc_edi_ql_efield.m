%
% Name
%   mms_sdc_proc_edi_ql_efield
%
% Purpose
%   Find all EDI L1A electric field-mode data files and calculate the electric field.
%
% Calling Sequence
%   SAVE_FILE = mms_sdc_bestarg()
%     Read and process all data required by Bestarg and save them to a
%     mat file named SAVE_FILE.
%
% Parameters
%
% Returns
%   SAVE_FILE       out, required, type=string
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-06-03      Written by Matthew Argall
%
function [] = mms_sdc_proc_edi_ql_efield()

%------------------------------------%
% Find Files                         %
%------------------------------------%
	% EDI L1A E-Field Data Files
	sc       = { 'mms1' 'mms2' 'mms3' 'mms4' };
	instr    = 'edi';
	mode     = 'slow';
	level    = 'l1a';
	optdesc  = 'efield';
	tstart   = '2015-03-01T00:00:00';
	tend     = '2015-06-18T24:00:00';

%------------------------------------%
% Process Data for Each Spacecraft   %
%------------------------------------%
	% Make sure we have a cell.
	if ischar(sc)
		sc = { sc };
	end

	% Step through each spacecraft
	for ii = 1 : length(sc)
		% Search for files.
		[edi_file, count, str] = mms_file_search(sc{ii}, instr, mode, level, ...
		                                         'TStart',    tstart, ...
		                                         'TEnd',      tend, ...
		                                         'OptDesc',   optdesc);
		if count == 0
			warning('SDC:EDI', ['EDI file not found: "' str '".']);
			continue
		end

	%------------------------------------%
	% Loop Through Each File             %
	%------------------------------------%
		for jj = 1 : count
		
			% Dissect the file name to discover which days efield data exists
			[~, ~, ~, ~, fstart] = mms_dissect_filename(edi_file{jj});
	
			% Convert time to yyyy-mm-ddTHH:MM:SS
			date_temp = MrTimeParser(fstart, '%Y%M%d', '%Y-%M-%d');
			fstart = strcat(date_temp, 'T00:00:00');
			fend   = strcat(date_temp, 'T24:00:00');
			
			% Create the data
			try
				cdf_file = mms_edi_create_ql_efield(sc{ii}, fstart, fend);
			catch ME
				fprintf('Unable to create file: %s %s %s\n', sc{ii}, fstart, fend);
			end
		end
	end