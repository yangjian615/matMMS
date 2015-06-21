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
	sc       = 'mms%([1-4]%)';
	instr    = 'edi';
	mode     = 'slow';
	level    = 'l1a';
	optdesc  = 'efield';
	sdc_root = '/nfs/';
	
	[edi_file, count, str] = mms_file_search(sc, instr, mode, level, ...
	                                         'OptDesc',   optdesc, ...
	                                         'SDCroot',   sdc_root);
	assert(count > 0, ['EDI file not found: "' str '".']);

%------------------------------------%
% Loop Through Each File             %
%------------------------------------%
	% Dissect the file name to discover which days efield data exists
	[sc, ~, ~, ~, fstart] = mms_dissect_filename(edi_file);
	
	% Convert time to yyyy-mm-ddTHH:MM:SS
	date_temp = MrTimeParser(fstart, '%Y%M%d', '%Y-%M-%d');
	tstart = strcat(date_temp, 'T00:00:00');
	tend   = strcat(date_temp, 'T24:00:00');

	% Loop through each file
	for ii = 1 : count
		try
			cdf_file = mms_sdc_edi_ql_efield(sc{ii}, tstart{ii}, tend{ii});
		catch ME
			fprintf('Unable to create file: %s %s %s\n', sc{ii}, tstart{ii}, tend{ii});
			rethrow(ME);
		end
	end
end