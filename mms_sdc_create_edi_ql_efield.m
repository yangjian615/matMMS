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
function [] = mms_sdc_create_edi_ql_efield()

%------------------------------------%
% Find Files                         %
%------------------------------------%
	% EDI L1A E-Field Data Files
	sc       = { 'mms4' }; % {'mms1' 'mms2' 'mms3' 'mms4'};
	instr    = 'edi';
	mode     = 'slow';
	level    = 'l1a';
	optdesc  = 'efield';
	theDate  = '2015-05-02';
	
% Error
%	theDate  = '2015-05-02';   % mms4
	
% SpinPeriod = 0
%	theDate  = '2015-04-22';   % mms3
%	theDate  = '2015-04-24';   % mms2

%------------------------------------%
% Create Data                        %
%------------------------------------%
	
	% The date three days ago.
	%   - Need to wait for the edi_l1a and dfg_ql files to be created.
	if isempty(theDate)
		theDate = datestr( now() - 3.0, 'yyyy-mm-dd' );
	end
	tstart  = [theDate 'T00:00:00'];
	tend    = [theDate 'T24:00:00'];
	
	% Step through each spacecraft.
	for ii = 1 : length(sc)
		% Search for the 
		[edi_file, count, str] = mms_file_search(sc{ii}, instr, mode, level, ...
		                                         'OptDesc', optdesc, ...
		                                         'TStart',  tstart, ...
		                                         'TEnd',    tend);

		% Skip files that were not found.
		if count == 0
			fprintf('EDI file not found: "%s".\n', str);
			continue
		end
		
		% Create the E-field file.
		try
			cdf_file = mms_edi_create_ql_efield(sc{ii}, tstart, tend);
		catch ME
			fprintf('Unable to create file: %s %s\n', sc{ii}, tstart);
			rethrow(ME);
		end
	end
end