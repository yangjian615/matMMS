%
% Name
%   mms_sdc_create_edi_ql_efield
%
% Purpose
%   Find all EDI L1A electric field-mode data files and calculate the electric field.
%
% Calling Sequence
%   SAVE_FILE = mms_sdc_create_edi_ql_efield()
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
function cdf_file = mms_sdc_create_edi_ql_efield(sc, date_start)

%------------------------------------%
% Find Files                         %
%------------------------------------%
	% EDI L1A E-Field Data Files
	if nargin < 2
		date_start = datestr( now() - 3.0, 'yyyy-mm-dd' );
	end
	if nargin < 1
		sc = {'mms1' 'mms2' 'mms3' 'mms4'};
	end
	
% Error
%	theDate  = '2015-05-02';   % mms4
	
% SpinPeriod = 0
%	theDate  = '2015-04-22';   % mms3
%	theDate  = '2015-04-24';   % mms2

%------------------------------------%
% Create Data                        %
%------------------------------------%
	if ischar(sc)
		sc = { sc };
	end
	
	% Process the entire day.
	tstart  = [date_start 'T00:00:00'];
	tend    = [date_start 'T24:00:00'];
	
	% Step through each spacecraft.
	for ii = 1 : length(sc)
		% Create the E-field file.
		try
			% Create the data
			cdf_file = mms_edi_create_ql_efield(sc{ii}, tstart, tend);
			
			% Notify where the file has been saved
			if nargout < 1
				fprintf('File written to %s\n', cdf_file);
			end
			
		% Catch error
		catch ME
			% Print error
			fprintf('Unable to create file: %s %s %s\n', sc{ii}, tstart, tend);
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