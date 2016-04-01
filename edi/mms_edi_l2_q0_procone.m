%
% Name
%   mms_edi_l2_q0_procone
%
% Purpose
%   Process one day of EDI E-field data to produce L2 Q0 data.
%
% Calling Sequence
%   CDF_FILE = mms_edi_l2_q0_procone()
%     Process data from all spacecraft gathered on the date three days
%     before the current date.
%
%   CDF_FILE = mms_edi_l2_q0_procone(SC)
%     Process data from spacecraft SC gathered on the date three days
%     before the current date. If SC is the empty string, all spacecraft
%     are used.
%
%   CDF_FILE = mms_edi_l2_q0_procone(..., DATE_START)
%     Process data from date DATE_START, formatted as "yyyy-mm-dd".
%
% Parameters
%   SC:             in, optional, type=char/cell, default={'mms1', 'mms2', 'mms3', 'mms4'}
%   DATE_START:     in, optional, type=char/cell, default=three days before current date
%
% Returns
%   CDF_FILE        out, required, type=string
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-09-10      Written by Matthew Argall
%
function cdf_file = mms_edi_l2_q0_procone(sc, mode, date_start)

%------------------------------------%
% Inputs                             %
%------------------------------------%
	% EDI L1A E-Field Data Files
	if nargin < 1 || isempty(sc)
		sc = {'mms1' 'mms2' 'mms3' 'mms4'};
	end
	if nargin < 2 || isempty(mode)
		mode = 'srvy';
	end
	if nargin < 3
		date_start = datestr( now() - 3.0, 'yyyy-mm-dd' );
	end

	% Constants
	create_log_file = true;
	save_dir        = '/nfs/edi/q0/';
	sdc_root        = '/nfs/';

%------------------------------------%
% Create Data                        %
%------------------------------------%
	if ischar(sc)
		sc = { sc };
	end
	
	% Process the entire day.
	tstart   = [date_start 'T00:00:00'];
	tend     = [date_start 'T24:00:00'];
	
	% Allocate memory
	nFiles   = length(sc);
	cdf_file = cell(nFiles);
	
	% Step through each spacecraft.
	for ii = 1 : nFiles
	%------------------------------------%
	% Try to Create the File             %
	%------------------------------------%
		try
			% Create the data
			cdf_file{ii} = mms_edi_l2_q0_create(sc{ii}, tstart, tend, ...
			                                    'CreateLogFile', create_log_file, ...
			                                    'Mode',          mode,            ...
			                                    'SaveDir',       save_dir);
			
			% Notify where the file has been saved
			if nargout < 1
				fprintf('File written to %s\n', cdf_file{ii});
				clear cdf_file
			elseif nFiles == 1
				cdf_file = cdf_file{1};
			end

	%------------------------------------%
	% Process Errors                     %
	%------------------------------------%
		catch ME
			logfile = mrstdlog();
			logfile.alert = true;
			
			% Print error
			mrfprintf('logerr', 'Unable to create file: %s %s %s', sc{ii}, tstart, tend);
			mrfprintf('logerr',  ME);
			
			% Turn alerts back off
			logfile.alert = false;
			
			% No CDF file made
			cdf_file = '';
		end
	end
end