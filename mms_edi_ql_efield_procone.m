%
% Name
%   mms_edi_ql_efield_procone
%
% Purpose
%   Process one day of EDI E-field data, changing L!A to QL data.
%
% Calling Sequence
%   CDF_FILE = mms_edi_ql_efield_procone()
%     Process data from all spacecraft gathered on the date three days
%     before the current date.
%
%   CDF_FILE = mms_edi_ql_efield_procone(SC)
%     Process data from spacecraft SC gathered on the date three days
%     before the current date. If SC is the empty string, all spacecraft
%     are used.
%
%   CDF_FILE = mms_edi_ql_efield_procone(..., DATE_START)
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
%   2015-06-03      Written by Matthew Argall
%   2015-08-23      Renamed from mms_sdc_create_edi_ql_efield to mms_edi_ql_efield_procone. - MRA
%
function cdf_file = mms_edi_ql_efield_procone(sc, date_start)

%------------------------------------%
% Find Files                         %
%------------------------------------%
	% EDI L1A E-Field Data Files
	if nargin < 2
		date_start = datestr( now() - 3.0, 'yyyy-mm-dd' );
	end
	if nargin < 1 || isempty(sc)
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
			cdf_file = mms_edi_ql_efield_create(sc{ii}, tstart, tend);
			
			% Notify where the file has been saved
			if nargout < 1
				fprintf('File written to %s\n', cdf_file);
				clear cdf_file
			end
			
		% Catch error
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