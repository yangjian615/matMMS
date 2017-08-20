%
% Name
%   mms_edi_ql_efield_procall
%
% Purpose
%   Process all EDI E-field, changing L1A data into QL data.
%
% Calling Sequence
%   FILES = mms_sdc_proc_edi_ql_efield()
%     Process all data from all spacecraft from the beginning of the mission
%     until three days before the current date.
%
%   FILES = mms_sdc_proc_edi_ql_efield(SC, TSTART, TEND)
%     Read and process all data from MMS spacecraft SC between the times
%     TSTART and TEND. Data is saved to CDF files named FILES.
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
%   2015-08-12      Look for fast & slow files and take union of dates. - MRA
%   2015-08-22      Renamed from mms_sdc_proc_edi_ql_efield to mms_edi_ql_efield_procall. - MRA
%
function files = mms_edi_ql_efield_procall(sc, mode, tstart, tend)

%------------------------------------%
% Find Files                         %
%------------------------------------%
	% EDI L1A E-Field Data Files
	if nargin < 4
		tend = datestr( now() - 3.0, 'yyyy-mm-dd' );
	end
	if nargin < 3 || isempty(tstart)
		tstart = '2015-03-17';
	end
	if nargin < 2 || isempty(tstart)
		mode = 'srvy';
	end
	if nargin < 1 || isempty(sc)
		sc = { 'mms1' 'mms2' 'mms3' 'mms4' };
	end
	
	% Constants
	beam_quality    = 3;
	create_log_file = true;
	dt              = 5.0;
	log_dir         = '/nfs/edi/logs/ql/';
	mode            = 'srvy';
	sdc_root        = '/nfs/';
	ql_dir          = '/nfs/edi/ql/';
	sl_dir          = '/nfs/edi/sl/';
	attitude_dir    = fullfile(sdc_root, 'ancillary', sc, 'defatt');
	hk_root         = fullfile(sdc_root, 'hk');

%------------------------------------%
% Process Data for Each Spacecraft   %
%------------------------------------%
	% Make sure we have a cell.
	if ischar(sc)
		sc = { sc };
	end
	
	% Full day
	tstart = [tstart 'T00:00:00'];
	tend   = [tend   'T24:00:00'];
	
	files = {};
	count = 0;

	% Step through each spacecraft
	for ii = 1 : length(sc)
	%------------------------------------%
	% Find Dates with Data               %
	%------------------------------------%
		% Find all slow files
		[slow_files, nSlow] = mms_file_search( sc{ii}, 'edi', 'slow', 'l1a', ...
		                                       'OptDesc', 'efield', ...
		                                       'TStart',  tstart,   ...
		                                       'TEnd',    tend );
		
		% Find all fast files
		[fast_files, nFast] = mms_file_search( sc{ii}, 'edi', 'fast', 'l1a', ...
		                                       'OptDesc', 'efield', ...
		                                       'TStart',  tstart,   ...
		                                       'TEnd',    tend );
		
		% Extract the start times
		if nSlow > 0 && nFast > 0
			files = [slow_files fast_files];
		elseif nSlow > 0
			files = slow_files;
		elseif nFast > 0
			files = fast_files;
		end
		
		% Dissect the file names
		[~, ~, ~, ~, fstart] = mms_dissect_filename( files );
	
		% All dates in which we have data
		fdates = unique( fstart );

		% Reformat them
		dates  = MrTimeParser(fdates, '%Y%M%d', '%Y-%M-%d');
		fstart = strcat( dates, 'T00:00:00' );
		fend   = strcat( dates, 'T24:00:00' );

	%------------------------------------%
	% Loop Through Each Date             %
	%------------------------------------%
		for jj = 1 : length(dates)
			try
				% Create the data
				files{count+1} = mms_edi_ql_efield_create(sc{ii}, fstart{jj}, fend{jj},     ...
				                                          'AttitudeDir',   attitude_dir{ii},...
				                                          'BeamQuality',   beam_quality,    ...
				                                          'CreateLogFile', create_log_file, ...
				                                          'dt',            dt,              ...
				                                          'HKdir',         hk_root,         ...
				                                          'Mode',          mode,            ...
				                                          'SDCroot',       sdc_root,        ...
				                                          'SlowLookDir',   sl_dir,          ...
				                                          'QuickLookDir',  ql_dir           ...
				                                        );
				
				% Count the file
				count = count + 1;
				
				if nargout < 1
					fprintf( ['File written to: "' files{count} '".\n'] );
				end
			catch ME
				% Alert me
				logfile = mrstdlog();
				logfile.alert = true;
			
				% Print error
				mrfprintf('logerr', 'Unable to create file: %s %s %s', sc{ii}, fstart{jj}, fend{jj});
				mrfprintf('logerr',  ME);
				
				% Turn alerts back off
				logfile.alert = false;
			end
		end
	end
	
	% Trim files
	files = files(1:count);
end
