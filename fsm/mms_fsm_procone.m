%
% Name
%   mms_edi_ql_efield_procone
%
% Purpose
%   Process one day of EDI E-field data, changing L!A to QL data.
%
% Calling Sequence
%   CDF_FILE = mms_fsm_ql_procone()
%     Process data from all spacecraft gathered on the date three days
%     before the current date.
%
%   CDF_FILE = mms_fsm_ql_procone(SC)
%     Process data from spacecraft SC gathered on the date three days
%     before the current date. If SC is the empty string, all spacecraft
%     are used.
%
%   CDF_FILE = mms_fsm_ql_procone(..., MODE)
%     Process data for telemetry mode MODE. If MODE is the empty string,
%     survey mode data is processed (as opposed to burst mode data).
%
%   CDF_FILE = mms_fsm_ql_procone(..., DATE_START)
%     Process data from date DATE_START, formatted as "yyyy-mm-dd".
%
% Parameters
%   SC:             in, optional, type=char/cell, default={'mms1', 'mms2', 'mms3', 'mms4'}
%   MODE:           in, optional, type=char/cell, default='srvy'
%   DATE_START:     in, optional, type=char/cell, default=three days before current date
%
% Returns
%   CDF_FILE        out, required, type=string
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-09-22      Written by Matthew Argall
%
function cdf_file = mms_fsm_ql_procone(sc, mode, date_start)

%------------------------------------%
% Inputs                             %
%------------------------------------%
	% EDI L1A E-Field Data Files
	if nargin < 3 || isempty(date_start)
		date_start = datestr( now() - 3.0, 'yyyy-mm-dd' );
	end
	if nargin < 2 || isempty(mode)
		mode = 'srvy';
	end
	if nargin < 1 || isempty(sc)
		sc = {'mms1' 'mms2' 'mms3' 'mms4'};
	end
	
	% Constants
	duration    = 300.0;
	f_max       = 5.0;
	f_min       = 1.0;
	fgm_instr   = 'dfg';
	fgm_mode    = mode;
	fgm_optdesc = '';
	save_dir    = '/nfs/fsm/temp';
	scm_mode    = mode;
	scm_optdesc = 'scm';
	scm_cal_dir = '/home/argall/data/mms/scm_cal/';

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
			slow_data = mms_fsm_create_ql( sc{ii}, tstart, tend,       ...
			                               'Duration',    duration,    ...
			                               'f_max',       f_max,       ...
			                               'f_min',       f_min,       ...
			                               'fgm_instr',   fgm_instr,   ...
			                               'fgm_mode',    'slow',      ...
			                               'fgm_optdesc', fgm_optdesc, ...
			                               'SaveDir',     save_dir,    ...
			                               'scm_mode',    'slow',      ...
			                               'scm_optdesc', 'scs',       ...
			                               'scm_cal_dir', scm_cal_dir  ...
			                             );

			% Fast survey data
			fast_data = mms_fsm_create_ql( sc{ii}, tstart, tend,       ...
			                               'Duration',    duration,    ...
			                               'f_max',       f_max,       ...
			                               'f_min',       f_min,       ...
			                               'fgm_instr',   fgm_instr,   ...
			                               'fgm_mode',    'fast',      ...
			                               'fgm_optdesc', fgm_optdesc, ...
			                               'SaveDir',     save_dir,    ...
			                               'scm_mode',    'fast',      ...
			                               'scm_optdesc', 'scf',       ...
			                               'scm_cal_dir', scm_cal_dir  ...
			                             );

			% Combine the data
			tt2000     = [ slow_data.tt2000;     fast_data.tt2000     ];
			tt2000_fgm = [ slow_data.tt2000_fgm; fast_data.tt2000_fgm ];
			b_fsm_dmpa = [ slow_data.b_fsm_dmpa; fast_data.b_fsm_dmpa ];
			b_fgm_dmpa = [ slow_data.b_fgm_dmpa; fast_data.b_fgm_dmpa ];
			
			% Sort data
			[tt2000,     ifsm] = sort(tt2000);
			[tt2000_fgm, ifgm] = sort(tt2000_fgm);
			b_fsm_dmpa         = b_fsm_dmpa(ifsm,:);
			b_fgm_dmpa         = b_fgm_dmpa(ifgm,:);
			clear ifsm ifgm

			% Parents
			parents = unique( [slow_data.parents fast_data.parents] );
			
			% Destroy old structures and create new
			clear slow_data fast_data
			fsm_data = struct( 'tt2000',      tt2000, ...
			                   'tt2000_fgm',  tt2000_fgm, ...
			                   'b_fsm_dmpa',  b_fsm_dmpa, ...
			                   'b_fgm_dmpa',  b_fgm_dmpa, ...
			                   'sc',          sc,         ...
			                   'mode',        'srvy',     ...
			                   'instr',       [fgm_instr '-scm'], ...
			                   'tstart',      tstart, ...
			                   'directory',   save_dir, ...
			                   'parents',     { parents } ...
			                 );

			% Write to file
			cdf_file = mms_fsm_write_ql( fsm_data );


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