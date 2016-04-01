%
% Name
%   mms_edi_create_ql_efield
%
% Purpose
%   Create a MATLAB save file of inputs needed for Bestarg.
%
% Calling Sequence
%   EDI_QL_FILE = mms_edi_create_ql_efield(SC, TSTART, TEND)
%     Compute the EDI electric field from MMS spacecraft SC during
%     the time interval [TSTART, TEND], and output a CDF file
%     named EDI_QL_FILE.
%
%   EDI_QL_FILE = mms_edi_create_ql_efield(..., 'ParamName', ParamValue)
%     Any parameter name-value pair given below.
%
% Parameters
%   SC:             in, required, type=string
%   TSTART:         in, required, type=string
%   TEND:           in, required, type=string
%   'SDCroot':      in, optional, type=string, default='/nfs/edi/'
%                   Location of the SDC root directory where data is saved.
%   'SaveDir':      in, optional, type=string, default='/nfs/edi/'
%                   Location where output file is saved.
%   'AttitudeDir':  in, optional, type=string, default=SDCroot/ancillary/SC/defatt
%                   Location of definitive attitude data.
%   'HKDir':        in, optional, type=string, default=SDCroot/hk/
%                   Location of housekeeping files.
%   'BeamQuality':  in, optional, type=integer/intarr, default=3
%                   Quality of beams to use. Options are 0, 1, 2, 3, or
%                     any combination of the four.
%   'dt':           in, optional, type=float, default=5.0
%                   Time interval (seconds) over which to average data.
%
% Returns
%   EDI_QL_FILE     out, required, type=string
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-09-10      Written by Matthew Argall
%
function q0_file = mms_edi_l2_q0_create(sc, tstart, tend, varargin)

	% Defaults
	create_log_file = true;
	log_dir         = '/nfs/edi/logs/q0/';
	mode            = 'srvy';
	sdc_root        = '/nfs/';
	save_dir        = '/nfs/edi/q0/';

	% Optional arguments
	nOptArgs = length(varargin);
	for ii = 1 : 2 : nOptArgs
		switch varargin{ii}
			case 'CreateLogFile'
				create_log_file = varargin{ii+1};
			case 'Mode'
				mode = varargin{ii+1};
			case 'SDCroot'
				sdc_root = varargin{ii+1};
			case 'SaveDir'
				save_dir = varargin{ii+1};
			otherwise
				error(['Unknown input parameter: "' varargin{ii} '".']);
		end
	end
	
	% Constants
	instr   = 'edi';
	level   = 'l2';
	optdesc = 'q0';
	version = 'v0.0.2';
	
%------------------------------------%
% Log File                           %
%------------------------------------%
	if create_log_file
		% Reformat the start time
		fstart = MrTimeParser(tstart, '%Y-%M-%dT%H:%m:%S', '%Y%M%d');
	
		% Create a file name
		log_name = mms_construct_filename( sc, instr, mode, level, ...
		                                   'OptDesc',   optdesc,   ...
		                                   'TStart',    fstart,     ...
		                                   'Version',   version );
		% Change extension to 'log'
		[~, log_name] = fileparts(log_name);
		log_name      = fullfile( log_dir, [log_name '.log']);
		
		% Create log object
		%   - Delete the log file if no errors occurred.
		logfile = mrstdlog(log_name);
		logfile.delete_on_destroy = true;
	else
		mrstdlog('');
	end

%------------------------------------%
% Find Files                         %
%------------------------------------%

	% EDI Slow L1A E-Field Data File
	%    - Find last, so file descriptors are saved.
	[edi_slow_file, slw_cnt, str] = mms_file_search(sc, instr, 'slow', 'l1a', ...
	                                                'TStart',    tstart, ...
	                                                'TEnd',      tend, ...
	                                                'OptDesc',   'efield', ...
	                                                'SDCroot',   sdc_root);
	
	% EDI Fast L1A E-Field Data File
	%    - Find last, so file descriptors are saved.
	[edi_fast_file, fst_cnt, str] = mms_file_search(sc, instr, 'fast', 'l1a', ...
	                                                'TStart',    tstart, ...
	                                                'TEnd',      tend, ...
	                                                'OptDesc',   'efield', ...
	                                                'SDCroot',   sdc_root);

	% Write names to log file
	if create_log_file
		mrfprintf('logtext', 'Parent files:');
		mrfprintf('logtext', strcat('  EDI slow: "', edi_slow_file, '"') );
		mrfprintf('logtext', strcat('  EDI fast: "', edi_fast_file, '"') );
	end

%------------------------------------%
% Read Data                          %
%------------------------------------%

	% EDI Fast data
	edi_slow = [];
	if slw_cnt > 0
		edi_slow = mms_edi_read_l1a_efield( edi_slow_file, tstart, tend );
	end

	% EDI Slow data
	edi_fast = [];
	if fst_cnt > 0
		%
		% Occasionally, a file will be found with zero records. This is a
		% problem with L1A processing which produces files when no data
		% is available.
		%
		try
			edi_fast = mms_edi_read_l1a_efield( edi_fast_file, tstart, tend );
		% Note the error and continue processing slow survey data.
		catch ME
			mrfprintf('logerr', ME);
			edi_fast_file = '';
			fst_cnt       = 0;
		end
	end
	
	% Check after we check for empty fast survey file.
	assert( slw_cnt + fst_cnt > 0, 'Unable to find fast or slow survey EDI files.' );

	% Combine slow and fast survey data
	%   - Early versions of L1A data had energy labelled incorrectly.
	try
		edi = mms_edi_srvy_combine( edi_slow, edi_fast );
		if fst_cnt == 0 && length(edi_slow.energy_gd12) ~= length(edi_slow.tt2000_gd12)
			error('Bad time tags for energy.');
		end
		if slw_cnt == 0 && length(edi_fast.energy_gd12) ~= length(edi_fast.tt2000_gd21)
			error('Bad time tags for energy.');
		end
	catch ME
		mrfprintf('logerr', ME);

		% Make energy empty
		if slw_cnt > 0
			edi_slow = rmfield( edi_slow, { 'energy_gd12', 'energy_gd21' } );
		end
		if fst_cnt > 0
			edi_fast = rmfield( edi_fast, { 'energy_gd12', 'energy_gd21' } );
		end
		edi                 = mms_edi_srvy_combine( edi_slow, edi_fast );
		edi.('energy_gd12') = [];
		edi.('energy_gd21') = [];
	end
		
	clear edi_slow edi_fast
	
%------------------------------------%
% Create Q0 Data                     %
%------------------------------------%
	% Find quality zero data
	iq0_gd12 = find(edi.quality_gd12 == 0);
	iq0_gd21 = find(edi.quality_gd21 == 0);

	% Extract times
	t_gd12 = edi.tt2000_gd12(iq0_gd12);
	t_gd21 = edi.tt2000_gd21(iq0_gd21);

	% Extract Energies
	if isempty(edi.energy_gd12)
		e_gd12 = uint16([]);
		e_gd21 = uint16([]);
	else
		e_gd12 = edi.energy_gd12(iq0_gd12);
		e_gd21 = edi.energy_gd21(iq0_gd21);
	end

	% Extract counts
	counts_gd12 = edi.word15_gd12(iq0_gd12);
	counts_gd21 = edi.word15_gd21(iq0_gd21);
	
	% Clear data we do not want anymore
	clear edi iq0_gd12 iq0_gd21
	
%------------------------------------%
% Gather Data                        %
%------------------------------------%
	% Parent files
	%   - Concatenate slow and fast survey file names
	%   - Remove filenames if any were empty
	%   - Strip off the directory.
	parents         = { edi_fast_file edi_slow_file };
	iempty          = find( cellfun(@isempty, parents) );
	parents(iempty) = [];
	[~, name, ext]  = cellfun(@fileparts, parents, 'UniformOutput', false);
	parents         = strcat(name, ext);
	
	% Metadata structure
	meta = struct( 'sc',           sc,         ...
	               'instr',        instr,      ...
	               'mode',         mode,       ...
	               'optdesc',      optdesc,    ...
	               'directory',    save_dir,   ...
	               'tstart',       tstart,     ...
	               'tend',         tend,       ...
	               'parents',      { parents } ...     % Prevent array of structures
	             );
	
	
	% E-field Beam Convergence
	q0_data = struct( 'tt2000_gd12',  t_gd12,      ...
	                  'tt2000_gd21',  t_gd21,      ...
	                  'energy_gd12',  e_gd12,      ...
	                  'energy_gd21',  e_gd21,      ...
	                  'counts_gd12',  counts_gd12, ...
	                  'counts_gd21',  counts_gd21  ...
	                );
	clear t_gd12 t_gd21 e_gd12 e_gd21 counts_gd12 counts_gd21

%------------------------------------%
% Write to File                      %
%------------------------------------%
	q0_file = mms_edi_l2_q0_write(meta, q0_data);
end
