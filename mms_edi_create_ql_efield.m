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
%   2015-06-03      Written by Matthew Argall
%   2015-06-20      Writing to CDF file now occurs in separate function. - MRA
%   2015-07-20      Renamed from mms_sdc_edi_ql_efield to mms_edi_create_ql_efield. - MRA
%
function edi_ql_file = mms_edi_create_ql_efield(sc, tstart, tend, varargin)

% MMS2: May 9, 2015  16:08 - 16:13
% MMS4: May 6, 2015  15:30 - 15:35

	% Defaults
	sdc_root        = '/nfs/';
	save_dir        = '/nfs/edi/temp/';
	log_dir         = '/nfs/edi/logs/';
	attitude_dir    = fullfile(sdc_root, 'ancillary', sc, 'defatt');
	hk_root         = fullfile(sdc_root, 'hk');
	beam_quality    = 3;
	create_log_file = true;
	dt              = 5.0;

	% Optional arguments
	nOptArgs = length(varargin);
	for ii = 1 : 2 : nOptArgs
		switch varargin{ii}
			case 'AttitudeDir'
				attitude_dir = varargin{ii+1};
			case 'BeamQuality'
				beam_quality = varargin{ii+1};
			case 'CreateLogFile'
				create_log_file = varargin{ii+1};
			case 'dt'
				dt = varargin{ii+1};
			case 'HKdir'
				hk_root = varargin{ii+1};
			case 'SDCroot'
				sdc_root = varargin{ii+1};
			case 'SaveDir'
				save_dir = varargin{ii+1};
			otherwise
				error(['Unknown input parameter: "' varargin{ii} '".']);
		end
	end

%------------------------------------%
% Log File                           %
%------------------------------------%
	if create_log_file
		% Reformat the start time
		fstart = MrTimeParser(tstart, '%Y-%M-%dT%H:%m:%S', '%Y%M%d');
	
		% Create a file name
		log_name = mms_construct_filename( sc, 'edi', 'srvy', 'ql', ...
		                                   'OptDesc',   'efield',   ...
		                                   'TStart',    fstart,     ...
		                                   'Version',   'v0.2.0' );
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

	% FG QL Data File
	instr   = 'dfg';
	mode    = 'srvy';
	level   = 'ql';
	optdesc = '';
	[fg_file, count, str] = mms_file_search(sc, instr, mode, level, ...
	                                        'TStart',    tstart, ...
	                                        'TEnd',      tend, ...
	                                        'OptDesc',   optdesc, ...
	                                        'SDCroot',   sdc_root);
	assert(count > 0, ['DFG file not found: "' str '".']);

	% FG L1A Data File
	instr   = 'dfg';
	mode    = 'srvy';
	level   = 'l1b';
	optdesc = '';
	[fg_l1b_file, count, str] = mms_file_search(sc, instr, mode, level, ...
	                                            'TStart',    tstart, ...
	                                            'TEnd',      tend, ...
	                                            'OptDesc',   optdesc, ...
	                                            'SDCroot',   sdc_root);
	assert(count > 0, ['DFG L1B file not found: "' str '".']);

	% Digital Sun Sensor L1B Data File
	instr   = 'fields';
	mode    = 'hk';
	level   = 'l1b';
	optdesc = '101';
	[dss_file, count, str] = mms_file_search(sc, instr, mode, level, ...
	                                         'TStart',    tstart, ...
	                                         'TEnd',      tend, ...
	                                         'OptDesc',   optdesc, ...
	                                         'SDCroot',   hk_root);
	assert(count > 0, ['DFG file not found: "' str '".']);

	% Attitude file
	%   - Do not throw errors for attitude files. They are used only
	%     to rotate from BCS to SMPA, which is very nearly a unitary
	%     transformation.
	%   - Sunpulse times are used to despin data.
	str = fullfile( attitude_dir, [upper(sc) '_DEFATT_%Y%D_%Y%D.V*'] );
	[att_file, att_cnt] = MrFile_Search( str, ...
	                                    'Closest',      true, ...
	                                    'TStart',       tstart, ...
	                                    'TEnd',         tend, ...
	                                    'TimeOrder',    '%Y%D', ...
	                                    'VersionRegex', 'V([0-9]{2})' );
	
	% EDI Slow L1A E-Field Data File
	%    - Find last, so file descriptors are saved.
	instr   = 'edi';
	mode    = 'slow';
	level   = 'l1a';
	optdesc = 'efield';
	[edi_slow_file, slw_cnt, str] = mms_file_search(sc, instr, mode, level, ...
	                                                'TStart',    tstart, ...
	                                                'TEnd',      tend, ...
	                                                'OptDesc',   optdesc, ...
	                                                'SDCroot',   sdc_root);
	
	% EDI Fast L1A E-Field Data File
	%    - Find last, so file descriptors are saved.
	instr   = 'edi';
	mode    = 'fast';
	level   = 'l1a';
	optdesc = 'efield';
	[edi_fast_file, fst_cnt, str] = mms_file_search(sc, instr, mode, level, ...
	                                                'TStart',    tstart, ...
	                                                'TEnd',      tend, ...
	                                                'OptDesc',   optdesc, ...
	                                                'SDCroot',   sdc_root);

	% Write names to log file
	if create_log_file
		mrfprintf('logtext', 'Parent files:');
		mrfprintf('logtext', strcat('  FG:       "', fg_file,  '"') );
		mrfprintf('logtext', strcat('  DSS:      "', dss_file, '"') );
		if att_cnt > 1
			mrfprintf('logtext', strcat('  Defatt:   "', strjoin(att_file, '\n            ') ) );
		else
			mrfprintf('logtext', strcat('  Defatt:   "', att_file, '"') );
		end
		mrfprintf('logtext', strcat('  EDI slow: "', edi_slow_file, '"') );
		mrfprintf('logtext', strcat('  EDI fast: "', edi_fast_file, '"') );
	end

%------------------------------------%
% Read Data                          %
%------------------------------------%
	% Attitude data
	if isempty(att_file)
		defatt  = [];
		zMPA    = [];
	else
		[defatt, att_hdr] = mms_fdoa_read_defatt(att_file, tstart, tend);
		zMPA = att_hdr.zMPA(:,1);
	end
	
	% Sunpulse data
	sunpulse = mms_dss_read_sunpulse(dss_file, tstart, tend, 'UniquePulse', true);

	% EDI Fast data
	if slw_cnt > 0
		edi_slow = mms_edi_create_l2pre_efield( edi_slow_file, tstart, tend, ...
		                                        'Attitude', defatt,          ...
		                                        'CS_GSE',   false,           ...
		                                        'CS_DMPA',  true,            ...
		                                        'CS_123',   true,            ...
		                                        'Quality',  beam_quality,    ...
		                                        'Sunpulse', sunpulse,        ...
		                                        'zMPA',     zMPA );
	end

	% EDI Slow data
	if fst_cnt > 0
		%
		% Occasionally, a file will be found with zero records. This is a
		% problem with L1A processing which produces files when no data
		% is available.
		%
		try
			edi_fast = mms_edi_create_l2pre_efield( edi_fast_file, tstart, tend, ...
			                                        'Attitude', defatt,          ...
			                                        'CS_GSE',   false,           ...
			                                        'CS_DMPA',  true,            ...
			                                        'CS_123',   true,            ...
			                                        'Quality',  beam_quality,    ...
			                                        'Sunpulse', sunpulse,        ...
			                                        'zMPA',     zMPA );
		% Not the error and continue processing slow survey data.
		catch ME
			mrfprintf('logerr', ME);
			fst_cnt = 0;
		end
	end
	
	% Check after we check for empty fast survey file.
	assert( slw_cnt + fst_cnt > 0, 'Unable to find fast or slow survey EDI files.' );

	% FGM data
	fg_ql  = mms_fg_read_ql(fg_file,      tstart, tend);
	fg_l1b = mms_fg_read_l1b(fg_l1b_file, tstart, tend);

%------------------------------------%
% Combine Slow and Fast Survey       %
%------------------------------------%
	% Combine slow and fast data
	%   - Already checked the fst_cnt == 0 && slw_cnt == 0 case
	if fst_cnt > 0 && slw_cnt > 0
		% Create the combined structure
		edi = struct();
		
		% Sort and combine times.
		[edi.('tt2000_gd12'), isort_gd12] = sort( [edi_slow.tt2000_gd12 edi_fast.tt2000_gd12] );
		[edi.('tt2000_gd21'), isort_gd21] = sort( [edi_slow.tt2000_gd21 edi_fast.tt2000_gd21] );
		
		% Find unique times
		[~, iuniq_gd12] = unique( edi.tt2000_gd12 );
		[~, iuniq_gd21] = unique( edi.tt2000_gd21 );
		
		% Sorted, unique order
		iorder_gd12 = isort_gd12( iuniq_gd12 );
		iorder_gd21 = isort_gd21( iuniq_gd21 );
		
		% Remove the time fields from the slow and fast structures
		rmfield(edi_slow, { 'tt2000_gd12', 'tt2000_gd21'} );
		rmfield(edi_fast, { 'tt2000_gd12', 'tt2000_gd21'} );
	
		% Get the names of each remaining field.
		fields  = fieldnames( edi_slow );
		nFields = length( fields );
		
		% Find fields related to gd12
		igd12 = find( ~cellfun( @isempty, regexp( fields, '(gd12|gun1)', 'once' ) ) );
		igd21 = find( ~cellfun( @isempty, regexp( fields, '(gd21|gun2)', 'once' ) ) );

		% Step through each field.
		%   - Both guns have the same fields.
		for ii = 1 : nFields / 2
			%  Combine data; select sorted, unique elements; add to EDI structure
			field_gd12       = fields{ igd12(ii) };
			tmp_data         = [ edi_slow.( field_gd12 ) edi_fast.( field_gd12 ) ];
			tmp_data         = tmp_data(:, iorder_gd12);
			edi.(field_gd12) = tmp_data;
			
			% Repeat for GD21
			field_gd21       = fields{ igd21(ii) };
			tmp_data         = [ edi_slow.( field_gd21 ) edi_fast.( field_gd21 ) ];
			tmp_data         = tmp_data(:, iorder_gd21);
			edi.(field_gd21) = tmp_data;
		end

		% Clear data that we no longer need
		clear edi_fast edi_slow tmp_data iuniq_gd12 iuniq_gd21 isort_gd12 isort_gd21 ...
		      iorder_gd12 iorder_gd21 fields nFields igd12 igd21 field_gd12 field_gd21
	
	% Take only fast survey data.
	elseif fst_cnt > 0
		edi = edi_fast;
		clear edi_fast
		
	% Take only slow survey data.
	else
		edi = edi_slow;
		clear edi_slow
	end

%------------------------------------%
% Compute E-Field                    %
%------------------------------------%

	% Beam Convergence Method
	[efield_bc, b_interp] = ...
		mms_edi_calc_efield_bc( fg_ql.tt2000, fg_ql.b_dmpa(1:3,:), ...
		                        edi.tt2000_gd12, edi.virtual_gun1_dmpa, edi.fv_gd12_dmpa, ...
		                        edi.tt2000_gd21, edi.virtual_gun2_dmpa, edi.fv_gd21_dmpa, dt );
	
	% Cost Function method
	efield_cf = mms_edi_calc_efield_cf( fg_l1b.tt2000, fg_l1b.b_bcs(1:3,:), ...
	                                    edi.tt2000_gd12, edi.virtual_gun1_123, edi.fv_gd12_123, ...
	                                    edi.tt2000_gd21, edi.virtual_gun2_123, edi.fv_gd21_123, dt );

	% Rotate from BCS to DMPA
	if ~isempty(zMPA)
		% Create the transformation and transform
		bcs2smpa = mms_fg_xbcs2smpa(zMPA);
		E_smpa   = mrvector_rotate( bcs2smpa, efield_cf.E_cf_bcs );
		v_smpa   = mrvector_rotate( bcs2smpa, efield_cf.v_cf_bcs );
		d_smpa   = mrvector_rotate( bcs2smpa, efield_cf.d_cf_bcs );
	else
		mrfprintf('logwarn', 'z-MPA axis not available. Cannot rotate cost function results to SMPA.');
		E_smpa = efield_cf.E_cf_bcs;
		v_smpa = efield_cf.v_cf_bcs;
		d_smpa = efield_cf.d_cf_bcs;
	end
	
	% Remove BCS fields
	efield_cf = rmfield( efield_cf, {'E_cf_bcs', 'v_cf_bcs', 'd_cf_bcs'} );
	
	% Despin and create structure tags for data
	smpa2dmpa = mms_dss_xdespin( sunpulse, efield_cf.tt2000_cf );
	efield_cf.('E_cf_dmpa') = single( mrvector_rotate( smpa2dmpa, E_smpa ) );
	efield_cf.('v_cf_dmpa') = single( mrvector_rotate( smpa2dmpa, v_smpa ) );
	efield_cf.('d_cf_dmpa') = single( mrvector_rotate( smpa2dmpa, d_smpa ) );

	% Clear unneeded data.
	clear fg_ql fg_l1b defatt att_hdr zMPA sunpulse E_smpa v_smpa d_smpa

%------------------------------------%
% Gather Data                        %
%------------------------------------%
	% Parent files
	if iscell(att_file)
		parents = [ fg_file dss_file att_file edi_fast_file edi_slow_file ];
	else
		parents = { fg_file dss_file att_file edi_fast_file edi_slow_file };
	end
	iempty  = find( cellfun(@isempty, parents) );
	parents(iempty) = [];
	
	% Metadata structure
	meta = struct( 'sc',           sc,         ...
	               'instr',        instr,      ...
	               'mode',         'srvy',     ...
	               'optdesc',      optdesc,    ...
	               'directory',    save_dir,   ...
	               'tstart',       tstart,     ...
	               'tend',         tend,       ...
	               'parents',      { parents } ...             % Prevent array of structures
	             );

	% Gather data to be written
	%   - Record varying dimension must be first
	beams = struct( 'tt2000_gd12',  edi.tt2000_gd12,               ...
	                'tt2000_gd21',  edi.tt2000_gd21,               ...
	                'pos_vg1_dmpa', single(edi.virtual_gun1_dmpa), ...
	                'pos_vg2_dmpa', single(edi.virtual_gun2_dmpa), ...
	                'fv_gd12_dmpa', single(edi.fv_gd12_dmpa),      ...
	                'fv_gd21_dmpa', single(edi.fv_gd21_dmpa),      ...
	                'quality_gd12', edi.quality_gd12,              ...
	                'quality_gd21', edi.quality_gd21               ...
	              );

	% Clear data
	clear edi b_avg

	% Write data to a file
	edi_ql_file = mms_edi_write_ql_efield(meta, beams, efield_cf, efield_bc, b_interp);
end