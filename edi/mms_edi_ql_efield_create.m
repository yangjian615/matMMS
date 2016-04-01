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
function edi_ql_file = mms_edi_ql_efield_create(sc, tstart, tend, varargin)

	% Defaults
	sdc_root        = '/nfs/';
	attitude_dir    = fullfile(sdc_root, 'ancillary', sc, 'defatt');
	beam_quality    = 3;
	create_log_file = true;
	dt              = 5.0;
	hk_root         = fullfile(sdc_root, 'hk');
	log_dir         = '/nfs/edi/logs/';
	mode            = 'srvy';
	ql_dir          = '/nfs/edi/ql/';
	sl_dir          = '/nfs/edi/sl/';

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
			case 'Mode'
				mode = varargin{ii+1};
			case 'SDCroot'
				sdc_root = varargin{ii+1};
			case 'SlowLookDir'
				ql_dir = varargin{ii+1};
			case 'QuickLookDir'
				ql_dir = varargin{ii+1};
			otherwise
				error(['Unknown input parameter: "' varargin{ii} '".']);
		end
	end
	
	% Constants
	instr   = 'edi';
	level   = 'ql';
	optdesc = 'efield';
	version = 'v0.4.0';
	
%------------------------------------%
% Log File                           %
%------------------------------------%
	if create_log_file
		% Reformat the start time
		fstart = MrTimeParser(tstart, '%Y-%M-%dT%H:%m:%S', '%Y%M%d');
	
		% Create a file name
		log_name = mms_construct_filename( sc, instr, mode, level,  ...
		                                   'OptDesc',   optdesc,    ...
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

	% FG L1B Data File
	[fg_l1b_file, count, str] = mms_file_search(sc, 'dfg', 'srvy', 'l1b', ...
	                                            'TStart',    tstart,      ...
	                                            'TEnd',      tend,        ...
	                                            'SDCroot',   sdc_root);
	assert(count > 0, ['DFG L1B file not found: "' str '".']);

	% FG L1B Data File
	[fg_ql_file, count, str] = mms_file_search(sc, 'dfg', 'srvy', 'ql', ...
	                                           'TStart',    tstart,     ...
	                                           'TEnd',      tend,       ...
	                                           'SDCroot',   sdc_root);
	assert(count > 0, ['DFG L1B file not found: "' str '".']);

	% Digital Sun Sensor L1B Data File
	[dss_file, count, str] = mms_file_search(sc, 'fields', 'hk', 'l1b', ...
	                                         'TStart',    tstart,       ...
	                                         'TEnd',      tend,         ...
	                                         'OptDesc',   '101',        ...
	                                         'SDCroot',   hk_root);
	assert(count > 0, ['DFG file not found: "' str '".']);

	% Attitude file
	%   - Do not throw errors for attitude files. They are used only
	%     to rotate from BCS to SMPA, which is very nearly a unitary
	%     transformation.
	%   - Sunpulse times are used to despin data.
	str = fullfile( attitude_dir, [upper(sc) '_DEFATT_%Y%D_%Y%D.V*'] );
	[att_file, att_cnt] = MrFile_Search( str, ...
	                                     'Closest',      true,   ...
	                                     'TStart',       tstart, ...
	                                     'TEnd',         tend,   ...
	                                     'TimeOrder',    '%Y%D', ...
	                                     'VersionRegex', 'V([0-9]{2})' );
	
	% EDI Slow L1A E-Field Data File
	%    - Find last, so file descriptors are saved.
	[edi_slow_file, slw_cnt, str] = mms_file_search(sc, instr, 'slow', 'l1a', ...
	                                                'TStart',    tstart,      ...
	                                                'TEnd',      tend,        ...
	                                                'OptDesc',   optdesc,     ...
	                                                'SDCroot',   sdc_root);
	
	% EDI Fast L1A E-Field Data File
	%    - Find last, so file descriptors are saved.
	[edi_fast_file, fst_cnt, str] = mms_file_search(sc, instr, 'fast', 'l1a', ...
	                                                'TStart',    tstart,      ...
	                                                'TEnd',      tend,        ...
	                                                'OptDesc',   optdesc,     ...
	                                                'SDCroot',   sdc_root);

	% Write names to log file
	if create_log_file
		mrfprintf('logtext', 'Parent files:');
		mrfprintf('logtext', strcat('  FG L1B:   "', fg_l1b_file,  '"') );
		mrfprintf('logtext', strcat('  FG QL:    "', fg_ql_file,  '"') );
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
		mrfprintf('logwarn', 'No attitude data. Assuming zMPA = [0; 0; 1].');
		defatt  = [];
		zMPA    = [0; 0; 1];
	else
		[defatt, att_hdr] = mms_fdoa_read_defatt(att_file, tstart, tend);
		zMPA = att_hdr.zMPA(:,1);
	end
	
	% Sunpulse data
	sunpulse = mms_dss_read_sunpulse(dss_file, tstart, tend, 'UniquePulse', true);

	% EDI Fast data
	edi_slow = [];
	if slw_cnt > 0
		edi_slow = mms_edi_l1b_efield_create( edi_slow_file, tstart, tend, ...
		                                        'CS_BCS',   true,          ...
		                                        'Quality',  beam_quality  );
		
		% Remove the energy field
		%   - It has incorrect time tags before v0.9.0
		edi_slow = rmfield(edi_slow, {'energy_gd12', 'energy_gd21'});
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
			edi_fast = mms_edi_l1b_efield_create( edi_fast_file, tstart, tend, ...
			                                        'CS_BCS',   true,          ...
			                                        'Quality',  beam_quality  );
		
			% Remove the energy field
			%   - It has incorrect time tags before v0.9.0
			edi_fast = rmfield(edi_fast, {'energy_gd12', 'energy_gd21'});
			
		% Not the error and continue processing slow survey data.
		catch ME
			mrfprintf('logerr', ME);
			edi_fast_file = '';
		end
	end
	
	% Check after we check for empty fast survey file.
	assert( slw_cnt + fst_cnt > 0, 'Unable to find fast or slow survey EDI files.' );

	% Combine slow and fast survey data
	edi = mms_edi_srvy_combine( edi_slow, edi_fast );
	clear edi_slow edi_fast

	% FGM data
	fg_l1b = mms_fg_read_l1b(fg_l1b_file, tstart, tend);
	fg_ql  = mms_fg_read_ql(fg_ql_file,   tstart, tend);

%------------------------------------%
% Compuate Average B                 %
%------------------------------------%
	% Time range that we have data
	if isempty(edi.tt2000_gd12) && isempty(edi.tt2000_gd12)
		error( 'No EDI E-field data available.' );
	elseif isempty(edi.tt2000_gd12)
		t0 = edi.tt2000_gd21(1);
		t1 = edi.tt2000_gd21(end);
	elseif isempty(edi.tt2000_gd21)
		t0 = edi.tt2000_gd12(1);
		t1 = edi.tt2000_gd12(end);
	else
		t0 = min( [edi.tt2000_gd12(1)   edi.tt2000_gd21(1)  ] );
		t1 = max( [edi.tt2000_gd12(end) edi.tt2000_gd21(end)] );
	end
	
	% Breakdown into time vectors
	tvec = MrCDF_Epoch_Breakdown( [t0, t1] );
	
	% Round down to the nearest DT seconds and recompute
	tvec(:,6)     = tvec(:,6) - mod(tvec(:,6), dt);
	tvec(:,7:end) = 0.0;
	tedge         = MrCDF_Epoch_Compute(tvec);
	t0            = tedge(1);
	t1            = tedge(2) + int64(dt * 1d9);

	% Find FGM data within this time interval
	%   - Prune out |B|
	ifg      = find( fg_l1b.tt2000 >= t0 & fg_l1b.tt2000 <= t1 );
	t_fg     = fg_l1b.tt2000(ifg);
	b_fg_bcs = fg_l1b.b_bcs(1:3, ifg);

	% Compute the averaged magnetic field
	b_avg_bcs = mms_edi_bavg(t_fg, b_fg_bcs, edi.tt2000_gd12, edi.tt2000_gd21, dt);

	% Clear unneeded data
	clear tvec tedge ifg t_fg b_fg_bcs

%------------------------------------%
% Beam Width                         %
%------------------------------------%
	% Number of beam hits per gun
	n12 = length( edi.tt2000_gd12 );
	n21 = length( edi.tt2000_gd21 );
	
	% Rotate interpolated B-field into EDI1 coordinate system
	bcs2edi1    = mms_instr_xxyz2instr('BCS',  'EDI1');
	b_gd12_edi1 = mrvector_rotate( bcs2edi1, b_avg_bcs.b_gd12 );
	b_gd21_edi1 = mrvector_rotate( bcs2edi1, b_avg_bcs.b_gd21 );
	
	% Combine data
	b_gdu_edi1 = [ b_gd12_edi1      b_gd21_edi1      ];
	fa_az      = [ edi.azimuth_gd12 edi.azimuth_gd21 ];
	fa_pol     = [ edi.polar_gd12   edi.polar_gd21   ];
	gun_id     = [ ones(1,n12)      ones(1,n21) + 1  ];
	
	% Beam width
	beam_width = mms_edi_beam_width( fa_az, fa_pol, b_gdu_edi1, gun_id );

	% Clear unneeded data
	clear bcs2edi1 b_gd12_edi1 b_gd21_edi1 fa_az fa_pol b_avg_bcs
	
%------------------------------------%
% Bavg in DMPA                       %
%------------------------------------%
	% Exctract data in time interval
	%   - Time stamps are always the same
	ifg       = find(fg_ql.tt2000 > t0 & fg_ql.tt2000 < t1);
	t_fg      = fg_ql.tt2000(ifg);
	b_fg_dmpa = fg_ql.b_dmpa(1:3,ifg);
	
	% Compute the average magnetic field
	b_avg_dmpa = mms_edi_bavg(t_fg, b_fg_dmpa, edi.tt2000_gd12, edi.tt2000_gd21, dt);
	
%------------------------------------%
% Rotate BCS to SMPA                 %
%------------------------------------%
	% BCS -> SMPA
	bcs2smpa = mms_fg_xbcs2smpa( zMPA );
	
	% BCS positions are not time dependent
	%   - Make them time-dependent by replicating them.
	%   - Added benefit of creating an empty array if no beams are in interval
	vg1_bcs = repmat( edi.virtual_gun1_bcs, 1, n12 );
	vg2_bcs = repmat( edi.virtual_gun2_bcs, 1, n21 );
	
	% Beam data
	%   - For scalar quantities, extract data.
	t_gd12       = edi.tt2000_gd12;
	t_gd21       = edi.tt2000_gd21;
	pos_vg1_smpa = mrvector_rotate( bcs2smpa, vg1_bcs );
	pos_vg2_smpa = mrvector_rotate( bcs2smpa, vg2_bcs );
	fv_gd12_smpa = mrvector_rotate( bcs2smpa, edi.fv_gd12_bcs );
	fv_gd21_smpa = mrvector_rotate( bcs2smpa, edi.fv_gd21_bcs );
	tof_gd12     = edi.tof_gd12;
	tof_gd21     = edi.tof_gd21;
	q_gd12       = edi.quality_gd12;
	q_gd21       = edi.quality_gd21;
	clear edi vg1_bcs vg2_bcs n12 n21
	
%------------------------------------%
% Despin SMPA to DMPA                %
%------------------------------------%
	% Despin using sunpulse times
	smpa2dmpa_gd12 = mms_dss_xdespin( sunpulse, t_gd12 );
	smpa2dmpa_gd21 = mms_dss_xdespin( sunpulse, t_gd21 );

	% Beam data
	%   - Clear data when finished
	pos_vg_dmpa = [ mrvector_rotate( smpa2dmpa_gd12, pos_vg1_smpa ), ...
	                mrvector_rotate( smpa2dmpa_gd21, pos_vg2_smpa ) ];
	fv_gdu_dmpa = [ mrvector_rotate( smpa2dmpa_gd12, fv_gd12_smpa ), ...
	                mrvector_rotate( smpa2dmpa_gd21, fv_gd21_smpa ) ];
	clear smpa2dmpa_gd12 smpa2dmpa_gd21 pos_vg1_smpa pos_vg2_dmpa fv_gd12_smpa fv_gd21_smpa
	
	% Averaged data
	%   - Clear data when finished
	t_avg      = b_avg_dmpa.t_avg;
	b_ave_dmpa = b_avg_dmpa.b_avg;
	b_std_dmpa = b_avg_dmpa.b_std;
	b_gdu_dmpa = [ b_avg_dmpa.b_gd12      b_avg_dmpa.b_gd21 ];
	recnum     = b_avg_dmpa.recnum;
	recnum_gdu = [ b_avg_dmpa.recnum_gd12 b_avg_dmpa.recnum_gd21 ];
	clear b_avg_dmpa
	
%------------------------------------%
% Compute E-Field                    %
%------------------------------------%

	% Cost Function
	[E_cf, v_ExB_cf, d_cf, d_delta_cf] = ...
		mms_edi_calc_efield_cf( b_ave_dmpa, pos_vg_dmpa, fv_gdu_dmpa, beam_width, recnum, recnum_gdu );

	% Beam convergence
	[E_bc, v_ExB_bc, d_bc, d_std_bc, q_bc] = ...
		mms_edi_calc_efield_bc( t_avg, b_ave_dmpa, pos_vg_dmpa, fv_gdu_dmpa, recnum, recnum_gdu, gun_id );

%------------------------------------%
% Gather Data                        %
%------------------------------------%
	% Parent files
	%   - Concatenate slow and fast survey file names
	%   - Remove filenames if any were empty
	%   - Strip off the directory.
	if iscell(att_file)
		parents = [ fg_l1b_file dss_file att_file edi_fast_file edi_slow_file ];
	else
		parents = { fg_l1b_file dss_file att_file edi_fast_file edi_slow_file };
	end
	iempty          = find( cellfun(@isempty, parents) );
	parents(iempty) = [];
	[~, name, ext]  = cellfun(@fileparts, parents, 'UniformOutput', false);
	parents         = strcat(name, ext);
	
	% Metadata structure
	meta = struct( 'sc',           sc,         ...
	               'instr',        instr,      ...
	               'mode',         'srvy',     ...
	               'optdesc',      optdesc,    ...
	               'directory',    sl_dir,     ...
	               'tstart',       tstart,     ...
	               'tend',         tend,       ...
	               'parents',      { parents } ...             % Prevent array of structures
	             );
	
	% Separate gun information
	tf_gd12 = gun_id == 1;
	tf_gd21 = gun_id == 2;
	clear gun_id
	
	% Averaged data
	b_interp_dmpa = struct( 't_avg',       t_avg,                            ...
	                        'dt_avg',      int64(dt*1e9),                    ...
	                        'b_avg',       single( b_ave_dmpa             ), ...
	                        'b_std',       single( b_std_dmpa             ), ...
	                        'b_gd12',      single( b_gdu_dmpa(:, tf_gd12) ), ...
	                        'b_gd21',      single( b_gdu_dmpa(:, tf_gd21) ), ...
	                        'recnum',      recnum,                           ...
	                        'recnum_gd12', recnum_gdu(tf_gd12),              ...
	                        'recnum_gd21', recnum_gdu(tf_gd21)               ...
	                      );
	clear b_ave_dmpa b_gdu_dmpa recnum recnum_gdu
	
	% E-field Cost Function
	efield_cf = struct( 'tt2000',  t_avg,                ...
	                    'E',       single( E_cf       ), ...
	                    'v_ExB',   single( v_ExB_cf   ), ...
	                    'd',       single( d_cf       ), ...
	                    'd_delta', single( d_delta_cf )  ...
	                  );
	clear E_cf v_ExB_cf d d_delta
	
	% E-field Beam Convergence
	efield_bc = struct( 'tt2000',  t_avg,              ...
	                    'E',       single( E_bc     ), ...
	                    'v_ExB',   single( v_ExB_bc ), ...
	                    'd',       single( d_bc     ), ...
	                    'd_std',   single( d_std_bc ), ...
	                    'quality', q_bc                ...
	                  );
	clear t_avg E_bc v_ExB_bc d_bc d_std_bc

	% Beam data
	beams = struct( 'tt2000_gd12',  t_gd12,                            ...
	                'tt2000_gd21',  t_gd21,                            ...
	                'pos_vg1_dmpa', single( pos_vg_dmpa(:, tf_gd12) ), ...
	                'pos_vg2_dmpa', single( pos_vg_dmpa(:, tf_gd21) ), ...
	                'fv_gd12_dmpa', single( fv_gdu_dmpa(:, tf_gd12) ), ...
	                'fv_gd21_dmpa', single( fv_gdu_dmpa(:, tf_gd21) ), ...
	                'tof_gd12',     tof_gd12,                          ...
	                'tof_gd21',     tof_gd21,                          ...
	                'quality_gd12', q_gd12,                            ...
	                'quality_gd21', q_gd21                             ...
	              );
	clear tf_gd12 tf_gd21 t_gd12 t_gd21 pos_vg_dmpa fv_gdu_dmpa q_gd12 q_gd21

%------------------------------------%
% Slow-Look                          %
%------------------------------------%
	edi_sl_file = mms_edi_sl_efield_write(meta, beams, efield_cf, efield_bc, b_interp_dmpa);

%------------------------------------%
% Quick-Look                         %
%------------------------------------%

	% Change output directory
	meta.directory = ql_dir;
	
	% Pick out quick-look data
	ql_data = struct( 'tt2000',     efield_cf.tt2000, ...
	                  'dt_avg',     b_interp_dmpa.dt_avg,  ...
	                  'E',          efield_cf.E,      ...
	                  'v_ExB',      efield_cf.v_ExB,  ...
	                  'E_bc',       efield_bc.E,      ...
	                  'v_ExB_bc',   efield_bc.v_ExB,  ...
	                  'quality_bc', efield_bc.quality ...
	                );
	clear efield_cf efield_bc beams b_interp_dmpa
	
	% Write the file
	edi_ql_file = mms_edi_ql_efield_write(meta, ql_data);
end
