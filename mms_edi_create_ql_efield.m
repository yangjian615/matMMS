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
	attitude_dir    = fullfile(sdc_root, 'ancillary', sc, 'defatt');
	hk_root         = fullfile(sdc_root, 'hk');
	beam_quality    = 3;
	create_log_file = false;
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
		% Create a file name
		log_name = mms_construct_fname( sc, 'edi', 'slow', 'ql', ...
		                                'OptDesc',   'efield', ...
		                                'TStart',    tstart, ...
		                                'Version',   'v0.1.1' );
		% Change extension to 'log'
		[~, log_name] = fileparts(log_name)
		log_name      = fullfile( save_dir, [log_name '.log']);
		
		% Create the log object
		oLog = MrErrorLogger(log_name, 'Timestamp', true)
	end

%------------------------------------%
% Find Files                         %
%------------------------------------%

	% FG L1B Data File
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
	[att_file, count] = MrFile_Search( str, ...
	                                   'Closest',      true, ...
	                                   'TStart',       tstart, ...
	                                   'TEnd',         tend, ...
	                                   'TimeOrder',    '%Y%D', ...
	                                   'VersionRegex', 'V([0-9]{2})' );
	
	% EDI L1A E-Field Data File
	%    - Find last, so file descriptors are saved.
	instr   = 'edi';
	mode    = 'slow';
	level   = 'l1a';
	optdesc = 'efield';
	[edi_file, count, str] = mms_file_search(sc, instr, mode, level, ...
	                                         'TStart',    tstart, ...
	                                         'TEnd',      tend, ...
	                                         'OptDesc',   optdesc, ...
	                                         'SDCroot',   sdc_root);
	assert(count > 0, ['EDI file not found: "' str '".']);

	% Write names to log file
	if create_log_file
		oLog.AddText( 'Parent files:');
		oLog.AddText( strcat('  FG:      "', fg_file,  '"') );
		oLog.AddText( strcat('  DSS:     "', dss_file, '"') );
		oLog.AddText( strcat('  DefAtt:  "', att_file, '"') );
		oLog.AddText( strcat('  EDI:     "', edi_file, '"') );
		oLog.AddText( '' );
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
		zMPA = att_hdr.zMPA(:,1)';
	end
	
	% Sunpulse data
	sunpulse = mms_dss_read_sunpulse(dss_file, tstart, tend, 'UniquePulse', true);

	% EDI data
	edi = mms_edi_create_l2pre_efield( edi_file, tstart, tend, ...
	                                   'Attitude', defatt, ...
	                                   'CS_GSE',   false, ...
	                                   'CS_DMPA',  true, ...
	                                   'Quality',  beam_quality, ...
	                                   'Sunpulse', sunpulse, ...
	                                   'zMPA',     zMPA );

	% FGM data
	fg_ql = mms_fg_read_ql(fg_file, tstart, tend);

%------------------------------------%
% Prepare Data                       %
%------------------------------------%
	% Time range that we have data
	if isempty(edi.tt2000_gd12) && isempty(edi.tt2000_gd21)
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
	
	% Round down to the nearest 5 seconds and recompute
	tvec(:,6)     = tvec(:,6) - mod(tvec(:,6), 5);
	tvec(:,7:end) = 0.0;
	tedge         = MrCDF_Epoch_Compute(tvec);
	t0            = tedge(1);
	t1            = tedge(2) + int64(dt * 1d9);

	% Find FGM data within this time interval
	ifg = find(fg_ql.tt2000 >= t0 & fg_ql.tt2000 <= t1);
	fg_ql.tt2000 = fg_ql.tt2000(ifg);
	fg_ql.b_dmpa = fg_ql.b_dmpa(:, ifg);

	% Compute the averaged magnetic field
	b_avg = mms_edi_bavg(fg_ql.tt2000, fg_ql.b_dmpa, edi.tt2000_gd12, edi.tt2000_gd21, dt);

%------------------------------------%
% Compute Drift Step                 %
%------------------------------------%
	% Allocate memory
	E_dmpa     = zeros( size(b_avg.b_avg), 'double');
	v_dmpa     = zeros( size(b_avg.b_avg), 'double');
	d_dmpa     = zeros( size(b_avg.b_avg), 'double');
	d_std_dmpa = zeros( size(b_avg.b_avg), 'double');
	quality    = zeros( 1, size(b_avg.b_avg, 2), 'uint8');

	% Declare global variables for the drift step function
	global plot_beams;      plot_beams      = false;
	global rotation_method; rotation_method = 2;
	global use_v10502;      use_v10502      = false;

	% Step through each interval
	for ii = 1 : length(b_avg.recnum)
		recnum = b_avg.recnum(ii);
		
		% B field data
		B_tt2000 = b_avg.t_avg(ii);
		B_dmpa   = b_avg.b_avg(1:3, ii);

		% GDU data that corresponds to the B field data: position and firing vectors
		iigd12_b_avgIntrp = find( b_avg.recnum_gd12 == recnum );
		iigd21_b_avgIntrp = find( b_avg.recnum_gd21 == recnum );

		% Virtual gun positions
		gd_virtual_dmpa = [ edi.virtual_gun1_dmpa(:, iigd12_b_avgIntrp), ...
		                    edi.virtual_gun2_dmpa(:, iigd21_b_avgIntrp) ];
		
		% Firing Vectors
		gd_fv_dmpa = [ edi.fv_gd12_dmpa(:, iigd12_b_avgIntrp), ...
		               edi.fv_gd21_dmpa(:, iigd21_b_avgIntrp) ];
		
		% Time of Flight
		gd_tof = [ edi.tof_gd12(iigd12_b_avgIntrp), ...
		           edi.tof_gd21(iigd21_b_avgIntrp) ];
		
		% Gun ID of each data point
		n_gd12                = size( iigd12_b_avgIntrp, 2 );
		gd_ID                 = ones( 1, size(gd_virtual_dmpa, 2) );
		gd_ID( n_gd12+1:end ) = 2;

		% Compute the electric field
		if size(gd_virtual_dmpa, 2) > 2
			[ d_temp d_std_temp v_temp E_temp q_temp] ...
				= edi_drift_step( '', ...
				                  B_tt2000, ...
				                  B_dmpa, ...
				                  gd_virtual_dmpa, ...
				                  gd_fv_dmpa, ...
				                  gd_ID );

			% Replace NaN with fill value
			if isnan(E_temp(1))
				E_temp     = [ -1e31; -1e31; -1e31 ];
				v_temp     = [ -1e31; -1e31; -1e31 ];
				d_temp     = [ -1e31; -1e31; -1e31 ];
				d_std_temp = [ -1e31; -1e31; -1e31 ];
				q_temp     = 0;
			end
		else
			E_temp     = [ -1e31; -1e31; -1e31 ];
			v_temp     = [ -1e31; -1e31; -1e31 ];
			d_temp     = [ -1e31; -1e31; -1e31 ];
			d_std_temp = [ -1e31; -1e31; -1e31 ];
			q_temp     = 0;
		end

		% Store results
		E_dmpa(:, ii)     = E_temp;
		v_dmpa(:, ii)     = v_temp;
		d_dmpa(:, ii)     = d_temp;
		d_std_dmpa(:, ii) = d_std_temp;
		quality(ii)       = q_temp;
	end

%------------------------------------%
% Gather Data                        %
%------------------------------------%

	% Gather data to be written
	%   - Record varying dimension must be first
	edi_ql = struct( 'sc',           sc, ...
	                 'instr',        instr, ...
	                 'mode',         mode, ...
	                 'optdesc',      optdesc, ...
	                 'directory',    save_dir, ...
	                 'tstart',       tstart, ...
	                 'tend',         tend, ...
	                 'parents',      { { fg_file, dss_file, att_file{:}, edi_file } }, ...
	                 'tt2000',       b_avg.t_avg', ...
	                 'dt',           int64( b_avg.dt_avg ) * int64(1e9), ... % s -> ns
	                 'B_dmpa',       b_avg.b_avg', ...
	                 'B_std_dmpa',   b_avg.b_std', ...
	                 'E_dmpa',       E_dmpa', ...
	                 'v_dmpa',       v_dmpa', ...
	                 'd_dmpa',       d_dmpa', ...
	                 'quality',      quality', ...
	                 'd_std_dmpa',   d_std_dmpa', ...
	                 'tt2000_gd12',  edi.tt2000_gd12', ...
	                 'tt2000_gd21',  edi.tt2000_gd21', ...
	                 'b_gd12_dmpa',  b_avg.b_gd12', ...
	                 'b_gd21_dmpa',  b_avg.b_gd21', ...
	                 'pos_vg1_dmpa', edi.virtual_gun1_dmpa', ...
	                 'pos_vg2_dmpa', edi.virtual_gun2_dmpa', ...
	                 'fv_gd12_dmpa', edi.fv_gd12_dmpa', ...
	                 'fv_gd21_dmpa', edi.fv_gd21_dmpa', ...
	                 'recnum',       b_avg.recnum', ...
	                 'recnum_gd12',  b_avg.recnum_gd12', ...
	                 'recnum_gd21',  b_avg.recnum_gd21', ...
	                 'quality_gd12', edi.quality_gd12', ...
	                 'quality_gd21', edi.quality_gd21' ...
	               );

	% Clear data
	clear edi b_avg

	% Write data to a file
	edi_ql_file = mms_edi_write_ql_efield(edi_ql);
end