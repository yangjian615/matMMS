%
% Name
%   mms_sdc_edi_ql_efield
%
% Purpose
%   Create a MATLAB save file of inputs needed for Bestarg.
%
% Calling Sequence
%   EDI_QL_FILE = mms_sdc_edi_ql_efield(SC, TSTART, TEND)
%     Compute the EDI electric field from MMS spacecraft SC during
%     the time interval [TSTART, TEND], and output a CDF file
%     named EDI_QL_FILE.
%
%   EDI_QL_FILE = mms_sdc_edi_ql_efield(..., 'ParamName', ParamValue)
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
%
function edi_ql_file = mms_sdc_edi_ql_efield(sc, tstart, tend, varargin)

% MMS2: May 9, 2015  16:08 - 16:13
% MMS4: May 6, 2015  15:30 - 15:35

	% Defaults
	sdc_root     = '/nfs/';
	save_dir     = '/nfs/edi/';
	attitude_dir = fullfile(sdc_root, 'ancillary', sc, 'defatt');
	hk_root      = fullfile(sdc_root, 'hk');
	beam_quality = 3;
	dt           = 5.0;

	% Optional arguments
	nOptArgs = length(varargin);
	for ii = 1 : 2 : nOptArgs
		switch varargin{ii}
			case 'SDCroot'
				sdc_root = varargin{ii+1};
			case 'SaveDir'
				save_dir = varargin{ii+1};
			case 'AttitudeDir'
				attitude_dir = varargin{ii+1};
			case 'HKdir'
				hk_root = varargin{ii+1};
			case 'BeamQuality'
				beam_quality = varargin{ii+1};
			case 'dt'
				dt = varargin{ii+1};
			otherwise
				error(['Unknown input parameter: "' varargin{ii} '".']);
		end
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
	edi = mms_edi_create_l2( edi_file, tstart, tend, ...
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
	t0 = min( [edi.tt2000_gd12(1)   edi.tt2000_gd21(1)  ] );
	t1 = max( [edi.tt2000_gd12(end) edi.tt2000_gd21(end)] );
	
	% Round down to the nearest 5-seconds
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
	E_dmpa = zeros( size(b_avg.b_avg), 'double');
	v_dmpa = zeros( size(b_avg.b_avg), 'double');
	d_dmpa = zeros( size(b_avg.b_avg), 'double');

	% Step through each interval
	for ii = 1 : size(b_avg.b_avg, 2)
		
		% B field data
		b_tt2000   = b_avg.t_avg(ii);
		b_avg_dmpa = b_avg.b_avg(1:3, ii);

		% GDU data that corresponds to the B field data: position and firing vectors
		iigd12_b_avgIntrp     = find( b_avg.inds_gd12 == ii );
		edi_gun1_virtual_dmpa = edi.virtual_gun1_dmpa(:, iigd12_b_avgIntrp);
		edi_gd12_fv_dmpa      = edi.fv_gd12_dmpa(:, iigd12_b_avgIntrp);
		edi_gd12_tof          = edi.tof_gd12 (iigd12_b_avgIntrp);

		iigd21_b_avgIntrp     = find( b_avg.inds_gd21 == ii );
		edi_gun2_virtual_dmpa = edi.virtual_gun2_dmpa(:, iigd21_b_avgIntrp);
		edi_gd21_fv_dmpa      = edi.fv_gd21_dmpa(:, iigd21_b_avgIntrp);
		edi_gd21_tof          = edi.tof_gd21(iigd21_b_avgIntrp);
		
		% Compute the electric field
		if (length (edi_gun1_virtual_dmpa) > 2) & (length (edi_gun2_virtual_dmpa) > 2)
			[E_temp, v_temp, d_temp] ...
				= edi_drift_step ( b_tt2000, ...
				                   b_avg_dmpa, ...
				                   edi_gun1_virtual_dmpa, ...
				                   edi_gd12_fv_dmpa, ...
				                   edi_gun2_virtual_dmpa, ...
				                   edi_gd21_fv_dmpa, ...
				                   edi_gd12_tof, ...
				                   edi_gd21_tof );
			
			E_dmpa(:, ii) = E_temp;
			v_dmpa(:, ii) = v_temp;
			d_dmpa(:, ii) = d_temp;
		else
			E_dmpa(:, ii) = [ -1e31; -1e31; -1e31 ];
			v_dmpa(:, ii) = [ -1e31; -1e31; -1e31 ];
			d_dmpa(:, ii) = [ -1e31; -1e31; -1e31 ];
		end
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
	                 'dt',           b_avg.dt_avg, ...
	                 'B_dmpa',       b_avg.b_avg', ...
	                 'B_std',        b_avg.b_std', ...
	                 'E_dmpa',       E_dmpa', ...
	                 'v_dmpa',       v_dmpa', ...
	                 'd_dmpa',       d_dmpa', ...
	                 'tt2000_gd12',  edi.tt2000_gd12', ...
	                 'tt2000_gd21',  edi.tt2000_gd21', ...
	                 'b_gd12_dmpa',  b_avg.b_gd12', ...
	                 'b_gd21_dmpa',  b_avg.b_gd21', ...
	                 'pos_vg1_dmpa', edi.virtual_gun1_dmpa', ...
	                 'pos_vg2_dmpa', edi.virtual_gun2_dmpa', ...
	                 'fv_gd12_dmpa', edi.fv_gd12_dmpa', ...
	                 'fv_gd21_dmpa', edi.fv_gd21_dmpa', ...
	                 'inds_gd12',    b_avg.inds_gd12', ...
	                 'inds_gd21',    b_avg.inds_gd12', ...
	                 'quality_gd12', edi.quality_gd12', ...
	                 'quality_gd21', edi.quality_gd21' ...
	               );

	% Clear data
	clear edi b_avg

	% Write data to a file
	file = mms_sdc_edi_write_ql_efield(edi_ql);
end