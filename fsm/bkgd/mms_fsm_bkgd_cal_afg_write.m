%
% Name
%   mms_fsm_l2plus_cal_afg_write
%
% Purpose
%   Create a MATLAB save file of inputs needed for Bestarg.
%
% Calling Sequence
%   EDI_QL_FILE = mms_edi_ql_efield_write(EDI_QL)
%     Write EDI quick-look data constained in the structure EDI_QL
%     and created by mms_edi_create_ql_efield.m to a CDF file named
%     EDI_QL_FILE.
%
% Parameters
%   EDI_QL:         in, required, type=string
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
function afg_cal_file = mms_fsm_l2plus_cal_afg_write(sc, mode, tstart, afg_bkgd, varargin)

	% Global variables
	%   - See mms_fsm_init.m
	global cal_path_root data_path_root dropbox_root hk_root log_path_root unh_data_root

%------------------------------------%
% Check Inputs                       %
%------------------------------------%

	% Defaults
	optdesc  = 'cal-afg';
	tf_brst  = false;
	tf_empty = false;
	parents  = {};

	% Check inputs
	nOptArgs = length(varargin);
	for ii = 1 : 2 : nOptArgs
		switch varargin{ii}
			case 'OptDesc'
				optdesc = varargin{ii+1};
			case 'Burst'
				tf_brst = varargin{ii+1};
			case 'EmptyFile'
				tf_empty = varargin{ii+1};
			case 'Parents'
				parents = varargin{ii+1};
			otherwise
				error(['Optional argument not recognized: "' varargin{ii} '".']);
		end
	end
	
	% Constants
	instr   = 'fsm';
	level   = 'l2plus';

%------------------------------------%
% Verify Data                        %
%------------------------------------%

	%
	% Check sizes
	%
	assert( isa(afg_bkgd.t,     'int64'),  'afg_bkgd.t should be int64.' );
	assert( isa(afg_bkgd.f,     'single'), 'afg_bkgd.f should be single.');
	assert( isa(afg_bkgd.amp,   'single'), 'afg_bkgd.amp should be single.');
	assert( isa(afg_bkgd.phase, 'single'), 'afg_bkgd.phase should be single.');
	assert( isa(afg_bkgd.psd,   'single'), 'afg_bkgd.psd should be single.');
	assert( isa(afg_bkgd.flag,  'uint8'),  'afg_bkgd.flag should be uint8.');
	
	% convert the start time to yyyy-mm-dd
%	tstart = MrTimeParser(tstart, '%Y-%M-%dT%H:%m:%S', '%Y%M%d');

%------------------------------------%
% Version History                    %
%------------------------------------%
	
	% Describe the modifications to each version
	mods    = {  'v0.0.0 -- First version.' ...
	          };
	
	version = regexp( mods{end}, '^v([0-9]+\.[0-9]+\.[0-9]+)', 'tokens' );
	version = version{1}{1};
%	version = regexp( mods{end}, '^v([0-9])+\.([0-9])+\.([0-9]+)', 'tokens' );
%	vx      = version{1}{1};
%	vy      = version{1}{2};
%	vz      = version{1}{3};

%------------------------------------%
% Create Output File Name            %
%------------------------------------%
	
	% Create the output filename
	afg_cal_file = mms_construct_filename( sc, instr, mode, level,    ...
	                                       'Directory', dropbox_root, ...
	                                       'OptDesc',   optdesc,      ...
	                                       'TStart',    tstart,       ...
	                                       'Version',   version );

%------------------------------------------------------
% Global Attributes                                   |
%------------------------------------------------------
	if isempty(optdesc)
		data_type      = [mode '_' level];
		logical_source = [instr '_' mode '_' level];
	else
		data_type      = [mode '_' level '_' optdesc];
		logical_source = [instr '_' mode '_' level '_' optdesc];
	end
	[~, logical_file_id, ext] = fileparts(afg_cal_file);
	source_name = ['MMS' sc(4) '>MMS Satellite Number ' sc(4)];

	% Define global attributes
	global_attrs = struct( 'Data_type',                  data_type, ...
	                       'Data_version',               version, ...
	                       'Descriptor',                 'FSM', ...
	                       'Discipline',                 'Space Physics>Magnetospheric Science', ...
	                       'File_naming_convention',     'source_descriptor_datatype_yyyyMMdd', ...
	                       'Generation_date',            datestr(now(), 'yyyymmdd'), ...
	                       'Instrument_type',            'Magnetic Fields (space)', ...
	                       'Logical_file_id',            logical_file_id, ...
	                       'Logical_source',             logical_source, ...
	                       'Logical_source_description', 'Level 2-Plus Fluxgate-Searchcoil Merged Magnetic Field Data', ...
	                       'Mission_group',              'MMS', ...
	                       'PI_affiliation',             'UNH, LPP, IWF, UCLA, UCLA, IWF, UNH', ...
	                       'PI_name',                    ['R.B. Torbert, O. LeContel, W. Magnes, C.T. Russell, R. Strangeway, ', ...
	                                                      'D. Fischer, M.R. Argall.'], ...
	                       'Project',                    'STP>Solar Terrestrial Physics', ...
	                       'Source_name',                source_name, ...
	                       'TEXT',                       ['AFG calibration data used while created the FSM data product. ' ...
	                                                      'A description of the method and of the data can be found at: ' ...
	                                                      ''], ...
	                       'HTTP_LINK',                  { {'http://mms-fields.unh.edu/' ...
	                                                        'http://mms.gsfc.nasa.gov/index.html'} }, ...
	                       'LINK_TEXT',                  { {'UNH FIELDS Home Page', ...
	                                                        'NASA MMS Home'} }, ...
	                       'MODS',                       { mods }, ...
	                       'Acknowledgements',           ' ', ...
	                       'Generated_by',               'UNH, IWF, LPP, UCLA', ...
	                       'Parents',                    { parents }, ...
	                       'Skeleton_version',           ' ', ...
	                       'Rules_of_use',               ' ', ...
	                       'Time_resolution',            ' '  ...
	                     );

%------------------------------------------------------
% Variables                                           |
%------------------------------------------------------
	% Variable naming convention
	%   scId_instrumentId_paramName[_coordSys][_paramQualifier][_subModeLevel][_mode][_level]
	prefix  = [sc '_' instr '_'];
	suffix  = ['_' mode '_' level];
	
	% Spectra
	t_vname         = 'Epoch';
	f_vname         = [prefix 'f'                suffix];
	c_vname         = [prefix 'component'        suffix];
	afg_amp_vname   = [prefix 'amp'   '_omb_afg' suffix];
	afg_phase_vname = [prefix 'phase' '_omb_afg' suffix];
	afg_psd_vname   = [prefix 'psd'   '_omb_afg' suffix];
	flag_vname      = [prefix 'flag'             suffix];
	flag_hist_vname = [prefix 'flag'  '_hist'    suffix];
	
	% AFG Histogram
	afg_hamp_vname   = [prefix 'amp'   '_hist_afg'  suffix];
	afg_hphase_vname = [prefix 'phase' '_hist_afg'  suffix];
	afg_hpsd_vname   = [prefix 'psd'   '_hist_afg'  suffix];
	afg_famp_vname   = [prefix 'amp'   '_floor_afg' suffix];
	afg_fphase_vname = [prefix 'phase' '_floor_afg' suffix];
	afg_fpsd_vname   = [prefix 'psd'   '_floor_afg' suffix];
	afg_bamp_vname   = [prefix 'amp'   '_bins_afg'  suffix];
	afg_bphase_vname = [prefix 'phase' '_bins_afg'  suffix];
	afg_bpsd_vname   = [prefix 'psd'   '_bins_afg'  suffix];
	
	
	f_labl_vname    = [prefix 'freq'  '_labels'  suffix];
	c_labl_vname    = [prefix 'xyz'   '_labels'  suffix];

	% Variable Data
	%   - Order as [ DEPEND_1,  DEPEND_2,  DEPEND_0 ]
	%   - Order as [ COMPONENT, FREQUENCY, TIME     ]
	var_list = { t_vname,         afg_bkgd.t',        ...
	             f_vname,         afg_bkgd.f',        ...
	             c_vname,         uint8([0 1 2]),      ...
	             afg_amp_vname,   afg_bkgd.amp,   ...
	             afg_phase_vname, afg_bkgd.phase, ...
	             afg_psd_vname,   afg_bkgd.psd,   ...
	             flag_vname,      afg_bkgd.flag',     ...
	             f_labl_vname,   'Frequency',          ...
	             c_labl_vname,   {'X', 'Y', 'Z'},      ...
	             flag_hist_vname,  afg_bkgd.hist_flag,   ...
	             afg_hamp_vname,   afg_bkgd.amp_hist,    ...
	             afg_hphase_vname, afg_bkgd.phase_hist,  ...
	             afg_hpsd_vname,   afg_bkgd.psd_hist,    ...
	             afg_famp_vname,   afg_bkgd.amp_floor,   ...
	             afg_fphase_vname, afg_bkgd.phase_floor, ...
	             afg_fpsd_vname,   afg_bkgd.psd_floor,   ...
	             afg_bamp_vname,   afg_bkgd.amp_bins,    ...
	             afg_bphase_vname, afg_bkgd.phase_bins,  ...
	             afg_bpsd_vname,   afg_bkgd.psd_bins     ...
	           };
	clear afg_bkgd

	% Record Variance
	recbound = { t_vname,         ...
	             afg_amp_vname,   ...
	             afg_phase_vname, ...
	             afg_psd_vname,   ...
	             flag_vname,      ...
	             flag_hist_vname,  ...
	             afg_hamp_vname,   ...
	             afg_hphase_vname, ...
	             afg_hpsd_vname,   ...
	             afg_famp_vname,   ...
	             afg_fphase_vname, ...
	             afg_fpsd_vname    ...
	           };

	% Data types
	vardatatypes = { t_vname,           'cdf_time_tt2000', ...
	                 f_vname,           'cdf_float', ...
	                 c_vname,           'cdf_uint1',  ...
	                 afg_amp_vname,     'cdf_float', ...
	                 afg_phase_vname,   'cdf_float', ...
	                 afg_psd_vname,     'cdf_float', ...
	                 flag_vname,        'cdf_uint1', ...
	                 f_labl_vname,      'cdf_char',  ...
	                 c_labl_vname,      'cdf_char',  ...
	                 flag_hist_vname,   'cdf_uint1', ...
	                 afg_hamp_vname,    'cdf_uint4', ...
	                 afg_hphase_vname,  'cdf_uint4', ...
	                 afg_hpsd_vname,    'cdf_uint4', ...
	                 afg_famp_vname,    'cdf_float', ...
	                 afg_fphase_vname,  'cdf_float', ...
	                 afg_fpsd_vname,    'cdf_float', ...
	                 afg_bamp_vname,    'cdf_float', ...
	                 afg_bphase_vname,  'cdf_float', ...
	                 afg_bpsd_vname,    'cdf_float'  ...
	               };
	
	% CDF Compression Level
	varcompress = { afg_amp_vname,     'gzip.6', ...
	                afg_phase_vname,   'gzip.6', ...
	                afg_psd_vname,     'gzip.6', ...
	                flag_vname,        'gzip.6',  ...
	                flag_hist_vname,   'gzip.6', ...
	                afg_hamp_vname,    'gzip.6', ...
	                afg_hphase_vname,  'gzip.6', ...
	                afg_hpsd_vname,    'gzip.6', ...
	                afg_famp_vname,    'gzip.6', ...
	                afg_fphase_vname,  'gzip.6', ...
	                afg_fpsd_vname,    'gzip.6', ...
	                afg_bamp_vname,    'gzip.6', ...
	                afg_bphase_vname,  'gzip.6', ...
	                afg_bpsd_vname,    'gzip.6'  ...
	              };
	

%------------------------------------------------------
% Variable Attributes                                 |
%------------------------------------------------------
	catdesc       = {};
	depend_0      = {};
	depend_1      = {};
	depend_2      = {};
	depend_3      = {};
	display_type  = {};
	fieldnam      = {};
	fillval       = {};
	format        = {};
	lablaxis      = {};
	labl_ptr_1    = {};
	labl_ptr_2    = {};
	si_conversion = {};
	time_base     = {};
	units         = {};
	validmin      = {};
	validmax      = {};
	var_type      = {};
	
	% EPOCH
	catdesc       = [ catdesc,       t_vname, 'Time tags for the merged magnetic field data.' ];
	fieldnam      = [ fieldnam,      t_vname, 'Time' ];
	fillval       = [ fillval,       t_vname, spdfcomputett2000([9999 12 31 23 59 59 999 999 999]) ];
	format        = [ format,        t_vname, 'I16' ];
	lablaxis      = [ lablaxis,      t_vname, 'TAI' ];
	si_conversion = [ si_conversion, t_vname, '1e-9>s' ];
	time_base     = [ time_base,     t_vname, 'J2000' ];
	units         = [ units,         t_vname, 'ns' ];
	validmin      = [ validmin,      t_vname, spdfcomputett2000([2015 3 12 0 0 0 0 0 0]) ];
	validmax      = [ validmax,      t_vname, spdfcomputett2000([2050 1  1 0 0 0 0 0 0]) ];
	var_type      = [ var_type,      t_vname, 'support_data' ];
	
	% Frequency
	catdesc       = [ catdesc,       f_vname, 'Center frequency bins.' ];
	fieldnam      = [ fieldnam,      f_vname, 'Frequency' ];
	fillval       = [ fillval,       f_vname, single(-1e31) ];
	format        = [ format,        f_vname, 'f6.1' ];
	lablaxis      = [ lablaxis,      f_vname, 'f' ];
	si_conversion = [ si_conversion, f_vname, '1e0>Hz' ];
	units         = [ units,         f_vname, 'Hz' ];
	validmin      = [ validmin,      f_vname, single(0.0) ];
	validmax      = [ validmax,      f_vname, single(1024.0) ];
	var_type      = [ var_type,      f_vname, 'support_data' ];
	
	% Component
	catdesc       = [ catdesc,       c_vname, 'Component index.' ];
	fieldnam      = [ fieldnam,      c_vname, 'Component' ];
	fillval       = [ fillval,       c_vname, uint8(255) ];
	format        = [ format,        c_vname, 'I1' ];
	lablaxis      = [ lablaxis,      c_vname, 'Component' ];
	validmin      = [ validmin,      c_vname, uint8(0) ];
	validmax      = [ validmax,      c_vname, uint8(2) ];
	var_type      = [ var_type,      c_vname, 'support_data' ];
	
	% AMP_afg
	catdesc       = [ catdesc,       afg_amp_vname, 'Wave amplitude of the FGM magnetic field.' ];
	depend_0      = [ depend_0,      afg_amp_vname, t_vname ];
	depend_1      = [ depend_1,      afg_amp_vname, f_vname ];
	depend_2      = [ depend_2,      afg_amp_vname, c_vname ];
	display_type  = [ display_type,  afg_amp_vname, 'spectrogram'];
	fieldnam      = [ fieldnam,      afg_amp_vname, 'Amplitude' ];
	fillval       = [ fillval,       afg_amp_vname, single(-1e31) ];
	format        = [ format,        afg_amp_vname, 'f11.4' ];
	lablaxis      = [ lablaxis,      afg_amp_vname, 'Amplitude'];
%	labl_ptr_1    = [ labl_ptr_1,    afg_amp_vname, f_labl_vname ];
%	labl_ptr_2    = [ labl_ptr_2,    afg_amp_vname, c_labl_vname ];
	si_conversion = [ si_conversion, afg_amp_vname, '1e-9>T' ];
	units         = [ units,         afg_amp_vname, 'nT' ];
	validmin      = [ validmin,      afg_amp_vname, single(0.0) ];
	validmax      = [ validmax,      afg_amp_vname, single(1e5) ];
	var_type      = [ var_type,      afg_amp_vname, 'data' ];
	
	% PHASE_afg
	catdesc       = [ catdesc,       afg_phase_vname, 'Wave phase of the FGM magnetic field.' ];
	depend_0      = [ depend_0,      afg_phase_vname, t_vname ];
	depend_1      = [ depend_1,      afg_phase_vname, f_vname ];
	display_type  = [ display_type,  afg_phase_vname, 'spectrogram'];
	fieldnam      = [ fieldnam,      afg_phase_vname, 'Phase' ];
	fillval       = [ fillval,       afg_phase_vname, single(-1e31) ];
	format        = [ format,        afg_phase_vname, 'f11.4' ];
	labl_ptr_1    = [ labl_ptr_1,    afg_phase_vname, c_labl_vname ];
	si_conversion = [ si_conversion, afg_phase_vname, '1e-9>T' ];
	units         = [ units,         afg_phase_vname, 'degrees' ];
	validmin      = [ validmin,      afg_phase_vname, single(0.0) ];
	validmax      = [ validmax,      afg_phase_vname, single(360.0) ];
	var_type      = [ var_type,      afg_phase_vname, 'data' ];
	
	% PSD_afg
	catdesc       = [ catdesc,       afg_psd_vname, 'Power spectral density of the FGM magnetic field.' ];
	depend_0      = [ depend_0,      afg_psd_vname, t_vname ];
	depend_1      = [ depend_1,      afg_psd_vname, f_vname ];
	display_type  = [ display_type,  afg_psd_vname, 'spectrogram'];
	fieldnam      = [ fieldnam,      afg_psd_vname, 'Phase' ];
	fillval       = [ fillval,       afg_psd_vname, single(-1e31) ];
	format        = [ format,        afg_psd_vname, 'f11.4' ];
	labl_ptr_1    = [ labl_ptr_1,    afg_psd_vname, c_labl_vname ];
	si_conversion = [ si_conversion, afg_psd_vname, '1e-9>T' ];
	units         = [ units,         afg_psd_vname, 'nT^2/Hz' ];
	validmin      = [ validmin,      afg_psd_vname, single(0.0) ];
	validmax      = [ validmax,      afg_psd_vname, single(1e4) ];
	var_type      = [ var_type,      afg_psd_vname, 'data' ];
	
	% FLAG
	catdesc       = [ catdesc,       flag_vname, 'Data flag' ];
	depend_0      = [ depend_0,      flag_vname, t_vname ];
	display_type  = [ display_type,  flag_vname, 'time_series'];
	fieldnam      = [ fieldnam,      flag_vname, 'Flag' ];
	fillval       = [ fillval,       flag_vname, uint8(255) ];
	format        = [ format,        flag_vname, 'I2' ];
	lablaxis      = [ lablaxis,      flag_vname, 'Phase' ];
	validmin      = [ validmin,      flag_vname, uint8(0) ];
	validmax      = [ validmax,      flag_vname, uint8(20) ];
	var_type      = [ var_type,      flag_vname, 'data' ];
	
	% COMPONENT LABEL
	catdesc       = [ catdesc,       c_labl_vname, 'Axis labels for vector components' ];
	fieldnam      = [ fieldnam,      c_labl_vname, 'Component labels' ];
	format        = [ format,        c_labl_vname, 'A1' ];
	var_type      = [ var_type,      c_labl_vname, 'metadata' ];
	
	% FREQUENCY LABEL
	catdesc       = [ catdesc,       f_labl_vname, 'Axis label for frequencies.' ];
	fieldnam      = [ fieldnam,      f_labl_vname, 'Frequency labels' ];
	format        = [ format,        f_labl_vname, 'A9' ];
	var_type      = [ var_type,      f_labl_vname, 'metadata' ];
	
%-------------------------------------------------------------------------
% FGM HISTOGRAM INFO /////////////////////////////////////////////////////
%-------------------------------------------------------------------------
	
	% HIST_AMP_afg
	catdesc       = [ catdesc,       afg_hamp_vname, 'FGM wave amplitude, Histogrammed in time at each frequency.' ];
	depend_0      = [ depend_0,      afg_hamp_vname, f_vname ];
	depend_1      = [ depend_1,      afg_hamp_vname, afg_bamp_vname ];
	depend_2      = [ depend_2,      afg_hamp_vname, flag_hist_vname ];
	depend_3      = [ depend_3,      afg_hamp_vname, c_vname ];
	display_type  = [ display_type,  afg_hamp_vname, 'spectrogram'];
	fieldnam      = [ fieldnam,      afg_hamp_vname, 'Histogrammed wave amplitude (FGM)' ];
	fillval       = [ fillval,       afg_hamp_vname, uint32(65535) ];
	format        = [ format,        afg_hamp_vname, 'f11.4' ];
	labl_ptr_1    = [ labl_ptr_1,    afg_hamp_vname, c_labl_vname ];
	si_conversion = [ si_conversion, afg_hamp_vname, '1e-9>T' ];
	units         = [ units,         afg_hamp_vname, 'degrees' ];
	validmin      = [ validmin,      afg_hamp_vname, uint32(0) ];
	validmax      = [ validmax,      afg_hamp_vname, uint32(65533) ];
	var_type      = [ var_type,      afg_hamp_vname, 'data' ];
	
	% HIST_PHASE_afg
	catdesc       = [ catdesc,       afg_hphase_vname, 'FGM wave phase, Histogrammed in time at each frequency.' ];
	depend_0      = [ depend_0,      afg_hphase_vname, f_vname ];
	depend_1      = [ depend_1,      afg_hphase_vname, afg_bphase_vname ];
	depend_2      = [ depend_2,      afg_hphase_vname, flag_hist_vname ];
	depend_3      = [ depend_3,      afg_hphase_vname, c_vname ];
	display_type  = [ display_type,  afg_hphase_vname, 'spectrogram'];
	fieldnam      = [ fieldnam,      afg_hphase_vname, 'Histogrammed wave phase (FGM)' ];
	fillval       = [ fillval,       afg_hphase_vname, uint32(65535) ];
	format        = [ format,        afg_hphase_vname, 'f11.4' ];
	labl_ptr_1    = [ labl_ptr_1,    afg_hphase_vname, c_labl_vname ];
	si_conversion = [ si_conversion, afg_hphase_vname, '1e-9>T' ];
	units         = [ units,         afg_hphase_vname, 'degrees' ];
	validmin      = [ validmin,      afg_hphase_vname, uint32(0) ];
	validmax      = [ validmax,      afg_hphase_vname, uint32(65533) ];
	var_type      = [ var_type,      afg_hphase_vname, 'data' ];
	
	% HIST_PSD_afg
	catdesc       = [ catdesc,       afg_hpsd_vname, 'FGM wave power, Histogrammed in time at each frequency.' ];
	depend_0      = [ depend_0,      afg_hpsd_vname, f_vname ];
	depend_1      = [ depend_1,      afg_hpsd_vname, afg_bpsd_vname ];
	depend_2      = [ depend_2,      afg_hpsd_vname, flag_hist_vname ];
	depend_3      = [ depend_3,      afg_hpsd_vname, c_vname ];
	display_type  = [ display_type,  afg_hpsd_vname, 'spectrogram'];
	fieldnam      = [ fieldnam,      afg_hpsd_vname, 'Histogrammed wave power (FGM)' ];
	fillval       = [ fillval,       afg_hpsd_vname, uint32(65535) ];
	format        = [ format,        afg_hpsd_vname, 'I5' ];
	labl_ptr_1    = [ labl_ptr_1,    afg_hpsd_vname, c_labl_vname ];
	si_conversion = [ si_conversion, afg_hpsd_vname, '1e-9>T' ];
	units         = [ units,         afg_hpsd_vname, 'degrees' ];
	validmin      = [ validmin,      afg_hpsd_vname, uint32(0) ];
	validmax      = [ validmax,      afg_hpsd_vname, uint32(65533) ];
	var_type      = [ var_type,      afg_hpsd_vname, 'data' ];
	
	
	% FLOOR_AMP_afg
	catdesc       = [ catdesc,       afg_famp_vname, ['Log base 10 of the most common wave amplitude, obtained by histogramming ' ...
	                                                  'wave amplitude at each frequency over the course of a day.'] ];
	depend_0      = [ depend_0,      afg_famp_vname, afg_bamp_vname ];
	depend_1      = [ depend_1,      afg_famp_vname, flag_hist_vname ];
	display_type  = [ display_type,  afg_famp_vname, 'time_series'];
	fieldnam      = [ fieldnam,      afg_famp_vname, 'Wave amplitude' ];
	fillval       = [ fillval,       afg_famp_vname, single(-1e31) ];
	format        = [ format,        afg_famp_vname, 'f11.4' ];
	labl_ptr_1    = [ labl_ptr_1,    afg_famp_vname, c_labl_vname ];
	si_conversion = [ si_conversion, afg_famp_vname, '1e-9>T' ];
	units         = [ units,         afg_famp_vname, 'degrees' ];
	validmin      = [ validmin,      afg_famp_vname, single(-10.0) ];
	validmax      = [ validmax,      afg_famp_vname, single(10.0) ];
	var_type      = [ var_type,      afg_famp_vname, 'data' ];
	
	% FLOOR_PHASE_afg
	catdesc       = [ catdesc,       afg_fphase_vname, ['The most common wave phase, obtained by histogramming ' ...
	                                                    'wave phase at each frequency over the course of a day.'] ];
	depend_0      = [ depend_0,      afg_fphase_vname, afg_bphase_vname ];
	depend_1      = [ depend_1,      afg_fphase_vname, flag_hist_vname ];
	depend_2      = [ depend_2,      afg_fphase_vname, c_vname ];
	display_type  = [ display_type,  afg_fphase_vname, 'time_series'];
	fieldnam      = [ fieldnam,      afg_fphase_vname, 'Wave phase' ];
	fillval       = [ fillval,       afg_fphase_vname, single(-1e31) ];
	format        = [ format,        afg_fphase_vname, 'f11.4' ];
	labl_ptr_1    = [ labl_ptr_1,    afg_fphase_vname, c_labl_vname ];
	si_conversion = [ si_conversion, afg_fphase_vname, '1e0>degrees' ];
	units         = [ units,         afg_fphase_vname, 'degrees' ];
	validmin      = [ validmin,      afg_fphase_vname, single(-180.0) ];
	validmax      = [ validmax,      afg_fphase_vname, single(180.0) ];
	var_type      = [ var_type,      afg_fphase_vname, 'data' ];
	
	
	% FLOOR_PSD_afg
	catdesc       = [ catdesc,       afg_fpsd_vname, ['Log base 10 of the most common wave power, obtained by histogramming ' ...
	                                                  'wave power at each frequency over the course of a day.'] ];
	depend_0      = [ depend_0,      afg_fpsd_vname, afg_bpsd_vname ];
	depend_1      = [ depend_1,      afg_fpsd_vname, flag_hist_vname ];
	depend_2      = [ depend_2,      afg_fpsd_vname, c_vname ];
	display_type  = [ display_type,  afg_fpsd_vname, 'time_series'];
	fieldnam      = [ fieldnam,      afg_fpsd_vname, 'Power spectral density' ];
	fillval       = [ fillval,       afg_fpsd_vname, single(-1e31) ];
	format        = [ format,        afg_fpsd_vname, 'f11.4' ];
	labl_ptr_1    = [ labl_ptr_1,    afg_fpsd_vname, c_labl_vname ];
	si_conversion = [ si_conversion, afg_fpsd_vname, '1e-18>T^2 Hz^-1' ];
	units         = [ units,         afg_fpsd_vname, 'nT^2 Hz^-1' ];
	validmin      = [ validmin,      afg_fpsd_vname, single(-15.0) ];
	validmax      = [ validmax,      afg_fpsd_vname, single(15.0) ];
	var_type      = [ var_type,      afg_fpsd_vname, 'data' ];
	
	
	% BINS_AMP_afg
	catdesc       = [ catdesc,       afg_bamp_vname, 'Histogram bins for wave amplitude.' ];
	fieldnam      = [ fieldnam,      afg_bamp_vname, 'Wave amplitude bins' ];
	fillval       = [ fillval,       afg_bamp_vname, single(-1e31) ];
	format        = [ format,        afg_bamp_vname, 'f11.4' ];
	si_conversion = [ si_conversion, afg_bamp_vname, '1e-9>T' ];
	units         = [ units,         afg_bamp_vname, 'nT' ];
	validmin      = [ validmin,      afg_bamp_vname, single(-10.0) ];
	validmax      = [ validmax,      afg_bamp_vname, single(10.0) ];
	var_type      = [ var_type,      afg_bamp_vname, 'support_data' ];
	
	% BINS_PHASE_afg
	catdesc       = [ catdesc,       afg_bphase_vname, 'Histogram bins for wave phase.' ];
	fieldnam      = [ fieldnam,      afg_bphase_vname, 'Wave amplitude bins' ];
	fillval       = [ fillval,       afg_bphase_vname, single(-1e31) ];
	format        = [ format,        afg_bphase_vname, 'f11.4' ];
	si_conversion = [ si_conversion, afg_bphase_vname, '1e0>degrees' ];
	units         = [ units,         afg_bphase_vname, 'degrees' ];
	validmin      = [ validmin,      afg_bphase_vname, single(-180.0) ];
	validmax      = [ validmax,      afg_bphase_vname, single(180.0) ];
	var_type      = [ var_type,      afg_bphase_vname, 'support_data' ];
	
	% BINS_PSD_afg
	catdesc       = [ catdesc,       afg_bpsd_vname, 'Histogram bins for wave amplitude.' ];
	fieldnam      = [ fieldnam,      afg_bpsd_vname, 'Wave amplitude bins' ];
	fillval       = [ fillval,       afg_bpsd_vname, single(-1e31) ];
	format        = [ format,        afg_bpsd_vname, 'f11.4' ];
	si_conversion = [ si_conversion, afg_bpsd_vname, '1e-18>T^2 Hz^-1' ];
	units         = [ units,         afg_bpsd_vname, 'nT^2 Hz^-1' ];
	validmin      = [ validmin,      afg_bpsd_vname, single(-10.0) ];
	validmax      = [ validmax,      afg_bpsd_vname, single(10.0) ];
	var_type      = [ var_type,      afg_bpsd_vname, 'support_data' ];
	
	
%-------------------------------------------------------------------------
% Collect Variable Attributes ////////////////////////////////////////////
%-------------------------------------------------------------------------
	% Collect into a structure
	var_attrs = struct( 'CATDESC',       { catdesc },       ...
	                    'DEPEND_0',      { depend_0 },      ...
	                    'DEPEND_1',      { depend_1 },      ...
	                    'DEPEND_2',      { depend_2 },      ...
	                    'DISPLAY_TYPE',  { display_type },  ...
	                    'FIELDNAM',      { fieldnam },      ...
	                    'FILLVAL',       { fillval },       ...
	                    'FORMAT',        { format },        ...
	                    'LABLAXIS',      { lablaxis },      ...
	                    'LABL_PTR_1',    { labl_ptr_1 },    ...
	                    'LABL_PTR_2',    { labl_ptr_2 },    ...
	                    'SI_CONVERSION', { si_conversion }, ...
	                    'TIME_BASE',     { time_base },     ...
	                    'UNITS',         { units },         ...
	                    'VALIDMIN',      { validmin },      ...
	                    'VALIDMAX',      { validmax },      ...
	                    'VAR_TYPE',      { var_type }       ...
	                  );

%------------------------------------------------------
% Write the File                                      |
%------------------------------------------------------
	spdfcdfwrite( afg_cal_file,                       ...
	              var_list,                           ...
	              'GlobalAttributes',   global_attrs, ...
	              'VariableAttributes', var_attrs,    ...
	              'VarDatatypes',       vardatatypes, ...
	              'VarCompress',        varcompress,  ...
	              'RecordBound',        recbound      ...
	            );
end