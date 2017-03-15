%
% Name
%   mms_fsm_l2plus_cal_dfg_write
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
function dfg_cal_file = mms_fsm_l2plus_cal_dfg_write(sc, mode, tstart, dfg_data, varargin)

	% Global variables
	%   - See mms_fsm_init.m
	global cal_path_root data_path_root dropbox_root hk_root log_path_root unh_data_root

%------------------------------------%
% Check Inputs                       %
%------------------------------------%

	% Defaults
	optdesc  = 'cal-dfg';
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
	assert( isa(dfg_data.t,     'int64'),  'dfg_data.t should be int64.' );
	assert( isa(dfg_data.f,     'single'), 'dfg_data.f should be single.');
	assert( isa(dfg_data.amp,   'single'), 'dfg_data.amp should be single.');
	assert( isa(dfg_data.phase, 'single'), 'dfg_data.phase should be single.');
	assert( isa(dfg_data.psd,   'single'), 'dfg_data.psd should be single.');
	assert( isa(dfg_data.flag,  'uint8'),  'dfg_data.flag should be uint8.');
	
	% convert the start time to yyyy-mm-dd
%	tstart = MrTimeParser(tstart, '%Y-%M-%dT%H:%m:%S', '%Y%M%d');

%------------------------------------%
% Version History                    %
%------------------------------------%
	
	% Describe the modifications to each version
	mods    = {  'v0.0.0 -- First version.'      ...
	             'v0.1.0 -- Use bigaussian fit.' ...
	          };

	parts = regexp( mods{end}, '^v([0-9]+)\.([0-9]+)\.([0-9]+)', 'tokens' );
	vx    = parts{1}{1};
	vy    = parts{1}{2};
	vz    = parts{1}{3};

%------------------------------------%
% Create Output File Name            %
%------------------------------------%
	% Output file
	version      = [vx '.' vy '.' vz];
	dfg_cal_file = mms_construct_filename(sc, instr, mode, level, ...
	                                      'TStart',  tstart,      ...
	                                      'Version', version,     ...
	                                      'OptDesc', optdesc);

	% Find the latest z-version
	%   - Look in both DROPBOX and DATA_PATH
	vz = mms_latest_zversion(dropbox_root, dfg_cal_file, 'RootDir', data_path_root);

	% Reform the file name
	version  = [vx '.' vy '.' sprintf('%0i', vz)];
	dfg_cal_file = mms_construct_filename(sc, instr, mode, level, ...
	                                      'TStart',  tstart,      ...
	                                      'Version', version,     ...
	                                      'OptDesc', optdesc);
	dfg_cal_file = fullfile(dropbox_root, dfg_cal_file);

	% Notify where file is located
	mrfprintf('logtext', ['Creating FSM dfg cal file at "' dfg_cal_file '".']);

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
	[~, logical_file_id, ext] = fileparts(dfg_cal_file);
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
	                       'TEXT',                       ['DFG calibration data used while created the FSM data product. ' ...
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
	f_vname         = [prefix 'f'                    suffix];
	c_index_vname   = [prefix 'component' '_index'   suffix];
	dfg_amp_vname   = [prefix 'amp'       '_omb_dfg' suffix];
	dfg_phase_vname = [prefix 'phase'     '_omb_dfg' suffix];
	dfg_psd_vname   = [prefix 'psd'       '_omb_dfg' suffix];
	flag_vname      = [prefix 'flag'                 suffix];
	flag_hist_vname = [prefix 'flag'      '_hist'    suffix];
	
	% FGM Histogram
	dfg_hamp_vname   = [prefix 'amp'   '_hist_dfg'  suffix];
	dfg_hphase_vname = [prefix 'phase' '_hist_dfg'  suffix];
	dfg_hpsd_vname   = [prefix 'psd'   '_hist_dfg'  suffix];
	dfg_famp_vname   = [prefix 'amp'   '_floor_dfg' suffix];
	dfg_fphase_vname = [prefix 'phase' '_floor_dfg' suffix];
	dfg_fpsd_vname   = [prefix 'psd'   '_floor_dfg' suffix];
	dfg_bamp_vname   = [prefix 'amp'   '_bins_dfg'  suffix];
	dfg_bphase_vname = [prefix 'phase' '_bins_dfg'  suffix];
	dfg_bpsd_vname   = [prefix 'psd'   '_bins_dfg'  suffix];
	
	% Metadata
	c_labl_vname      = [prefix 'component' '_labl'      suffix];
	bamp_labl_vname   = [prefix 'amp'       '_bins_labl' suffix];
	bphase_labl_vname = [prefix 'phase'     '_bins_labl' suffix];
	bpsd_labl_vname   = [prefix 'psd'       '_bins_labl' suffix];
	f_labl_vname      = [prefix 'f'         '_labl'      suffix];
	hflag_labl_vname  = [prefix 'flag'      '_hist_labl' suffix];

%------------------------------------------------------
% Metadata: Labels & Indices                          |
%------------------------------------------------------
	
	% Amplitude bins Label
	len       = length(dfg_data.amp_bins);
	powten    = floor( log10( len ) ) + 1;
	fmt       = sprintf('%%0%ii', powten);
	bamp_labl = strcat( {'amp'}, num2str( [1:len]', fmt) )';
	
	% Component label
	c_labl = {'X', 'Y', 'Z'};
	
	% Frequency Labels
	len    = length(dfg_data.f);
	powten = floor( log10( len ) ) + 1;
	fmt    = sprintf('%%0%ii', powten);
	f_labl = strcat( {'freq'}, num2str( [1:len]', fmt) )';
	
	% Hist Flag Label
	len        = length(dfg_data.hist_flag);
	powten     = floor( log10( len ) ) + 1;
	fmt        = sprintf('%%0%ii', powten);
	hflag_labl = strcat( {'histflag'}, num2str( [1:len]', fmt) )';
	
	% Phase bins Label
	len         = length(dfg_data.phase_bins);
	powten      = floor( log10( len ) ) + 1;
	fmt         = sprintf('%%0%ii', powten);
	bphase_labl = strcat( {'phase'}, num2str( [1:len]', fmt) )';
	
	% PSD bins Label
	len       = length(dfg_data.psd_bins);
	powten    = floor( log10( len ) ) + 1;
	fmt       = sprintf('%%0%ii', powten);
	bpsd_labl = strcat( {'psd'}, num2str( [1:len]', fmt) )';
	
	% Indices
	%   - CDF is 0-based
	c_index = uint8(0:2);

%------------------------------------------------------
% Variable Data                                       |
%------------------------------------------------------

	% Variable Data
	%   - Order as [ DEPEND_1,  DEPEND_2,  DEPEND_0 ]
	%   - Order as [ COMPONENT, FREQUENCY, TIME     ]
	var_list = { t_vname,           dfg_data.t,           ...
	             f_vname,           dfg_data.f,           ...
	             dfg_amp_vname,     permute( dfg_data.amp,   [3,1,2] ), ...
	             dfg_phase_vname,   permute( dfg_data.phase, [3,1,2] ), ...
	             dfg_psd_vname,     permute( dfg_data.psd,   [3,1,2] ), ...
	             flag_vname,        dfg_data.flag',       ...
	             bamp_labl_vname,   bamp_labl,            ...
	             bphase_labl_vname, bphase_labl,          ...
	             bpsd_labl_vname,   bpsd_labl,            ...
	             c_index_vname,     c_index,              ...
	             c_labl_vname,      c_labl,               ...
	             f_labl_vname,      f_labl,               ...
	             hflag_labl_vname,  hflag_labl,           ...
	             flag_hist_vname,   dfg_data.hist_flag,   ...
	             dfg_hamp_vname,    permute( dfg_data.amp_hist,   [3,2,1,4] ), ...
	             dfg_hphase_vname,  permute( dfg_data.phase_hist, [3,2,1,4] ), ...
	             dfg_hpsd_vname,    permute( dfg_data.psd_hist,   [3,2,1,4] ), ...
	             dfg_famp_vname,    permute( dfg_data.amp_floor,   [2,1,3] ), ...
	             dfg_fphase_vname,  permute( dfg_data.phase_floor, [2,1,3] ), ...
	             dfg_fpsd_vname,    permute( dfg_data.psd_floor,   [2,1,3] ), ...
	             dfg_bamp_vname,    dfg_data.amp_bins,    ...
	             dfg_bphase_vname,  dfg_data.phase_bins,  ...
	             dfg_bpsd_vname,    dfg_data.psd_bins     ...
	           };
	clear dfg_data

	% Record Variance
	recbound = { t_vname,          ...
	             dfg_amp_vname,    ...
	             dfg_phase_vname,  ...
	             dfg_psd_vname,    ...
	             flag_vname,       ...
	             dfg_hamp_vname,   ...
	             dfg_hphase_vname, ...
	             dfg_hpsd_vname,   ...
	             dfg_famp_vname,   ...
	             dfg_fphase_vname, ...
	             dfg_fpsd_vname    ...
	           };

	% Data types
	vardatatypes = { t_vname,           'cdf_time_tt2000', ...
	                 f_vname,           'cdf_float', ...
	                 dfg_amp_vname,     'cdf_float', ...
	                 dfg_phase_vname,   'cdf_float', ...
	                 dfg_psd_vname,     'cdf_float', ...
	                 flag_vname,        'cdf_uint1', ...
	                 bamp_labl_vname,   'cdf_char',  ...
	                 bphase_labl_vname, 'cdf_char',  ...
	                 bpsd_labl_vname,   'cdf_char',  ...
	                 c_index_vname,     'cdf_uint1', ...
	                 c_labl_vname,      'cdf_char',  ...
	                 f_labl_vname,      'cdf_char',  ...
	                 hflag_labl_vname,  'cdf_char',  ...
	                 flag_hist_vname,   'cdf_uint1', ...
	                 dfg_hamp_vname,    'cdf_uint4', ...
	                 dfg_hphase_vname,  'cdf_uint4', ...
	                 dfg_hpsd_vname,    'cdf_uint4', ...
	                 dfg_famp_vname,    'cdf_float', ...
	                 dfg_fphase_vname,  'cdf_float', ...
	                 dfg_fpsd_vname,    'cdf_float', ...
	                 dfg_bamp_vname,    'cdf_float', ...
	                 dfg_bphase_vname,  'cdf_float', ...
	                 dfg_bpsd_vname,    'cdf_float'  ...
	               };
	
	% CDF Compression Level
	varcompress = { dfg_amp_vname,     'gzip.6', ...
	                dfg_phase_vname,   'gzip.6', ...
	                dfg_psd_vname,     'gzip.6', ...
	                flag_vname,        'gzip.6', ...
	                flag_hist_vname,   'gzip.6', ...
	                dfg_hamp_vname,    'gzip.6', ...
	                dfg_hphase_vname,  'gzip.6', ...
	                dfg_hpsd_vname,    'gzip.6', ...
	                dfg_famp_vname,    'gzip.6', ...
	                dfg_fphase_vname,  'gzip.6', ...
	                dfg_fpsd_vname,    'gzip.6', ...
	                dfg_bamp_vname,    'gzip.6', ...
	                dfg_bphase_vname,  'gzip.6', ...
	                dfg_bpsd_vname,    'gzip.6'  ...
	              };
	
	singleton = { dfg_amp_vname,    ...
	              dfg_phase_vname,  ...
	              dfg_psd_vname,    ...
	              dfg_hamp_vname,   ...
	              dfg_hphase_vname, ...
	              dfg_hpsd_vname,   ...
	              dfg_famp_vname,   ...
	              dfg_fphase_vname, ...
	              dfg_fpsd_vname    ...
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
	labl_ptr_3    = {};
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
	
	% Component Index
	catdesc       = [ catdesc,       c_index_vname, 'Component index.' ];
	fieldnam      = [ fieldnam,      c_index_vname, 'Component index' ];
	fillval       = [ fillval,       c_index_vname, uint8(255) ];
	format        = [ format,        c_index_vname, 'I1' ];
	lablaxis      = [ lablaxis,      c_index_vname, 'Comp Index' ];
	validmin      = [ validmin,      c_index_vname, uint8(0) ];
	validmax      = [ validmax,      c_index_vname, uint8(2) ];
	var_type      = [ var_type,      c_index_vname, 'support_data' ];
	
	% AMP
	catdesc       = [ catdesc,       dfg_amp_vname, 'Wave amplitude of the FGM magnetic field.' ];
	depend_0      = [ depend_0,      dfg_amp_vname, t_vname ];
	depend_1      = [ depend_1,      dfg_amp_vname, f_vname ];
	depend_2      = [ depend_2,      dfg_amp_vname, c_index_vname ];
	display_type  = [ display_type,  dfg_amp_vname, 'spectrogram'];
	fieldnam      = [ fieldnam,      dfg_amp_vname, 'Amplitude' ];
	fillval       = [ fillval,       dfg_amp_vname, single(-1e31) ];
	format        = [ format,        dfg_amp_vname, 'f11.4' ];
	labl_ptr_1    = [ labl_ptr_1,    dfg_amp_vname, f_labl_vname ];
	labl_ptr_2    = [ labl_ptr_2,    dfg_amp_vname, c_labl_vname ];
	si_conversion = [ si_conversion, dfg_amp_vname, '1e-9>T' ];
	units         = [ units,         dfg_amp_vname, 'nT' ];
	validmin      = [ validmin,      dfg_amp_vname, single(0.0) ];
	validmax      = [ validmax,      dfg_amp_vname, single(1e5) ];
	var_type      = [ var_type,      dfg_amp_vname, 'data' ];
	
	% PHASE
	catdesc       = [ catdesc,       dfg_phase_vname, 'Wave phase of the FGM magnetic field.' ];
	depend_0      = [ depend_0,      dfg_phase_vname, t_vname ];
	depend_1      = [ depend_1,      dfg_phase_vname, f_vname ];
	depend_2      = [ depend_2,      dfg_phase_vname, c_index_vname ];
	display_type  = [ display_type,  dfg_phase_vname, 'spectrogram'];
	fieldnam      = [ fieldnam,      dfg_phase_vname, 'Phase' ];
	fillval       = [ fillval,       dfg_phase_vname, single(-1e31) ];
	format        = [ format,        dfg_phase_vname, 'f11.4' ];
	labl_ptr_1    = [ labl_ptr_1,    dfg_phase_vname, f_labl_vname ];
	labl_ptr_2    = [ labl_ptr_2,    dfg_phase_vname, c_labl_vname ];
	si_conversion = [ si_conversion, dfg_phase_vname, '1e-9>T' ];
	units         = [ units,         dfg_phase_vname, 'degrees' ];
	validmin      = [ validmin,      dfg_phase_vname, single(0.0) ];
	validmax      = [ validmax,      dfg_phase_vname, single(360.0) ];
	var_type      = [ var_type,      dfg_phase_vname, 'data' ];
	
	% PSD
	catdesc       = [ catdesc,       dfg_psd_vname, 'Power spectral density of the FGM magnetic field.' ];
	depend_0      = [ depend_0,      dfg_psd_vname, t_vname ];
	depend_1      = [ depend_1,      dfg_psd_vname, f_vname ];
	depend_2      = [ depend_2,      dfg_psd_vname, c_index_vname ];
	display_type  = [ display_type,  dfg_psd_vname, 'spectrogram'];
	fieldnam      = [ fieldnam,      dfg_psd_vname, 'Phase' ];
	fillval       = [ fillval,       dfg_psd_vname, single(-1e31) ];
	format        = [ format,        dfg_psd_vname, 'f11.4' ];
	labl_ptr_1    = [ labl_ptr_1,    dfg_psd_vname, f_labl_vname ];
	labl_ptr_2    = [ labl_ptr_2,    dfg_psd_vname, c_labl_vname ];
	si_conversion = [ si_conversion, dfg_psd_vname, '1e-9>T' ];
	units         = [ units,         dfg_psd_vname, 'nT^2/Hz' ];
	validmin      = [ validmin,      dfg_psd_vname, single(0.0) ];
	validmax      = [ validmax,      dfg_psd_vname, single(1e4) ];
	var_type      = [ var_type,      dfg_psd_vname, 'data' ];
	
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
	
%-------------------------------------------------------------------------
% FGM HISTOGRAM INFO /////////////////////////////////////////////////////
%-------------------------------------------------------------------------
	
	% HIST_AMP
	catdesc       = [ catdesc,       dfg_hamp_vname, 'FGM wave amplitude, Histogrammed in time at each frequency.' ];
	depend_0      = [ depend_0,      dfg_hamp_vname, f_vname ];
	depend_1      = [ depend_1,      dfg_hamp_vname, dfg_bamp_vname ];
	depend_2      = [ depend_2,      dfg_hamp_vname, flag_hist_vname ];
	depend_3      = [ depend_3,      dfg_hamp_vname, c_index_vname ];
	display_type  = [ display_type,  dfg_hamp_vname, 'spectrogram'];
	fieldnam      = [ fieldnam,      dfg_hamp_vname, 'Histogrammed wave amplitude (FGM)' ];
	fillval       = [ fillval,       dfg_hamp_vname, uint32(4294967295) ];
	format        = [ format,        dfg_hamp_vname, 'I10' ];
	labl_ptr_1    = [ labl_ptr_1,    dfg_hamp_vname, bamp_labl_vname ];
	labl_ptr_2    = [ labl_ptr_2,    dfg_hamp_vname, hflag_labl_vname ];
	labl_ptr_3    = [ labl_ptr_3,    dfg_hamp_vname, c_labl_vname ];
	si_conversion = [ si_conversion, dfg_hamp_vname, '1e-9>T' ];
	units         = [ units,         dfg_hamp_vname, 'degrees' ];
	validmin      = [ validmin,      dfg_hamp_vname, uint32(0) ];
	validmax      = [ validmax,      dfg_hamp_vname, uint32(4294967294) ];
	var_type      = [ var_type,      dfg_hamp_vname, 'data' ];
	
	% HIST_PHASE
	catdesc       = [ catdesc,       dfg_hphase_vname, 'FGM wave phase, Histogrammed in time at each frequency.' ];
	depend_0      = [ depend_0,      dfg_hphase_vname, f_vname ];
	depend_1      = [ depend_1,      dfg_hphase_vname, dfg_bphase_vname ];
	depend_2      = [ depend_2,      dfg_hphase_vname, flag_hist_vname ];
	depend_3      = [ depend_3,      dfg_hphase_vname, c_index_vname ];
	display_type  = [ display_type,  dfg_hphase_vname, 'spectrogram'];
	fieldnam      = [ fieldnam,      dfg_hphase_vname, 'Histogrammed wave phase (FGM)' ];
	fillval       = [ fillval,       dfg_hphase_vname, uint32(4294967295) ];
	format        = [ format,        dfg_hphase_vname, 'I10' ];
	labl_ptr_1    = [ labl_ptr_1,    dfg_hphase_vname, bphase_labl_vname ];
	labl_ptr_2    = [ labl_ptr_2,    dfg_hphase_vname, hflag_labl_vname ];
	labl_ptr_3    = [ labl_ptr_3,    dfg_hphase_vname, c_labl_vname ];
	si_conversion = [ si_conversion, dfg_hphase_vname, '1e-9>T' ];
	units         = [ units,         dfg_hphase_vname, 'degrees' ];
	validmin      = [ validmin,      dfg_hphase_vname, uint32(0) ];
	validmax      = [ validmax,      dfg_hphase_vname, uint32(4294967294) ];
	var_type      = [ var_type,      dfg_hphase_vname, 'data' ];
	
	% HIST_PSD
	catdesc       = [ catdesc,       dfg_hpsd_vname, 'FGM wave power, Histogrammed in time at each frequency.' ];
	depend_0      = [ depend_0,      dfg_hpsd_vname, f_vname ];
	depend_1      = [ depend_1,      dfg_hpsd_vname, dfg_bpsd_vname ];
	depend_2      = [ depend_2,      dfg_hpsd_vname, flag_hist_vname ];
	depend_3      = [ depend_3,      dfg_hpsd_vname, c_index_vname ];
	display_type  = [ display_type,  dfg_hpsd_vname, 'spectrogram'];
	fieldnam      = [ fieldnam,      dfg_hpsd_vname, 'Histogrammed wave power (FGM)' ];
	fillval       = [ fillval,       dfg_hpsd_vname, uint32(4294967295) ];
	format        = [ format,        dfg_hpsd_vname, 'I10' ];
	labl_ptr_1    = [ labl_ptr_1,    dfg_hpsd_vname, bpsd_labl_vname ];
	labl_ptr_2    = [ labl_ptr_2,    dfg_hpsd_vname, hflag_labl_vname ];
	labl_ptr_3    = [ labl_ptr_3,    dfg_hpsd_vname, c_labl_vname ];
	si_conversion = [ si_conversion, dfg_hpsd_vname, '1e-9>T' ];
	units         = [ units,         dfg_hpsd_vname, 'degrees' ];
	validmin      = [ validmin,      dfg_hpsd_vname, uint32(0) ];
	validmax      = [ validmax,      dfg_hpsd_vname, uint32(4294967294) ];
	var_type      = [ var_type,      dfg_hpsd_vname, 'data' ];
	
	% HIST_FLAG
	catdesc       = [ catdesc,       flag_hist_vname, 'Bitflag indicating regime. Bit=unset/set: 1=slow/fast(brst), 2=hi/lo-range, 3=Deck32/64' ];
	fieldnam      = [ fieldnam,      flag_hist_vname, 'Histogram bitflag' ];
	fillval       = [ fillval,       flag_hist_vname, uint32(255) ];
	format        = [ format,        flag_hist_vname, 'I3' ];
	lablaxis      = [ lablaxis,      flag_hist_vname, 'Hist flag' ];
	validmin      = [ validmin,      flag_hist_vname, uint8(0) ];
	validmax      = [ validmax,      flag_hist_vname, uint8(254) ];
	var_type      = [ var_type,      flag_hist_vname, 'support_data' ];
	
%-------------------------------------------------------------------------
% FGM NOISE FLOOR INFO ///////////////////////////////////////////////////
%-------------------------------------------------------------------------
	
	% FLOOR_AMP
	catdesc       = [ catdesc,       dfg_famp_vname, ['Log base 10 of the most common wave amplitude, obtained by histogramming ' ...
	                                                  'wave amplitude at each frequency over the course of a day.'] ];
	depend_0      = [ depend_0,      dfg_famp_vname, f_vname ];
	depend_1      = [ depend_1,      dfg_famp_vname, flag_hist_vname ];
	depend_2      = [ depend_2,      dfg_famp_vname, c_index_vname ];
	display_type  = [ display_type,  dfg_famp_vname, 'time_series'];
	fieldnam      = [ fieldnam,      dfg_famp_vname, 'Wave amplitude' ];
	fillval       = [ fillval,       dfg_famp_vname, single(-1e31) ];
	format        = [ format,        dfg_famp_vname, 'f11.4' ];
	labl_ptr_1    = [ labl_ptr_1,    dfg_famp_vname, hflag_labl_vname ];
	labl_ptr_2    = [ labl_ptr_2,    dfg_famp_vname, c_labl_vname ];
	si_conversion = [ si_conversion, dfg_famp_vname, '1e-9>T' ];
	units         = [ units,         dfg_famp_vname, 'degrees' ];
	validmin      = [ validmin,      dfg_famp_vname, single(-10.0) ];
	validmax      = [ validmax,      dfg_famp_vname, single(10.0) ];
	var_type      = [ var_type,      dfg_famp_vname, 'data' ];
	
	% FLOOR_PHASE
	catdesc       = [ catdesc,       dfg_fphase_vname, ['The most common wave phase, obtained by histogramming ' ...
	                                                    'wave phase at each frequency over the course of a day.'] ];
	depend_0      = [ depend_0,      dfg_fphase_vname, f_vname ];
	depend_1      = [ depend_1,      dfg_fphase_vname, flag_hist_vname ];
	depend_2      = [ depend_2,      dfg_fphase_vname, c_index_vname ];
	display_type  = [ display_type,  dfg_fphase_vname, 'time_series'];
	fieldnam      = [ fieldnam,      dfg_fphase_vname, 'Wave phase' ];
	fillval       = [ fillval,       dfg_fphase_vname, single(-1e31) ];
	format        = [ format,        dfg_fphase_vname, 'f11.4' ];
	labl_ptr_1    = [ labl_ptr_1,    dfg_fphase_vname, hflag_labl_vname ];
	labl_ptr_2    = [ labl_ptr_2,    dfg_fphase_vname, c_labl_vname ];
	si_conversion = [ si_conversion, dfg_fphase_vname, '1e0>degrees' ];
	units         = [ units,         dfg_fphase_vname, 'degrees' ];
	validmin      = [ validmin,      dfg_fphase_vname, single(-180.0) ];
	validmax      = [ validmax,      dfg_fphase_vname, single(180.0) ];
	var_type      = [ var_type,      dfg_fphase_vname, 'data' ];
	
	% FLOOR_PSD
	catdesc       = [ catdesc,       dfg_fpsd_vname, ['Log base 10 of the most common wave power, obtained by histogramming ' ...
	                                                  'wave power at each frequency over the course of a day.'] ];
	depend_0      = [ depend_0,      dfg_fpsd_vname, f_vname ];
	depend_1      = [ depend_1,      dfg_fpsd_vname, flag_hist_vname ];
	depend_2      = [ depend_2,      dfg_fpsd_vname, c_index_vname ];
	display_type  = [ display_type,  dfg_fpsd_vname, 'time_series'];
	fieldnam      = [ fieldnam,      dfg_fpsd_vname, 'Power spectral density' ];
	fillval       = [ fillval,       dfg_fpsd_vname, single(-1e31) ];
	format        = [ format,        dfg_fpsd_vname, 'f11.4' ];
	labl_ptr_1    = [ labl_ptr_1,    dfg_fpsd_vname, hflag_labl_vname ];
	labl_ptr_2    = [ labl_ptr_2,    dfg_fpsd_vname, c_labl_vname ];
	si_conversion = [ si_conversion, dfg_fpsd_vname, '1e-18>T^2 Hz^-1' ];
	units         = [ units,         dfg_fpsd_vname, 'nT^2 Hz^-1' ];
	validmin      = [ validmin,      dfg_fpsd_vname, single(-15.0) ];
	validmax      = [ validmax,      dfg_fpsd_vname, single(15.0) ];
	var_type      = [ var_type,      dfg_fpsd_vname, 'data' ];
	
%-------------------------------------------------------------------------
% FGM BIN INFO ///////////////////////////////////////////////////////////
%-------------------------------------------------------------------------
	
	% BINS_AMP
	catdesc       = [ catdesc,       dfg_bamp_vname, 'Histogram bins for wave amplitude.' ];
	fieldnam      = [ fieldnam,      dfg_bamp_vname, 'Wave amplitude bins' ];
	fillval       = [ fillval,       dfg_bamp_vname, single(-1e31) ];
	format        = [ format,        dfg_bamp_vname, 'f11.4' ];
	lablaxis      = [ lablaxis,      dfg_bamp_vname, 'Amp Bins' ];
	si_conversion = [ si_conversion, dfg_bamp_vname, '1e-9>T' ];
	units         = [ units,         dfg_bamp_vname, 'nT' ];
	validmin      = [ validmin,      dfg_bamp_vname, single(-10.0) ];
	validmax      = [ validmax,      dfg_bamp_vname, single(10.0) ];
	var_type      = [ var_type,      dfg_bamp_vname, 'support_data' ];
	
	% BINS_PHASE
	catdesc       = [ catdesc,       dfg_bphase_vname, 'Histogram bins for wave phase.' ];
	fieldnam      = [ fieldnam,      dfg_bphase_vname, 'Wave amplitude bins' ];
	fillval       = [ fillval,       dfg_bphase_vname, single(-1e31) ];
	format        = [ format,        dfg_bphase_vname, 'f11.4' ];
	lablaxis      = [ lablaxis,      dfg_bphase_vname, 'Phase Bins' ];
	si_conversion = [ si_conversion, dfg_bphase_vname, '1e0>degrees' ];
	units         = [ units,         dfg_bphase_vname, 'degrees' ];
	validmin      = [ validmin,      dfg_bphase_vname, single(-180.0) ];
	validmax      = [ validmax,      dfg_bphase_vname, single(180.0) ];
	var_type      = [ var_type,      dfg_bphase_vname, 'support_data' ];
	
	% BINS_PSD
	catdesc       = [ catdesc,       dfg_bpsd_vname, 'Histogram bins for wave amplitude.' ];
	fieldnam      = [ fieldnam,      dfg_bpsd_vname, 'Wave amplitude bins' ];
	fillval       = [ fillval,       dfg_bpsd_vname, single(-1e31) ];
	format        = [ format,        dfg_bpsd_vname, 'f11.4' ];
	lablaxis      = [ lablaxis,      dfg_bpsd_vname, 'PSD Bins' ];
	si_conversion = [ si_conversion, dfg_bpsd_vname, '1e-18>T^2 Hz^-1' ];
	units         = [ units,         dfg_bpsd_vname, 'nT^2 Hz^-1' ];
	validmin      = [ validmin,      dfg_bpsd_vname, single(-10.0) ];
	validmax      = [ validmax,      dfg_bpsd_vname, single(10.0) ];
	var_type      = [ var_type,      dfg_bpsd_vname, 'support_data' ];
	
%-------------------------------------------------------------------------
% FGM METADATA INFO //////////////////////////////////////////////////////
%-------------------------------------------------------------------------
	
	% AMPLITUDE BIN LABEL
	catdesc       = [ catdesc,       bamp_labl_vname, 'Axis labels for amplitude bins' ];
	fieldnam      = [ fieldnam,      bamp_labl_vname, 'Amplitude labels' ];
	format        = [ format,        bamp_labl_vname, 'A7' ];
	var_type      = [ var_type,      bamp_labl_vname, 'metadata' ];
	
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
	
	% HISTOGRAM FLAG LABEL
	catdesc       = [ catdesc,       hflag_labl_vname, 'Axis label for histogram flag.' ];
	fieldnam      = [ fieldnam,      hflag_labl_vname, 'Hist flag labels' ];
	format        = [ format,        hflag_labl_vname, 'A13' ];
	var_type      = [ var_type,      hflag_labl_vname, 'metadata' ];
	
	% PHASE BINS LABEL
	catdesc       = [ catdesc,       bphase_labl_vname, 'Axis label for phase bins.' ];
	fieldnam      = [ fieldnam,      bphase_labl_vname, 'Phase labels' ];
	format        = [ format,        bphase_labl_vname, 'A13' ];
	var_type      = [ var_type,      bphase_labl_vname, 'metadata' ];
	
	% PSD BINS LABEL
	catdesc       = [ catdesc,       bpsd_labl_vname, 'Axis label for psd bins.' ];
	fieldnam      = [ fieldnam,      bpsd_labl_vname, 'PSD labels' ];
	format        = [ format,        bpsd_labl_vname, 'A13' ];
	var_type      = [ var_type,      bpsd_labl_vname, 'metadata' ];
	
%-------------------------------------------------------------------------
% Collect Variable Attributes ////////////////////////////////////////////
%-------------------------------------------------------------------------
	% Collect into a structure
	var_attrs = struct( 'CATDESC',       { catdesc },       ...
	                    'DEPEND_0',      { depend_0 },      ...
	                    'DEPEND_1',      { depend_1 },      ...
	                    'DEPEND_2',      { depend_2 },      ...
	                    'DEPEND_3',      { depend_3 },      ...
	                    'DISPLAY_TYPE',  { display_type },  ...
	                    'FIELDNAM',      { fieldnam },      ...
	                    'FILLVAL',       { fillval },       ...
	                    'FORMAT',        { format },        ...
	                    'LABLAXIS',      { lablaxis },      ...
	                    'LABL_PTR_1',    { labl_ptr_1 },    ...
	                    'LABL_PTR_2',    { labl_ptr_2 },    ...
	                    'LABL_PTR_3',    { labl_ptr_3 },    ...
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
	spdfcdfwrite( dfg_cal_file,                       ...
	              var_list,                           ...
	              'GlobalAttributes',   global_attrs, ...
	              'VariableAttributes', var_attrs,    ...
	              'VarDatatypes',       vardatatypes, ...
	              'VarCompress',        varcompress,  ...
	              'RecordBound',        recbound,     ...
	              'Singleton',          singleton     ...
	            );
end