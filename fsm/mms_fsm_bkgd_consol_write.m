%
% Name
%   mms_fsm_l2plus_cal_write
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
function fsm_cal_file = mms_fsm_bkgd_consol_write(sc, mode, optdesc, tstart, data, varargin)

	% Global variables
	%   - See mms_fsm_init.m
	global cal_path_root data_path_root dropbox_root hk_root log_path_root unh_data_root

%------------------------------------%
% Check Inputs                       %
%------------------------------------%

	% Defaults
	tf_empty = false;
	parents  = {};

	% Check inputs
	nOptArgs = length(varargin);
	for ii = 1 : 2 : nOptArgs
		switch varargin{ii}
			case 'EmptyFile'
				tf_empty = varargin{ii+1};
			case 'Parents'
				parents = varargin{ii+1};
			otherwise
				error(['Optional argument not recognized: "' varargin{ii} '".']);
		end
	end
	
	% Constants
	instr = 'fsm';
	level = 'l2plus';

%------------------------------------%
% Verify Data                        %
%------------------------------------%

	%
	% Check sizes
	%
	assert( isa( data.f,          'single'), 'data.f should be single.');
	assert( isa( data.flag,       'uint8'),  'data.flag should be uint8.');
	assert( isa( data.comp,       'uint8'),  'data.comp should be uint8.');
	assert( isa( data.amp_bins,   'single'), 'data.amp_bins should be single.');
	assert( isa( data.phase_bins, 'single'), 'data.phase_bins should be single.');
	assert( isa( data.psd_bins,   'single'), 'data.psd_bins should be single.');
	assert( isa( data.amp_hist,   'uint32'), 'data.amp_hist should be uint32.');
	assert( isa( data.phase_hist, 'uint32'), 'data.phase_hist should be uint32.');
	assert( isa( data.psd_hist,   'uint32'), 'data.psd_hist should be uint32.');
	assert( isa( data.psd_floor,  'single'), 'data.psd_floor should be single.');

%------------------------------------%
% Version History                    %
%------------------------------------%
	
	% Describe the modifications to each version
	mods    = {  'v0.0.0 -- First version.' ...
	          };

	parts = regexp( mods{end}, '^v([0-9]+)\.([0-9]+)\.([0-9]+)', 'tokens' );
	vx    = parts{1}{1};
	vy    = parts{1}{2};
	vz    = parts{1}{3};

%------------------------------------%
% Create Output File Name            %
%------------------------------------%
	% Output file
	version = [vx '.' vy '.' vz];
	fsm_cal_file = mms_construct_filename(sc, instr, mode, level, ...
	                                      'TStart',  tstart,      ...
	                                      'Version', version,     ...
	                                      'OptDesc', optdesc);

	% Find the latest z-version
	%   - Look in both DROPBOX and DATA_PATH
	vz = mms_latest_zversion(dropbox_root, fsm_cal_file, 'RootDir', data_path_root);

	% Reform the file name
	version  = [vx '.' vy '.' sprintf('%0i', vz)];
	fsm_cal_file = mms_construct_filename(sc, instr, mode, level, ...
	                                      'TStart',  tstart,      ...
	                                      'Version', version,     ...
	                                      'OptDesc', optdesc);
	fsm_cal_file = fullfile(dropbox_root, fsm_cal_file);

	% Notify where file is located
	mrfprintf('logtext', ['Creating FSM cal file at "' fsm_cal_file '".']);

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
	[~, logical_file_id, ext] = fileparts(fsm_cal_file);
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
	                       'Logical_source_description', 'Level 2-Plus Fluxgate-Searchcoil Merged Magnetic Field Calibration Data', ...
	                       'Mission_group',              'MMS', ...
	                       'PI_affiliation',             'UNH, LPP, IWF, UCLA, UCLA, IWF, UNH', ...
	                       'PI_name',                    ['R.B. Torbert, O. LeContel, W. Magnes, C.T. Russell, R. Strangeway, ', ...
	                                                      'D. Fischer, M.R. Argall.'], ...
	                       'Project',                    'STP>Solar Terrestrial Physics', ...
	                       'Source_name',                source_name, ...
	                       'TEXT',                       ['Magnetometer calibration data used while created the FSM data product. ' ...
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
	t_vname          = 'Epoch';
	f_vname          = [prefix 'f'                     suffix];
	c_index_vname    = [prefix 'component'  '_index'   suffix];
	flag_vname       = [prefix 'flag'                  suffix];
	
	% Histogram
	amp_hist_vname    = [prefix 'amp'   '_hist' suffix];
	phase_hist_vname  = [prefix 'phase' '_hist' suffix];
	psd_hist_vname    = [prefix 'psd'   '_hist' suffix];
	
	% Histogram bins
	amp_bins_vname   = [prefix 'amp'   '_bins'  suffix];
	phase_bins_vname = [prefix 'phase' '_bins'  suffix];
	psd_bins_vname   = [prefix 'psd'   '_bins'  suffix];
	
	% Noise floor
	amp_floor_vname   = [prefix 'amp'   '_floor' suffix];
	phase_floor_vname = [prefix 'phase' '_floor' suffix];
	psd_floor_vname   = [prefix 'psd'   '_floor' suffix];
	
	% Metadata
	f_labl_vname         = [prefix 'f'         '_labl'      suffix];
	c_labl_vname         = [prefix 'component' '_labl'      suffix];
	flag_labl_vname      = [prefix 'flag'      '_hist_labl' suffix];
	amp_bin_labl_vname   = [prefix 'amp'       '_bins_labl' suffix];
	phase_bin_labl_vname = [prefix 'phase'     '_bins_labl' suffix];
	psd_bin_labl_vname   = [prefix 'psd'       '_bins_labl' suffix];

%------------------------------------------------------
% Metadata: Labels & Indices                          |
%------------------------------------------------------
	
	% Component label & indices
	%   - CDF is 0-based
	c_labl  = {'X', 'Y', 'Z'};
	c_index = uint8(0:2);
	
	% Frequency Labels
	len    = length(data.f);
	powten = floor( log10( len ) ) + 1;
	fmt    = sprintf('%%0%ii', powten);
	f_labl = strcat( {'freq'}, num2str( [1:len]', fmt) )';
	
	% Amplitude bins Label
	len          = length(data.amp_bins);
	powten       = floor( log10( len ) ) + 1;
	fmt          = sprintf('%%0%ii', powten);
	amp_bin_labl = strcat( {'amp'}, num2str( [1:len]', fmt) )';
	
	% Hist Flag Label
	len       = length(data.flag);
	powten    = floor( log10( len ) ) + 1;
	fmt       = sprintf('%%0%ii', powten);
	flag_labl = strcat( {'flag'}, num2str( [1:len]', fmt) )';
	
	% Phase bins Label
	len            = length(data.phase_bins);
	powten         = floor( log10( len ) ) + 1;
	fmt            = sprintf('%%0%ii', powten);
	phase_bin_labl = strcat( {'phase'}, num2str( [1:len]', fmt) )';
	
	% PSD bins Label
	len          = length(data.psd_bins);
	powten       = floor( log10( len ) ) + 1;
	fmt          = sprintf('%%0%ii', powten);
	psd_bin_labl = strcat( {'psd'}, num2str( [1:len]', fmt) )';

%------------------------------------------------------
% Variable Data                                       |
%------------------------------------------------------

	% Variable Data
	%   - Order as [ DEPEND_1,  DEPEND_2,  DEPEND_0 ]
	%   - Order as [ COMPONENT, FREQUENCY, TIME     ]
	var_list = { t_vname,              [],              ...
	             f_vname,              data.f,          ...
	             c_index_vname,        c_index,         ...
	             flag_vname,           data.flag,       ...
	             c_labl_vname,         c_labl,          ...
	             f_labl_vname,         f_labl,          ...
	             flag_labl_vname,      flag_labl,       ...
	             amp_bin_labl_vname,   amp_bin_labl,    ...
	             phase_bin_labl_vname, phase_bin_labl,  ...
	             psd_bin_labl_vname,   psd_bin_labl,    ...
	             amp_bins_vname,       data.amp_bins,   ...
	             phase_bins_vname,     data.phase_bins, ...
	             psd_bins_vname,       data.psd_bins    ...
	             amp_hist_vname,       permute( data.amp_hist,   [3,2,1,4] ), ...
	             phase_hist_vname,     permute( data.phase_hist, [3,2,1,4] ), ...
	             psd_hist_vname,       permute( data.psd_hist,   [3,2,1,4] ), ...
	             amp_floor_vname,      [], ... % permute( data.amp_floor,   [2,1,3] ),  ...
	             phase_floor_vname,    [], ... % permute( data.phase_floor, [2,1,3] ),  ...
	             psd_floor_vname,      permute( data.psd_floor,   [2,1,3] )   ...
	           };
	clear data

	% Record Variance
	recbound = { t_vname,           ...
	             amp_hist_vname,    ...
	             phase_hist_vname,  ...
	             psd_hist_vname,    ...
	             amp_floor_vname,   ...
	             phase_floor_vname, ...
	             psd_floor_vname    ...
	           };

	% Data types
	vardatatypes = { t_vname,              'cdf_time_tt2000', ...
	                 f_vname,              'cdf_float', ...
	                 c_index_vname,        'cdf_uint1', ...
	                 flag_vname,           'cdf_uint1', ...
	                 amp_bin_labl_vname,   'cdf_char',  ...
	                 phase_bin_labl_vname, 'cdf_char',  ...
	                 psd_bin_labl_vname,   'cdf_char',  ...
	                 c_labl_vname,         'cdf_char',  ...
	                 f_labl_vname,         'cdf_char',  ...
	                 flag_vname,           'cdf_uint1', ...
	                 amp_bins_vname,       'cdf_float', ...
	                 phase_bins_vname,     'cdf_float', ...
	                 psd_bins_vname        'cdf_float', ...
	                 amp_hist_vname,       'cdf_uint4', ...
	                 phase_hist_vname,     'cdf_uint4', ...
	                 psd_hist_vname,       'cdf_uint4', ...
	                 amp_floor_vname,      'cdf_float', ...
	                 phase_floor_vname,    'cdf_float', ...
	                 psd_floor_vname,      'cdf_float'  ...
	               };
	
	% CDF Compression Level
	varcompress = { flag_vname,        'gzip.6', ...
	                amp_hist_vname,    'gzip.6', ...
	                phase_hist_vname,  'gzip.6', ...
	                psd_hist_vname,    'gzip.6', ...
	                amp_floor_vname,   'gzip.6', ...
	                phase_floor_vname, 'gzip.6', ...
	                psd_floor_vname,   'gzip.6', ...
	                amp_bins_vname,    'gzip.6', ...
	                phase_bins_vname,  'gzip.6', ...
	                psd_bins_vname     'gzip.6'  ...
	              };
	
	singleton = { amp_hist_vname,    ...
	              phase_hist_vname,  ...
	              psd_hist_vname,    ...
	              amp_floor_vname,   ...
	              phase_floor_vname, ...
	              psd_floor_vname    ...
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
% HISTOGRAM INFO /////////////////////////////////////////////////////////
%-------------------------------------------------------------------------
	
	% AMPLITUDE
	catdesc       = [ catdesc,       amp_hist_vname, 'Wave amplitude, Histogrammed in time at each frequency.' ];
	depend_0      = [ depend_0,      amp_hist_vname, f_vname ];
	depend_1      = [ depend_1,      amp_hist_vname, amp_bins_vname ];
	depend_2      = [ depend_2,      amp_hist_vname, flag_vname ];
	depend_3      = [ depend_3,      amp_hist_vname, c_index_vname ];
	display_type  = [ display_type,  amp_hist_vname, 'spectrogram'];
	fieldnam      = [ fieldnam,      amp_hist_vname, 'Histogrammed wave amplitude' ];
	fillval       = [ fillval,       amp_hist_vname, uint32(4294967295) ];
	format        = [ format,        amp_hist_vname, 'I10' ];
	labl_ptr_1    = [ labl_ptr_1,    amp_hist_vname, amp_bin_labl_vname ];
	labl_ptr_2    = [ labl_ptr_2,    amp_hist_vname, flag_labl_vname ];
	labl_ptr_3    = [ labl_ptr_3,    amp_hist_vname, c_labl_vname ];
	si_conversion = [ si_conversion, amp_hist_vname, '1e-9>T' ];
	units         = [ units,         amp_hist_vname, 'degrees' ];
	validmin      = [ validmin,      amp_hist_vname, uint32(0) ];
	validmax      = [ validmax,      amp_hist_vname, uint32(4294967294) ];
	var_type      = [ var_type,      amp_hist_vname, 'data' ];
	
	% PHASE
	catdesc       = [ catdesc,       phase_hist_vname, 'Wave phase, Histogrammed in time at each frequency.' ];
	depend_0      = [ depend_0,      phase_hist_vname, f_vname ];
	depend_1      = [ depend_1,      phase_hist_vname, phase_bins_vname ];
	depend_2      = [ depend_2,      phase_hist_vname, flag_vname ];
	depend_3      = [ depend_3,      phase_hist_vname, c_index_vname ];
	display_type  = [ display_type,  phase_hist_vname, 'spectrogram'];
	fieldnam      = [ fieldnam,      phase_hist_vname, 'Histogrammed wave phase' ];
	fillval       = [ fillval,       phase_hist_vname, uint32(4294967295) ];
	format        = [ format,        phase_hist_vname, 'I10' ];
	labl_ptr_1    = [ labl_ptr_1,    phase_hist_vname, phase_bin_labl_vname ];
	labl_ptr_2    = [ labl_ptr_2,    phase_hist_vname, flag_labl_vname ];
	labl_ptr_3    = [ labl_ptr_3,    phase_hist_vname, c_labl_vname ];
	si_conversion = [ si_conversion, phase_hist_vname, '1e-9>T' ];
	units         = [ units,         phase_hist_vname, 'degrees' ];
	validmin      = [ validmin,      phase_hist_vname, uint32(0) ];
	validmax      = [ validmax,      phase_hist_vname, uint32(4294967294) ];
	var_type      = [ var_type,      phase_hist_vname, 'data' ];
	
	% POWER SPECTRAL DENSITY
	catdesc       = [ catdesc,       psd_hist_vname, 'Wave power, Histogrammed in time at each frequency.' ];
	depend_0      = [ depend_0,      psd_hist_vname, f_vname ];
	depend_1      = [ depend_1,      psd_hist_vname, psd_bins_vname ];
	depend_2      = [ depend_2,      psd_hist_vname, flag_vname ];
	depend_3      = [ depend_3,      psd_hist_vname, c_index_vname ];
	display_type  = [ display_type,  psd_hist_vname, 'spectrogram'];
	fieldnam      = [ fieldnam,      psd_hist_vname, 'Histogrammed wave power' ];
	fillval       = [ fillval,       psd_hist_vname, uint32(4294967295) ];
	format        = [ format,        psd_hist_vname, 'I10' ];
	labl_ptr_1    = [ labl_ptr_1,    psd_hist_vname, psd_bin_labl_vname ];
	labl_ptr_2    = [ labl_ptr_2,    psd_hist_vname, flag_labl_vname ];
	labl_ptr_3    = [ labl_ptr_3,    psd_hist_vname, c_labl_vname ];
	si_conversion = [ si_conversion, psd_hist_vname, '1e-9>T' ];
	units         = [ units,         psd_hist_vname, 'degrees' ];
	validmin      = [ validmin,      psd_hist_vname, uint32(0) ];
	validmax      = [ validmax,      psd_hist_vname, uint32(4294967294) ];
	var_type      = [ var_type,      psd_hist_vname, 'data' ];
	
%-------------------------------------------------------------------------
% NOISE FLOOR INFO ///////////////////////////////////////////////////////
%-------------------------------------------------------------------------
	
	% FLOOR_AMP
	catdesc       = [ catdesc,       amp_floor_vname, ['Log base 10 of the most common wave amplitude, obtained by histogramming ' ...
	                                                  'wave amplitude at each frequency over the course of a day.'] ];
	depend_0      = [ depend_0,      amp_floor_vname, f_vname ];
	depend_1      = [ depend_1,      amp_floor_vname, flag_vname ];
	depend_2      = [ depend_2,      amp_floor_vname, c_index_vname ];
	display_type  = [ display_type,  amp_floor_vname, 'time_series'];
	fieldnam      = [ fieldnam,      amp_floor_vname, 'Wave amplitude' ];
	fillval       = [ fillval,       amp_floor_vname, single(-1e31) ];
	format        = [ format,        amp_floor_vname, 'f11.4' ];
	labl_ptr_1    = [ labl_ptr_1,    amp_floor_vname, flag_labl_vname ];
	labl_ptr_2    = [ labl_ptr_2,    amp_floor_vname, c_labl_vname ];
	si_conversion = [ si_conversion, amp_floor_vname, '1e-9>T' ];
	units         = [ units,         amp_floor_vname, 'degrees' ];
	validmin      = [ validmin,      amp_floor_vname, single(-10.0) ];
	validmax      = [ validmax,      amp_floor_vname, single(10.0) ];
	var_type      = [ var_type,      amp_floor_vname, 'data' ];
	
	% FLOOR_PHASE
	catdesc       = [ catdesc,       phase_floor_vname, ['The most common wave phase, obtained by histogramming ' ...
	                                                    'wave phase at each frequency over the course of a day.'] ];
	depend_0      = [ depend_0,      phase_floor_vname, f_vname ];
	depend_1      = [ depend_1,      phase_floor_vname, flag_vname ];
	depend_2      = [ depend_2,      phase_floor_vname, c_index_vname ];
	display_type  = [ display_type,  phase_floor_vname, 'time_series'];
	fieldnam      = [ fieldnam,      phase_floor_vname, 'Wave phase' ];
	fillval       = [ fillval,       phase_floor_vname, single(-1e31) ];
	format        = [ format,        phase_floor_vname, 'f11.4' ];
	labl_ptr_1    = [ labl_ptr_1,    phase_floor_vname, flag_labl_vname ];
	labl_ptr_2    = [ labl_ptr_2,    phase_floor_vname, c_labl_vname ];
	si_conversion = [ si_conversion, phase_floor_vname, '1e0>degrees' ];
	units         = [ units,         phase_floor_vname, 'degrees' ];
	validmin      = [ validmin,      phase_floor_vname, single(-180.0) ];
	validmax      = [ validmax,      phase_floor_vname, single(180.0) ];
	var_type      = [ var_type,      phase_floor_vname, 'data' ];
	
	
	% FLOOR_PSD
	catdesc       = [ catdesc,       psd_floor_vname, ['Log base 10 of the most common wave power, obtained by histogramming ' ...
	                                                  'wave power at each frequency over the course of a day.'] ];
	depend_0      = [ depend_0,      psd_floor_vname, f_vname ];
	depend_1      = [ depend_1,      psd_floor_vname, flag_vname ];
	depend_2      = [ depend_2,      psd_floor_vname, c_index_vname ];
	display_type  = [ display_type,  psd_floor_vname, 'time_series'];
	fieldnam      = [ fieldnam,      psd_floor_vname, 'Power spectral density' ];
	fillval       = [ fillval,       psd_floor_vname, single(-1e31) ];
	format        = [ format,        psd_floor_vname, 'f11.4' ];
	labl_ptr_1    = [ labl_ptr_1,    psd_floor_vname, flag_labl_vname ];
	labl_ptr_2    = [ labl_ptr_2,    psd_floor_vname, c_labl_vname ];
	si_conversion = [ si_conversion, psd_floor_vname, '1e-18>T^2 Hz^-1' ];
	units         = [ units,         psd_floor_vname, 'nT^2 Hz^-1' ];
	validmin      = [ validmin,      psd_floor_vname, single(-15.0) ];
	validmax      = [ validmax,      psd_floor_vname, single(15.0) ];
	var_type      = [ var_type,      psd_floor_vname, 'data' ];
	
%-------------------------------------------------------------------------
% HISTOGRAM BINS INFO ////////////////////////////////////////////////////
%-------------------------------------------------------------------------
	
	% BINS_AMP
	catdesc       = [ catdesc,       amp_bins_vname, 'Histogram bins for wave amplitude.'];
	fieldnam      = [ fieldnam,      amp_bins_vname, 'Wave amplitude bins' ];
	fillval       = [ fillval,       amp_bins_vname, single(-1e31) ];
	format        = [ format,        amp_bins_vname, 'f11.4' ];
	lablaxis      = [ lablaxis,      amp_bins_vname, 'Amp bins' ];
	si_conversion = [ si_conversion, amp_bins_vname, '1e-9>T' ];
	units         = [ units,         amp_bins_vname, 'nT' ];
	validmin      = [ validmin,      amp_bins_vname, single(-10.0) ];
	validmax      = [ validmax,      amp_bins_vname, single(10.0) ];
	var_type      = [ var_type,      amp_bins_vname, 'support_data' ];
	
	% BINS_PHASE
	catdesc       = [ catdesc,       phase_bins_vname, 'Histogram bins for wave phase.'];
	fieldnam      = [ fieldnam,      phase_bins_vname, 'Wave amplitude bins' ];
	fillval       = [ fillval,       phase_bins_vname, single(-1e31) ];
	format        = [ format,        phase_bins_vname, 'f11.4' ];
	lablaxis      = [ lablaxis,      phase_bins_vname, 'Phase bins' ];
	si_conversion = [ si_conversion, phase_bins_vname, '1e0>degrees' ];
	units         = [ units,         phase_bins_vname, 'degrees' ];
	validmin      = [ validmin,      phase_bins_vname, single(-180.0) ];
	validmax      = [ validmax,      phase_bins_vname, single(180.0) ];
	var_type      = [ var_type,      phase_bins_vname, 'support_data' ];
	
	% BINS_PSD
	catdesc       = [ catdesc,       psd_bins_vname, 'Histogram bins for wave amplitude.' ];
	fieldnam      = [ fieldnam,      psd_bins_vname, 'Wave amplitude bins' ];
	fillval       = [ fillval,       psd_bins_vname, single(-1e31) ];
	format        = [ format,        psd_bins_vname, 'f11.4' ];
	lablaxis      = [ lablaxis,      psd_bins_vname, 'PSD bins' ];
	si_conversion = [ si_conversion, psd_bins_vname, '1e-18>T^2 Hz^-1' ];
	units         = [ units,         psd_bins_vname, 'nT^2 Hz^-1' ];
	validmin      = [ validmin,      psd_bins_vname, single(-10.0) ];
	validmax      = [ validmax,      psd_bins_vname, single(10.0) ];
	var_type      = [ var_type,      psd_bins_vname, 'support_data' ];
	
%-------------------------------------------------------------------------
% METADATA INFO //////////////////////////////////////////////////////////
%-------------------------------------------------------------------------
	
	% AMPLITUDE BIN LABEL
	catdesc       = [ catdesc,       amp_bin_labl_vname, 'Axis labels for amplitude bins' ];
	fieldnam      = [ fieldnam,      amp_bin_labl_vname, 'Amplitude labels' ];
	format        = [ format,        amp_bin_labl_vname, 'A7' ];
	var_type      = [ var_type,      amp_bin_labl_vname, 'metadata' ];
	
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
	catdesc       = [ catdesc,       flag_labl_vname, 'Axis label for histogram flag.' ];
	fieldnam      = [ fieldnam,      flag_labl_vname, 'Hist flag labels' ];
	format        = [ format,        flag_labl_vname, 'A13' ];
	var_type      = [ var_type,      flag_labl_vname, 'metadata' ];
	
	% PHASE BINS LABEL
	catdesc       = [ catdesc,       phase_bin_labl_vname, 'Axis label for phase bins.' ];
	fieldnam      = [ fieldnam,      phase_bin_labl_vname, 'Phase labels' ];
	format        = [ format,        phase_bin_labl_vname, 'A13' ];
	var_type      = [ var_type,      phase_bin_labl_vname, 'metadata' ];
	
	% PSD BINS LABEL
	catdesc       = [ catdesc,       psd_bin_labl_vname, 'Axis label for psd bins.' ];
	fieldnam      = [ fieldnam,      psd_bin_labl_vname, 'PSD labels' ];
	format        = [ format,        psd_bin_labl_vname, 'A13' ];
	var_type      = [ var_type,      psd_bin_labl_vname, 'metadata' ];
	
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
	spdfcdfwrite( fsm_cal_file,                       ...
	              var_list,                           ...
	              'GlobalAttributes',   global_attrs, ...
	              'VariableAttributes', var_attrs,    ...
	              'VarDatatypes',       vardatatypes, ...
	              'VarCompress',        varcompress,  ...
	              'RecordBound',        recbound,     ...
	              'Singleton',          singleton     ...
	            );
end