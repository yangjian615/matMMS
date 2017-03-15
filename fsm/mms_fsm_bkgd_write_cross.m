%
% Name
%   mms_fsm_bkgd_write_cross
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
function fsm_cal_file = mms_fsm_bkgd_write_cross(sc, mode, optdesc, tstart, data, varargin)

	% Global variables
	%   - See mms_fsm_init.m
	global cal_path_root data_path_root dropbox_root hk_root log_path_root unh_data_root

%------------------------------------%
% Check Inputs                       %
%------------------------------------%

	% Defaults
	tf_brst  = false;
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
	mag   = regexp(optdesc, '-', 'split');
	mag   = upper(mag{2});

%------------------------------------%
% Verify Data                        %
%------------------------------------%

	%
	% Check sizes
	%
	assert( isa(data.t,           'int64'),  'data.t should be int64.' );
	assert( isa(data.f,           'single'), 'data.f should be single.');
	assert( isa(data.gain,        'single'), 'data.gain should be single.');
	assert( isa(data.phase_shift, 'single'), 'data.phase_shift should be single.');
	assert( isa(data.psd_rat,     'single'), 'data.psd_rat should be single.');
	assert( isa(data.flag,        'uint8'),  'data.flag should be uint8.');
	
	% convert the start time to yyyy-mm-dd
%	tstart = MrTimeParser(tstart, '%Y-%M-%dT%H:%m:%S', '%Y%M%d');

%------------------------------------%
% Version History                    %
%------------------------------------%
	
	% Describe the modifications to each version
	mods    = {   'v0.0.0 -- First version.' ...
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
	                       'Logical_source_description', ['Level 3 Fluxgate-Searchcoil Merged Magnetic Field Calibration Data.'], ...
	                       'Mission_group',              'MMS', ...
	                       'PI_affiliation',             'UNH, LPP, IWF, UCLA, UCLA, IWF, UNH', ...
	                       'PI_name',                    ['R.B. Torbert, O. LeContel, W. Magnes, C.T. Russell, R. Strangeway, ', ...
	                                                      'D. Fischer, M.R. Argall.'], ...
	                       'Project',                    'STP>Solar Terrestrial Physics', ...
	                       'Source_name',                source_name, ...
	                       'TEXT',                       ['Calibration data used while created the FSM data product. ' ...
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
	f_vname         = [prefix 'f'                   suffix];
	c_index_vname   = [prefix 'component'  '_index' suffix];
	gain_vname      = [prefix 'gain'       '_omb'   suffix];
	offset_vname    = [prefix 'phaseshift' '_omb'   suffix];
	psdrat_vname    = [prefix 'psdrat'     '_omb'   suffix];
	flag_vname      = [prefix 'flag'                suffix];
	flag_hist_vname = [prefix 'flag'       '_hist'  suffix];
	
	% Histogram
	gain_hist_vname   = [prefix 'gain'       '_hist'  suffix];
	offset_hist_vname = [prefix 'phaseshift' '_hist'  suffix];
	psdrat_hist_vname = [prefix 'psdrat'     '_hist'  suffix];
	
	% Histogram bins
	gain_bins_vname   = [prefix 'gain'       '_bins'  suffix];
	offset_bins_vname = [prefix 'phaseshift' '_bins'  suffix];
	psdrat_bins_vname = [prefix 'psdrat'     '_bins'  suffix];
	
	% Metadata
	c_labl_vname           = [prefix 'component'  '_labl'      suffix];
	gain_bins_labl_vname   = [prefix 'gain'       '_bins_labl' suffix];
	offset_bins_labl_vname = [prefix 'phaseshift' '_bins_labl' suffix];
	psdrat_bins_labl_vname = [prefix 'psdrat'     '_bins_labl' suffix];
	f_labl_vname           = [prefix 'f'          '_labl'      suffix];
	flag_hist_labl_vname   = [prefix 'flag'       '_hist_labl' suffix];

%------------------------------------------------------
% Metadata: Labels & Indices                          |
%------------------------------------------------------
	
	% Component label & index
	%   - CDF is 0-based
	c_labl  = {'X', 'Y', 'Z'};
	c_index = uint8(0:2);
	
	% Frequency Labels
	len    = length(data.f);
	powten = floor( log10( len ) ) + 1;
	fmt    = sprintf('%%0%ii', powten);
	f_labl = strcat( {'freq'}, num2str( [1:len]', fmt) )';
	
	% Gain bins Label
	len            = length(data.gain_bins);
	powten         = floor( log10( len ) ) + 1;
	fmt            = sprintf('%%0%ii', powten);
	gain_bins_labl = strcat( {'amp'}, num2str( [1:len]', fmt) )';
	
	% Hist Flag Label
	len            = length(data.hist_flag);
	powten         = floor( log10( len ) ) + 1;
	fmt            = sprintf('%%0%ii', powten);
	flag_hist_labl = strcat( {'histflag'}, num2str( [1:len]', fmt) )';
	
	% Phase bins Label
	len             = length(data.phase_shift_bins);
	powten          = floor( log10( len ) ) + 1;
	fmt             = sprintf('%%0%ii', powten);
	offset_bins_labl = strcat( {'phase'}, num2str( [1:len]', fmt) )';
	
	% PSD bins Label
	len              = length(data.psd_rat_bins);
	powten           = floor( log10( len ) ) + 1;
	fmt              = sprintf('%%0%ii', powten);
	psdrat_bins_labl = strcat( {'psd'}, num2str( [1:len]', fmt) )';

%------------------------------------------------------
% Variable Data                                       |
%------------------------------------------------------

	% Variable Data
	%   - Order as [ DEPEND_1,  DEPEND_2,  DEPEND_0 ]
	%   - Order as [ COMPONENT, FREQUENCY, TIME     ]
	var_list = { t_vname,                data.t,                ...
	             f_vname,                data.f,                ...
	             c_index_vname,          c_index,               ...
	             gain_vname,             permute( data.gain,        [3,1,2] ), ...
	             offset_vname,           permute( data.phase_shift, [3,1,2] ), ...
	             psdrat_vname,           permute( data.psd_rat,     [3,1,2] ), ...
	             flag_vname,             data.flag,             ...
	             gain_bins_labl_vname,   gain_bins_labl,        ...
	             offset_bins_labl_vname, offset_bins_labl,      ...
	             psdrat_bins_labl_vname, psdrat_bins_labl,      ...
	             c_labl_vname,           c_labl,                ...
	             f_labl_vname,           f_labl,                ...
	             flag_hist_labl_vname,   flag_hist_labl,        ...
	             flag_hist_vname,        data.hist_flag,        ...
	             gain_hist_vname,        permute( data.gain_hist,        [3,2,1,4] ), ...
	             offset_hist_vname,      permute( data.phase_shift_hist, [3,2,1,4] ), ...
	             psdrat_hist_vname,      permute( data.psd_rat_hist,     [3,2,1,4] ), ...
	             gain_bins_vname,        data.gain_bins,        ...
	             offset_bins_vname,      data.phase_shift_bins, ...
	             psdrat_bins_vname,      data.psd_rat_bins      ...
	           };
	clear data

	% Record Variance
	recbound = { t_vname,           ...
	             flag_vname,        ...
	             gain_vname,        ...
	             offset_vname,      ...
	             psdrat_vname,      ...
	             gain_hist_vname,   ...
	             offset_hist_vname, ...
	             psdrat_hist_vname, ...
	           };

	% Data types
	vardatatypes = { t_vname,                'cdf_time_tt2000', ...
	                 f_vname,                'cdf_float', ...
	                 c_index_vname,          'cdf_uint1', ...
	                 flag_vname,             'cdf_uint1', ...
	                 gain_bins_labl_vname,   'cdf_char',  ...
	                 offset_bins_labl_vname, 'cdf_char',  ...
	                 psdrat_bins_labl_vname, 'cdf_char',  ...
	                 c_labl_vname,           'cdf_char',  ...
	                 f_labl_vname,           'cdf_char',  ...
	                 flag_hist_labl_vname,   'cdf_char',  ...
	                 flag_hist_vname,        'cdf_uint1', ...
	                 gain_vname,             'cdf_float', ...
	                 offset_vname,           'cdf_float', ...
	                 psdrat_vname,           'cdf_float', ...
	                 gain_hist_vname,        'cdf_uint4', ...
	                 offset_hist_vname,      'cdf_uint4', ...
	                 psdrat_hist_vname,      'cdf_uint4', ...
	                 gain_bins_vname,        'cdf_float', ...
	                 offset_bins_vname,      'cdf_float', ...
	                 psdrat_bins_vname,      'cdf_float'  ...
	               };
	
	% CDF Compression Level
	varcompress = { flag_vname,        'gzip.6', ...
	                flag_hist_vname,   'gzip.6', ...
	                gain_vname,        'gzip.6', ...
	                offset_vname,      'gzip.6', ...
	                psdrat_vname,      'gzip.6', ...
	                gain_hist_vname,   'gzip.6', ...
	                offset_hist_vname, 'gzip.6', ...
	                psdrat_hist_vname, 'gzip.6', ...
	                gain_bins_vname    'gzip.6', ...
	                offset_bins_vname  'gzip.6', ...
	                psdrat_bins_vname  'gzip.6'  ...
	              };
	
	singleton = { gain_vname,        ...
	              offset_vname,      ...
	              psdrat_vname,      ...
	              gain_hist_vname,   ...
	              offset_hist_vname, ...
	              psdrat_hist_vname  ...
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
	
	% GAIN
	catdesc       = [ catdesc,       gain_vname, ['Ratio of ' mag ' to SCM wave amplitude.'] ];
	depend_0      = [ depend_0,      gain_vname, t_vname ];
	depend_1      = [ depend_1,      gain_vname, f_vname ];
	display_type  = [ display_type,  gain_vname, 'spectrogram'];
	fieldnam      = [ fieldnam,      gain_vname, 'Gain' ];
	fillval       = [ fillval,       gain_vname, single(-1e31) ];
	format        = [ format,        gain_vname, 'f11.4' ];
	labl_ptr_1    = [ labl_ptr_1,    gain_vname, c_labl_vname ];
	si_conversion = [ si_conversion, gain_vname, ' ' ];
	units         = [ units,         gain_vname, ' ' ];
	validmin      = [ validmin,      gain_vname, single(0.0) ];
	validmax      = [ validmax,      gain_vname, single(1e4) ];
	var_type      = [ var_type,      gain_vname, 'data' ];
	
	% PHASE OFFSET
	catdesc       = [ catdesc,       offset_vname, ['Difference of ' mag ' to SCM phase.'] ];
	depend_0      = [ depend_0,      offset_vname, t_vname ];
	depend_1      = [ depend_1,      offset_vname, f_vname ];
	display_type  = [ display_type,  offset_vname, 'spectrogram'];
	fieldnam      = [ fieldnam,      offset_vname, 'Phase Diff.' ];
	fillval       = [ fillval,       offset_vname, single(-1e31) ];
	format        = [ format,        offset_vname, 'f11.4' ];
	labl_ptr_1    = [ labl_ptr_1,    offset_vname, c_labl_vname ];
	si_conversion = [ si_conversion, offset_vname, '1e0>degrees' ];
	units         = [ units,         offset_vname, 'degrees' ];
	validmin      = [ validmin,      offset_vname, single(-360.0) ];
	validmax      = [ validmax,      offset_vname, single(360.0) ];
	var_type      = [ var_type,      offset_vname, 'data' ];
	
	% PSD RATIO
	catdesc       = [ catdesc,       psdrat_vname, ['Ratio of ' mag ' to SCM power spectral density.'] ];
	depend_0      = [ depend_0,      psdrat_vname, t_vname ];
	depend_1      = [ depend_1,      psdrat_vname, f_vname ];
	display_type  = [ display_type,  psdrat_vname, 'spectrogram'];
	fieldnam      = [ fieldnam,      psdrat_vname, 'PSD' ];
	fillval       = [ fillval,       psdrat_vname, single(-1e31) ];
	format        = [ format,        psdrat_vname, 'f11.4' ];
	labl_ptr_1    = [ labl_ptr_1,    psdrat_vname, c_labl_vname ];
	si_conversion = [ si_conversion, psdrat_vname, '1e-18>T^2/Hz' ];
	units         = [ units,         psdrat_vname, 'nT^2/Hz' ];
	validmin      = [ validmin,      psdrat_vname, single(1e4) ];
	validmax      = [ validmax,      psdrat_vname, single(0) ];
	var_type      = [ var_type,      psdrat_vname, 'data' ];
	
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
	
	% HIST_FLAG
	catdesc       = [ catdesc,       flag_hist_vname, 'Bitflag indicating regime. Bit=unset/set: 1=slow/fast(brst), 2=hi/lo-range, 3=Deck32/64' ];
	fieldnam      = [ fieldnam,      flag_hist_vname, 'Histogram bitflag' ];
	fillval       = [ fillval,       flag_hist_vname, uint32(255) ];
	format        = [ format,        flag_hist_vname, 'I3' ];
	lablaxis      = [ lablaxis,      flag_hist_vname, 'Hist flag' ];
	validmin      = [ validmin,      flag_hist_vname, uint8(0) ];
	validmax      = [ validmax,      flag_hist_vname, uint8(254) ];
	var_type      = [ var_type,      flag_hist_vname, 'support_data' ];
	
	% GAIN
	catdesc       = [ catdesc,       gain_hist_vname, ['Ratio of ' mag ' to SCM wave amplitude, histogrammed in time at each frequency.'] ];
	depend_0      = [ depend_0,      gain_hist_vname, f_vname ];
	depend_1      = [ depend_1,      gain_hist_vname, gain_bins_vname ];
	depend_2      = [ depend_2,      gain_hist_vname, flag_hist_vname ];
	depend_3      = [ depend_3,      gain_hist_vname, c_index_vname ];
	display_type  = [ display_type,  gain_hist_vname, 'spectrogram'];
	fieldnam      = [ fieldnam,      gain_hist_vname, 'Histogrammed gain' ];
	fillval       = [ fillval,       gain_hist_vname, uint32(65535) ];
	format        = [ format,        gain_hist_vname, 'I5' ];
	labl_ptr_1    = [ labl_ptr_1,    gain_hist_vname, c_labl_vname ];
	si_conversion = [ si_conversion, gain_hist_vname, ' ' ];
	units         = [ units,         gain_hist_vname, ' ' ];
	validmin      = [ validmin,      gain_hist_vname, uint32(0) ];
	validmax      = [ validmax,      gain_hist_vname, uint32(65533) ];
	var_type      = [ var_type,      gain_hist_vname, 'data' ];
	
	% PHASE SHIFT
	catdesc       = [ catdesc,       offset_hist_vname, ['Difference between ' mag ' to SCM wave phase, histogrammed in time at each frequency.'] ];
	depend_0      = [ depend_0,      offset_hist_vname, f_vname ];
	depend_1      = [ depend_1,      offset_hist_vname, offset_bins_vname ];
	depend_2      = [ depend_2,      offset_hist_vname, flag_hist_vname ];
	depend_3      = [ depend_3,      offset_hist_vname, c_index_vname ];
	display_type  = [ display_type,  offset_hist_vname, 'spectrogram'];
	fieldnam      = [ fieldnam,      offset_hist_vname, 'Histogrammed phase shift' ];
	fillval       = [ fillval,       offset_hist_vname, uint32(65535) ];
	format        = [ format,        offset_hist_vname, 'I5' ];
	labl_ptr_1    = [ labl_ptr_1,    offset_hist_vname, c_labl_vname ];
	si_conversion = [ si_conversion, offset_hist_vname, '1e0>degrees' ];
	units         = [ units,         offset_hist_vname, 'degrees' ];
	validmin      = [ validmin,      offset_hist_vname, uint32(0) ];
	validmax      = [ validmax,      offset_hist_vname, uint32(65533) ];
	var_type      = [ var_type,      offset_hist_vname, 'data' ];
	
	% PSDRAT
	catdesc       = [ catdesc,       psdrat_hist_vname, ['Ratio of ' mag ' to SCM power spectral density, histogrammed in time at each frequency.'] ];
	depend_0      = [ depend_0,      psdrat_hist_vname, f_vname ];
	depend_1      = [ depend_1,      psdrat_hist_vname, psdrat_bins_vname ];
	depend_2      = [ depend_2,      psdrat_hist_vname, flag_hist_vname ];
	depend_3      = [ depend_3,      psdrat_hist_vname, c_index_vname ];
	display_type  = [ display_type,  psdrat_hist_vname, 'spectrogram'];
	fieldnam      = [ fieldnam,      psdrat_hist_vname, 'Histogrammed phase shift' ];
	fillval       = [ fillval,       psdrat_hist_vname, uint32(65535) ];
	format        = [ format,        psdrat_hist_vname, 'I5' ];
	labl_ptr_1    = [ labl_ptr_1,    psdrat_hist_vname, c_labl_vname ];
	si_conversion = [ si_conversion, psdrat_hist_vname, '1e0>degrees' ];
	units         = [ units,         psdrat_hist_vname, 'degrees' ];
	validmin      = [ validmin,      psdrat_hist_vname, uint32(0) ];
	validmax      = [ validmax,      psdrat_hist_vname, uint32(65533) ];
	var_type      = [ var_type,      psdrat_hist_vname, 'data' ];
	
%-------------------------------------------------------------------------
% HISTOGRAM BINS INFO ////////////////////////////////////////////////////
%-------------------------------------------------------------------------
	
	% GAIN BINS
	catdesc       = [ catdesc,       gain_bins_vname, 'Histogram bins for gain factor.' ];
	fieldnam      = [ fieldnam,      gain_bins_vname, 'Gain bins' ];
	fillval       = [ fillval,       gain_bins_vname, single(-1e31) ];
	format        = [ format,        gain_bins_vname, 'f11.4' ];
	si_conversion = [ si_conversion, gain_bins_vname, ' ' ];
	units         = [ units,         gain_bins_vname, ' ' ];
	validmin      = [ validmin,      gain_bins_vname, single(-10.0) ];
	validmax      = [ validmax,      gain_bins_vname, single(10.0) ];
	var_type      = [ var_type,      gain_bins_vname, 'support_data' ];
	
	% PHASE SHIFT BINS
	catdesc       = [ catdesc,       offset_bins_vname, 'Histogram bins for phase shift.' ];
	fieldnam      = [ fieldnam,      offset_bins_vname, 'Phase shift bins' ];
	fillval       = [ fillval,       offset_bins_vname, single(-1e31) ];
	format        = [ format,        offset_bins_vname, 'f11.4' ];
	si_conversion = [ si_conversion, offset_bins_vname, '1e0>degrees' ];
	units         = [ units,         offset_bins_vname, 'degrees' ];
	validmin      = [ validmin,      offset_bins_vname, single(-180.0) ];
	validmax      = [ validmax,      offset_bins_vname, single(180.0) ];
	var_type      = [ var_type,      offset_bins_vname, 'support_data' ];
	
	% PSD RATIO BINS
	catdesc       = [ catdesc,       psdrat_bins_vname, 'Histogram bins for PSD ratio.' ];
	fieldnam      = [ fieldnam,      psdrat_bins_vname, 'PSD ratio bins' ];
	fillval       = [ fillval,       psdrat_bins_vname, single(-1e31) ];
	format        = [ format,        psdrat_bins_vname, 'f11.4' ];
	si_conversion = [ si_conversion, psdrat_bins_vname, ' ' ];
	units         = [ units,         psdrat_bins_vname, ' ' ];
	validmin      = [ validmin,      psdrat_bins_vname, single(-10.0) ];
	validmax      = [ validmax,      psdrat_bins_vname, single(10.0) ];
	var_type      = [ var_type,      psdrat_bins_vname, 'support_data' ];
	
%-------------------------------------------------------------------------
% METADATA INFO //////////////////////////////////////////////////////
%-------------------------------------------------------------------------
	
	% GAIN BIN LABEL
	catdesc       = [ catdesc,       gain_bins_labl_vname, 'Axis labels for gain bins' ];
	fieldnam      = [ fieldnam,      gain_bins_labl_vname, 'Gain labels' ];
	format        = [ format,        gain_bins_labl_vname, 'A7' ];
	var_type      = [ var_type,      gain_bins_labl_vname, 'metadata' ];
	
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
	catdesc       = [ catdesc,       flag_hist_labl_vname, 'Axis label for histogram flag.' ];
	fieldnam      = [ fieldnam,      flag_hist_labl_vname, 'Hist flag labels' ];
	format        = [ format,        flag_hist_labl_vname, 'A13' ];
	var_type      = [ var_type,      flag_hist_labl_vname, 'metadata' ];
	
	% PHASE OFFSET BINS LABEL
	catdesc       = [ catdesc,       offset_bins_labl_vname, 'Axis label for phase offset bins.' ];
	fieldnam      = [ fieldnam,      offset_bins_labl_vname, 'Phase offset labels' ];
	format        = [ format,        offset_bins_labl_vname, 'A13' ];
	var_type      = [ var_type,      offset_bins_labl_vname, 'metadata' ];
	
	% PSD RATIO BINS LABEL
	catdesc       = [ catdesc,       psdrat_bins_labl_vname, 'Axis label for psd-ratio bins.' ];
	fieldnam      = [ fieldnam,      psdrat_bins_labl_vname, 'PSD-ratio labels' ];
	format        = [ format,        psdrat_bins_labl_vname, 'A13' ];
	var_type      = [ var_type,      psdrat_bins_labl_vname, 'metadata' ];
	
%-------------------------------------------------------------------------
% Collect Variable Attributes ////////////////////////////////////////////
%-------------------------------------------------------------------------
	% Collect into a structure
	var_attrs = struct( 'CATDESC',       { catdesc       }, ...
	                    'DEPEND_0',      { depend_0      }, ...
	                    'DEPEND_1',      { depend_1      }, ...
	                    'DEPEND_2',      { depend_2      }, ...
	                    'DEPEND_3',      { depend_3      }, ...
	                    'DISPLAY_TYPE',  { display_type  }, ...
	                    'FIELDNAM',      { fieldnam      }, ...
	                    'FILLVAL',       { fillval       }, ...
	                    'FORMAT',        { format        }, ...
	                    'LABLAXIS',      { lablaxis      }, ...
	                    'LABL_PTR_1',    { labl_ptr_1    }, ...
	                    'LABL_PTR_2',    { labl_ptr_2    }, ...
	                    'LABL_PTR_3',    { labl_ptr_3    }, ...
	                    'SI_CONVERSION', { si_conversion }, ...
	                    'TIME_BASE',     { time_base     }, ...
	                    'UNITS',         { units         }, ...
	                    'VALIDMIN',      { validmin      }, ...
	                    'VALIDMAX',      { validmax      }, ...
	                    'VAR_TYPE',      { var_type      }  ...
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