%
% Name
%   mms_fsm_nf_write
%
% Purpose
%   Create a MATLAB save file of inputs needed for Bestarg.
%
% Calling Sequence
%   EDI_QL_FILE = mms_fsm_nf_write(EDI_QL)
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
function nf_file = mms_fsm_nf_write(sc, mode, optdesc, tstart, data, varargin)

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
	instr     = 'fsm';
	level     = 'l2plus';
	mag_instr = regexp( optdesc, '-', 'split' );
	mag_instr = mag_instr{2};

%------------------------------------%
% Verify Data                        %
%------------------------------------%

	%
	% Check sizes
	%
	assert( isa(data.t,    'int64'),  'data.t should be int64.' );
	assert( isa(data.dt,   'int64'),  'data.dt should be int64.' );
	assert( isa(data.f,    'single'), 'data.f should be single.' );
	assert( isa(data.flag, 'uint8'),  'data.flag should be uint8.');
	assert( isa(data.nf,   'single'), 'data.nf should be single.');
	assert( isa(data.std,  'single'), 'data.std should be single.');

%------------------------------------%
% Version History                    %
%------------------------------------%
	
	% Describe the modifications to each version
	mods    = { 'v0.0.0 -- First version.' ...
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
	nf_file = mms_construct_filename(sc, instr, mode, level, ...
	                                 'TStart',  tstart,      ...
	                                 'Version', version,     ...
	                                 'OptDesc', optdesc);

	% Find the latest z-version
	%   - Look in both DROPBOX and DATA_PATH
	vz = mms_latest_zversion(dropbox_root, nf_file, 'RootDir', data_path_root);

	% Reform the file name
	version  = [vx '.' vy '.' sprintf('%0i', vz)];
	nf_file = mms_construct_filename(sc, instr, mode, level, ...
	                                 'TStart',  tstart,      ...
	                                 'Version', version,     ...
	                                 'OptDesc', optdesc);
	nf_file = fullfile(dropbox_root, nf_file);

	% Notify where file is located
	mrfprintf('logtext', ['Creating FSM NF file at "' nf_file '".']);

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
	[~, logical_file_id, ext] = fileparts(nf_file);
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
	                       'Logical_source_description', ['Level 2-Plus Fluxgate-Searchcoil Merged Magnetic Field Calibration ' ...
	                                                      'noise floor data from ' upper(mag_instr) '.'], ...
	                       'Mission_group',              'MMS', ...
	                       'PI_affiliation',             'UNH, LPP, IWF, UCLA, UCLA, IWF, UNH', ...
	                       'PI_name',                    ['R.B. Torbert, O. LeContel, W. Magnes, C.T. Russell, R. Strangeway, ', ...
	                                                      'D. Fischer, M.R. Argall.'], ...
	                       'Project',                    'STP>Solar Terrestrial Physics', ...
	                       'Source_name',                source_name, ...
	                       'TEXT',                       ['Noise floor data used while created the FSM data product. ' ...
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
	
	% Data
	t_vname       = 'Epoch';
	dt_vname      = [prefix 'dt'            suffix];
	c_index_vname = [prefix 'comp' '_index' suffix];
	f_vname       = [prefix 'f'             suffix];
	flag_vname    = [prefix 'flag'          suffix];
	nf_vname      = [prefix 'nf'            suffix];
	std_vname     = [prefix 'std'           suffix];
	
	% Metadata
	c_labl_vname    = [prefix 'comp' '_labl' suffix];
	f_labl_vname    = [prefix 'f'    '_labl' suffix];
	flag_labl_vname = [prefix 'flag' '_labl' suffix];

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
	
	% Hist Flag Label
	len       = length(data.flag);
	powten    = floor( log10( len ) ) + 1;
	fmt       = sprintf('%%0%ii', powten);
	flag_labl = strcat( {'flag'}, num2str( [1:len]', fmt) )';

%------------------------------------------------------
% Variable Data                                       |
%------------------------------------------------------

	% Variable Data
	%   - Order as [ DEPEND_1,  DEPEND_2,  DEPEND_0 ]
	%   - Order as [ COMPONENT, FREQUENCY, TIME     ]
	var_list = { t_vname,         data.t,    ...
	             dt_vname,        data.dt,   ...
	             f_vname,         data.f,    ...
	             c_index_vname,   c_index,   ...
	             flag_vname,      data.flag, ...
	             nf_vname,        data.nf,   ...
	             std_vname,       data.std,  ...
	             c_labl_vname,    c_labl,    ...
	             f_labl_vname,    f_labl,    ...
	             flag_labl_vname, flag_labl  ...
	           };

	% Record Variance
	%   - If data is for a single record, NTIMES=1
	%   - NF = [nFreqs, nComp, nFlags, nTimes]
	%   - MATLAB removes trailing shallow dimensions
	%   - A work-around is to remove record variance
	if ndims( data.nf ) == 4
		recbound = { t_vname,    ...
		             flag_vname, ...
		             nf_vname,   ...
		             std_vname   ...
		           };
		
		singleton = { t_vname,    ...
		              flag_vname, ...
		              nf_vname,   ...
		              std_vname   ...
		            };
		
	% Variables always have at least 2 dimensions
	%   - FLAG is [nFlags, nTimes] = Nx1
	%   - If we set singleton, it will become Nx1x1
	%   - The FLAG dimension is the one we want to be singleton
	else
		recbound = {};
		singleton = { t_vname,    ...
		              nf_vname,   ...
		              std_vname   ...
		            };
	end

	% Data types
	vardatatypes = { t_vname,         'cdf_time_tt2000', ...
	                 dt_vname,        'cdf_int8',  ...
	                 f_vname,         'cdf_float', ...
	                 c_index_vname,   'cdf_uint1', ...
	                 flag_vname,      'cdf_uint1', ...
	                 nf_vname,        'cdf_float', ...
	                 std_vname,       'cdf_float', ...
	                 c_labl_vname,    'cdf_char',  ...
	                 f_labl_vname,    'cdf_char',  ...
	                 flag_labl_vname, 'cdf_char'   ...
	               };
	
	% CDF Compression Level
	varcompress = { f_vname,    'gzip.6', ...
	                flag_vname, 'gzip.6', ...
	                nf_vname,   'gzip.6', ...
	                std_vname,  'gzip.6'  ...
	              };
	
	% Free up memory
	clear data

%------------------------------------------------------
% Variable Attributes                                 |
%------------------------------------------------------
	catdesc         = {};
	depend_0        = {};
	depend_1        = {};
	depend_2        = {};
	depend_3        = {};
	delta_plus_var  = {};
	delta_minus_var = {};
	display_type    = {};
	fieldnam        = {};
	fillval         = {};
	format          = {};
	lablaxis        = {};
	labl_ptr_1      = {};
	labl_ptr_2      = {};
	labl_ptr_3      = {};
	si_conversion   = {};
	time_base       = {};
	units           = {};
	validmin        = {};
	validmax        = {};
	var_type        = {};
	
	% EPOCH
	catdesc        = [ catdesc,        t_vname, 'Time tags for the merged magnetic field data.' ];
	delta_plus_var = [ delta_plus_var, t_vname, dt_vname ];
	fieldnam       = [ fieldnam,       t_vname, 'Time' ];
	fillval        = [ fillval,        t_vname, spdfcomputett2000([9999 12 31 23 59 59 999 999 999]) ];
	format         = [ format,         t_vname, 'I16' ];
	lablaxis       = [ lablaxis,       t_vname, 'TAI' ];
	si_conversion  = [ si_conversion,  t_vname, '1e-9>s' ];
	time_base      = [ time_base,      t_vname, 'J2000' ];
	units          = [ units,          t_vname, 'ns' ];
	validmin       = [ validmin,       t_vname, spdfcomputett2000([2015 3 12 0 0 0 0 0 0]) ];
	validmax       = [ validmax,       t_vname, spdfcomputett2000([2050 1  1 0 0 0 0 0 0]) ];
	var_type       = [ var_type,       t_vname, 'support_data' ];
	
	% DT
	catdesc       = [ catdesc,       dt_vname, 'Time tags for the merged magnetic field data.' ];
	fieldnam      = [ fieldnam,      dt_vname, 'Time' ];
	fillval       = [ fillval,       dt_vname, int64(  4294967295 ) ];
	format        = [ format,        dt_vname, 'I16' ];
	lablaxis      = [ lablaxis,      dt_vname, 'TAI' ];
	si_conversion = [ si_conversion, dt_vname, '1e-9>s' ];
	units         = [ units,         dt_vname, 'ns' ];
	validmin      = [ validmin,      dt_vname, int64( 31.0 * 86400.0 * 1e9 ) ];
	validmax      = [ validmax,      dt_vname, int64(  7.0 * 86400.0 * 1e9 ) ];
	var_type      = [ var_type,      dt_vname, 'support_data' ];
	
	% Frequency
	catdesc       = [ catdesc,       f_vname, 'Center frequency bins.' ];
	fieldnam      = [ fieldnam,      f_vname, 'Frequency' ];
	fillval       = [ fillval,       f_vname, single(-1e31) ];
	format        = [ format,        f_vname, 'f6.1' ];
	lablaxis      = [ lablaxis,      f_vname, 'f' ];
	si_conversion = [ si_conversion, f_vname, '1e0>Hz' ];
	units         = [ units,         f_vname, 'Hz' ];
	validmin      = [ validmin,      f_vname, single(0.0) ];
	validmax      = [ validmax,      f_vname, single(4096.0) ];
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
	catdesc       = [ catdesc,       flag_vname, 'Instrument operations flag. Bit values' ];
	depend_0      = [ depend_0,      flag_vname, t_vname ];
	display_type  = [ display_type,  flag_vname, 'time_series'];
	fieldnam      = [ fieldnam,      flag_vname, 'Flag' ];
	fillval       = [ fillval,       flag_vname, uint8(255) ];
	format        = [ format,        flag_vname, 'I2' ];
	lablaxis      = [ lablaxis,      flag_vname, 'Flag' ];
	validmin      = [ validmin,      flag_vname, uint8(0) ];
	validmax      = [ validmax,      flag_vname, uint8(254) ];
	var_type      = [ var_type,      flag_vname, 'support_data' ];
	
	% NF
	catdesc         = [ catdesc,         nf_vname, ['Standard deviation of the noise floor of the ' mag_instr ' instrument.'] ];
	delta_plus_var  = [ delta_plus_var,  nf_vname, std_vname ];
	delta_minus_var = [ delta_minus_var, nf_vname, std_vname ];
	depend_0        = [ depend_0,        nf_vname, t_vname ];
	depend_1        = [ depend_1,        nf_vname, f_vname ];
	depend_2        = [ depend_2,        nf_vname, c_index_vname ];
	depend_3        = [ depend_3,        nf_vname, flag_vname ];
	display_type    = [ display_type,    nf_vname, 'time_series'];
	fieldnam        = [ fieldnam,        nf_vname, 'Standard deviation' ];
	fillval         = [ fillval,         nf_vname, single(-1e31) ];
	format          = [ format,          nf_vname, 'f11.4' ];
	labl_ptr_1      = [ labl_ptr_1,      nf_vname, f_labl_vname ];
	labl_ptr_2      = [ labl_ptr_2,      nf_vname, c_labl_vname ];
	labl_ptr_3      = [ labl_ptr_3,      nf_vname, flag_labl_vname ];
	si_conversion   = [ si_conversion,   nf_vname, '1e-9>T' ];
	units           = [ units,           nf_vname, 'nT/sqrt(Hz)' ];
	validmin        = [ validmin,        nf_vname, single(0) ];
	validmax        = [ validmax,        nf_vname, single(1e5) ];
	var_type        = [ var_type,        nf_vname, 'data' ];

%-------------------------------------------------------------------------
% METADATA INFO //////////////////////////////////////////////////////
%-------------------------------------------------------------------------
	
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
	
	% FLAG LABEL
	catdesc       = [ catdesc,       flag_labl_vname, 'Axis label for histogram flag.' ];
	fieldnam      = [ fieldnam,      flag_labl_vname, 'Hist flag labels' ];
	format        = [ format,        flag_labl_vname, 'A13' ];
	var_type      = [ var_type,      flag_labl_vname, 'metadata' ];
	
%-------------------------------------------------------------------------
% Collect Variable Attributes ////////////////////////////////////////////
%-------------------------------------------------------------------------
	% Collect into a structure
	var_attrs = struct( 'CATDESC',         { catdesc },         ...
	                    'DELTA_MINUS_VAR', { delta_minus_var }, ...
	                    'DELTA_PLUS_VAR',  { delta_plus_var },  ...
	                    'DEPEND_0',        { depend_0 },        ...
	                    'DEPEND_1',        { depend_1 },        ...
	                    'DEPEND_2',        { depend_2 },        ...
	                    'DEPEND_3',        { depend_3 },        ...
	                    'DISPLAY_TYPE',    { display_type },    ...
	                    'FIELDNAM',        { fieldnam },        ...
	                    'FILLVAL',         { fillval },         ...
	                    'FORMAT',          { format },          ...
	                    'LABLAXIS',        { lablaxis },        ...
	                    'LABL_PTR_1',      { labl_ptr_1 },      ...
	                    'LABL_PTR_2',      { labl_ptr_2 },      ...
	                    'LABL_PTR_3',      { labl_ptr_3 },      ...
	                    'SI_CONVERSION',   { si_conversion },   ...
	                    'TIME_BASE',       { time_base },       ...
	                    'UNITS',           { units },           ...
	                    'VALIDMIN',        { validmin },        ...
	                    'VALIDMAX',        { validmax },        ...
	                    'VAR_TYPE',        { var_type }         ...
	                  );

%------------------------------------------------------
% Write the File                                      |
%------------------------------------------------------
	spdfcdfwrite( nf_file,                            ...
	              var_list,                           ...
	              'GlobalAttributes',   global_attrs, ...
	              'VariableAttributes', var_attrs,    ...
	              'VarDatatypes',       vardatatypes, ...
	              'VarCompress',        varcompress,  ...
	              'RecordBound',        recbound,     ...
	              'Singleton',          singleton     ...
	            );
end