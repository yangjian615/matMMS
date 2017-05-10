%--------------------------------------------------------------------------
% NAME
%   mms_fsm_l3_write
%
% PURPOSE
%   Write L3 merged magnetometer data to a file.
%
% CALLING SEQUENCE:
%   FSM_FILE = mms_fsm_l3_write( SC MODE TSTART DATA )
%       Write merged magnetic field data DATA for MMS spacecraft SC in data rate mode MODE
%       with start time TSTART to a CDF file named FSM_FILE.
%
%   FSM_FILE = mms_fsm_l3_write( ..., 'ParamName', ParamValue)
%       Also accept any of the parameter name-value pairs given below.
%
% Inputs:
%   SC:         in, required, type=char
%   MODE:       in, required, type=char
%   TSTART:     in, required, type=char
%   DATA:       in, required, type=struct
%   'OptDesc':  in, optional, type=struct
%               An optional descriptor to be included in the file name.
%   'Parents':  in, optional, type=struct
%               The parent files from which source and/or dependent data was retrieved.
%
% History:
%  2017-03-15  - Written by Matthew Argall
%
%--------------------------------------------------------------------------
function fsm_file = mms_fsm_l2plus_write( sc, mode, tstart, data, varargin )

	% Global variables
	%   - See mms_fsm_init.m
	global cal_path_root data_path_root dropbox_root hk_root log_path_root unh_data_root

%------------------------------------%
% Check Inputs                       %
%------------------------------------%

	% Defaults
	tf_empty = false;
	parents  = {};
	optdesc  = '';

	% Check inputs
	nOptArgs = length(varargin);
	for ii = 1 : 2 : nOptArgs
		switch varargin{ii}
			case 'OptDesc'
				optdesc = varargin{ii+1};
			case 'Parents'
				parents = varargin{ii+1};
			otherwise
				error(['Optional argument not recognized: "' varargin{ii} '".']);
		end
	end
	
	% Constants
	instr   = 'fsm';
	level   = 'l3';

%------------------------------------%
% Verify Data                        %
%------------------------------------%

	%
	% Check sizes
	%
	assert( isa(data.tt2000, 'int64'),  'data.tt2000 should be int64.' );
	assert( isa(data.b,      'single'), 'data.b should be single.');
	assert( isa(data.b_omb,  'single'), 'data.b_bcs should be single.');

%------------------------------------%
% Version History                    %
%------------------------------------%
	
	% Describe the modifications to each version
	mods    = { 'v0.0.0 -- First version.', ...
	            'v0.1.0 -- Added |B|. Updated DISPLAY_TYPE and FIELDNAM var attrs.', ...
	            'v1.0.0 -- Data is merged in time domain with "step" filter.', ...
	            'v1.1.0 -- Fix filter delay and settling time.' ...
	          };

	% Extract the version numbers
	parts = regexp( mods{end}, '^v([0-9]+)\.([0-9]+)\.([0-9]+)', 'tokens' );
	vx    = parts{1}{1};
	vy    = parts{1}{2};
	vz    = parts{1}{3};

%------------------------------------%
% Create Output File Name            %
%------------------------------------%
	%
	% TODO: Link version number to cal file
	%
	
	% Output file
	version = [vx '.' vy '.' vz];
	fsm_file = mms_construct_filename(sc, instr, mode, level, ...
	                                  'TStart',  tstart,      ...
	                                  'Version', version,     ...
	                                  'OptDesc', optdesc);

	% Find the latest z-version
	%   - Look in both DROPBOX and DATA_PATH
	vz = mms_latest_zversion(dropbox_root, fsm_file, 'RootDir', data_path_root);

	% Reform the file name
	version  = [vx '.' vy '.' sprintf('%0i', vz)];
	fsm_file = mms_construct_filename(sc, instr, mode, level, ...
	                                  'TStart',  tstart,      ...
	                                  'Version', version,     ...
	                                  'OptDesc', optdesc);
	fsm_file = fullfile(dropbox_root, fsm_file);

	% Notify where file is located
	mrfprintf('logtext', ['Creating FSM file at "' fsm_file '".']);

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
	[~, logical_file_id, ext] = fileparts(fsm_file);
	logical_file_id = [logical_file_id ext];
	
	%   - Instrument Type (1+)
	%           Electric Fields (space)
	%           Magnetic Fields (space)
	%           Particles (space)
	%           Plasma and Solar Wind
	%           Spacecraft Potential Control
	global_attrs = struct( 'Data_type',                  data_type, ...
	                       'Data_version',               version, ...
	                       'Descriptor',                 instr, ...
	                       'Discipline',                 'Space Physics>Magnetospheric Science', ...
	                       'Generation_date',            datestr(now(), 'yyyymmdd'), ...
	                       'Instrument_type',            'Magnetic Fields (space)', ...
	                       'Logical_file_id',            logical_file_id, ...
	                       'Logical_source',             logical_source, ...
	                       'Logical_source_description', ' ', ...
	                       'Mission_group',              'MMS', ...
	                       'PI_affiliation',             'SwRI, UNH, LPP, UCLA, IWF', ...
	                       'PI_name',                    'R.B. Torbert, M.R. Argall, O. LeContel, L. Mirioni, C. Russell, R. Strangeway, W. Magnes, D. Fischer', ...
	                       'Project',                    'STP>Solar Terrestrial Physics', ...
	                       'Source_name',                ['MMS' sc(4) '>MMS Satellite Number ' sc(4)], ...
	                       'TEXT',                       ['The merged magnetic field ' ...
	        'dataset is a combination of the FGM and SCM magnetometer data. ' ...
	        'Instrument papers for FGM and SCM, as well as their data products ' ...
	        'guides, can be found at the following links: ' ...
	        'http://dx.doi.org/10.1007/s11214-014-0057-3, ' ...
	        'http://dx.doi.org/10.1007/s11214-014-0096-9, ' ...
	        'https://lasp.colorado.edu/mms/sdc/public/datasets/fields/'], ...
	                       'HTTP_LINK',                  { {'http://mms-fields.unh.edu/' ...
	                                                        'http://mms.gsfc.nasa.gov/index.html'} }, ...
	                       'LINK_TEXT',                  { {'UNH FIELDS Home Page', ...
	                                                        'NASA MMS Home'} }, ...
	                       'MODS',                       { mods }, ...
	                       'Acknowledgements',           ' ', ...
	                       'Generated_by',               'University of New Hampshire', ...
	                       'Parents',                    { parents }, ...
	                       'Skeleton_version',           ' ', ...
	                       'Rules_of_use',               ' ', ...
	                       'Time_resolution',            ' '  ...
	                     );

%------------------------------------------------------
% Variable Names                                      |
%------------------------------------------------------
	% Variable naming convention
	%   scId_instrumentId_paramName[_coordSys][_paramQualifier][_subModeLevel][_mode][_level]
	prefix  = [sc '_' instr '_'];
	suffix  = ['_' mode '_' level];
	
	t_vname      = 'Epoch';
	b_vname      = [prefix 'b' '_mag'    suffix];
	b_omb_vname  = [prefix 'b' '_omb'    suffix];
	b_labl_vname = [prefix 'b' '_labels' suffix];

%------------------------------------------------------
% Metadata                                             |
%------------------------------------------------------
	% Labels
	b_labl = {'Bx', 'By', 'Bz'};

%------------------------------------------------------
% Variable Data                                       |
%------------------------------------------------------

	% Variables
	var_list = { t_vname,     data.tt2000, ...
	             b_vname,     data.b,      ...
	             b_omb_vname, data.b_omb,  ...
	             b_labl_vname, b_labl      ...
	           };

	recbound = { t_vname,    ...
	             b_vname,    ...
	             b_omb_vname ...
	           };

	% Variable data types
	vardatatypes = { t_vname,      'cdf_time_tt2000', ...
	                 b_vname,      'cdf_float',       ...
	                 b_omb_vname,  'cdf_float',       ...
	                 b_labl_vname, 'cdf_char'         ...
	               };
	
	% Variable compression
	varcompress = { b_vname,      'gzip.6', ...
	                b_omb_vname,  'gzip.6', ...
	              };

%------------------------------------------------------
% Variable Attributes                                 |
%------------------------------------------------------
	catdesc       = {};
	depend_0      = {};
	display_type  = {};
	fieldnam      = {};
	fillval       = {};
	format        = {};
	lablaxis      = {};
	labl_ptr_1    = {};
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
	
	% |B|
	catdesc       = [ catdesc,       b_vname, 'Magnetic field magnitude, computed from the vector in DMPA coordinates.' ];
	depend_0      = [ depend_0,      b_vname, t_vname ];
	display_type  = [ display_type,  b_vname, 'time_series'];
	fieldnam      = [ fieldnam,      b_vname, '|B|' ];
	fillval       = [ fillval,       b_vname, single(-1e31) ];
	format        = [ format,        b_vname, 'f11.4' ];
	si_conversion = [ si_conversion, b_vname, '1e-9>T' ];
	units         = [ units,         b_vname, 'nT' ];
	validmin      = [ validmin,      b_vname, single(0.0) ];
	validmax      = [ validmax,      b_vname, single(1e5) ];
	var_type      = [ var_type,      b_vname, 'data' ];
	
	% B_OMB
	catdesc       = [ catdesc,       b_omb_vname, 'Three components of the magnetic field in OMB coordinates.' ];
	depend_0      = [ depend_0,      b_omb_vname, t_vname ];
	display_type  = [ display_type,  b_omb_vname, 'time_series'];
	fieldnam      = [ fieldnam,      b_omb_vname, 'B' ];
	fillval       = [ fillval,       b_omb_vname, single(-1e31) ];
	format        = [ format,        b_omb_vname, 'f11.4' ];
	labl_ptr_1    = [ labl_ptr_1,    b_omb_vname, b_labl_vname ];
	si_conversion = [ si_conversion, b_omb_vname, '1e-9>T' ];
	units         = [ units,         b_omb_vname, 'nT' ];
	validmin      = [ validmin,      b_omb_vname, single(0.0) ];
	validmax      = [ validmax,      b_omb_vname, single(1e5) ];
	var_type      = [ var_type,      b_omb_vname, 'data' ];
	
	% B_LABL
	catdesc       = [ catdesc,  b_labl_vname, 'Axis labels for the magnetic field.' ];
	fieldnam      = [ fieldnam, b_labl_vname, 'Magnetic field labels' ];
	format        = [ format,   b_labl_vname, 'A2' ];
	var_type      = [ var_type, b_labl_vname, 'metadata' ];
	
%-------------------------------------------------------------------------
% Collect Variable Attributes ////////////////////////////////////////////
%-------------------------------------------------------------------------
	% Collect into a structure
	var_attrs = struct( 'CATDESC',       { catdesc },       ...
	                    'DEPEND_0',      { depend_0 },      ...
	                    'DISPLAY_TYPE',  { display_type },  ...
	                    'FIELDNAM',      { fieldnam },      ...
	                    'FILLVAL',       { fillval },       ...
	                    'FORMAT',        { format },        ...
	                    'LABLAXIS',      { lablaxis },      ...
	                    'LABL_PTR_1',    { labl_ptr_1 },    ...
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
	spdfcdfwrite( fsm_file, ...
	              var_list, ...
	              'GlobalAttributes',   global_attrs, ...
	              'RecordBound',        recbound,     ...
	              'VariableAttributes', var_attrs,    ...
	              'VarDatatypes',       vardatatypes, ...
	              'VarCompress',        varcompress   ...
	            );
	
	% If the file name is not output, print location to command window.
	if nargout == 0
		clear fsm_file
		disp(['File written to: "' fsm_file '".']);
	end
end