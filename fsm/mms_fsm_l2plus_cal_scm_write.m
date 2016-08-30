%
% Name
%   mms_fsm_l2plus_cal_scm_write
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
function fsm_cal_file = mms_fsm_l2plus_cal_scm_write(sc, mode, tstart, scm_bkgd, varargin)

	% Global variables
	%   - See mms_fsm_init.m
	global cal_path_root data_path_root dropbox_root hk_root log_path_root unh_data_root

%------------------------------------%
% Check Inputs                       %
%------------------------------------%

	% Defaults
	optdesc  = 'cal-scm';
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
	assert( isa(scm_bkgd.t,     'int64'),  'scm_bkgd.t should be int64.' );
	assert( isa(scm_bkgd.f,     'single'), 'scm_bkgd.f should be single.');
	assert( isa(scm_bkgd.amp,   'single'), 'scm_bkgd.amp should be single.');
	assert( isa(scm_bkgd.phase, 'single'), 'scm_bkgd.phase should be single.');
	assert( isa(scm_bkgd.psd,   'single'), 'scm_bkgd.psd should be single.');
	assert( isa(scm_bkgd.flag,  'uint8'),  'scm_bkgd.flag should be uint8.');
	
	% convert the start time to yyyy-mm-dd
%	tstart = MrTimeParser(tstart, '%Y-%M-%dT%H:%m:%S', '%Y%M%d');

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
	mrfprintf('logtext', ['Creating FSM scm cal file at "' fsm_cal_file '".']);

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
	                       'Logical_source_description', 'Level 2-Plus Fluxgate-Searchcoil Merged Magnetic Field Data', ...
	                       'Mission_group',              'MMS', ...
	                       'PI_affiliation',             'UNH, LPP, IWF, UCLA, UCLA, IWF, UNH', ...
	                       'PI_name',                    ['R.B. Torbert, O. LeContel, W. Magnes, C.T. Russell, R. Strangeway, ', ...
	                                                      'D. Fischer, M.R. Argall.'], ...
	                       'Project',                    'STP>Solar Terrestrial Physics', ...
	                       'Source_name',                source_name, ...
	                       'TEXT',                       ['SCM calibration data used while created the FSM data product. ' ...
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
	scm_amp_vname    = [prefix 'amp'        '_omb_scm' suffix];
	scm_phase_vname  = [prefix 'phase'      '_omb_scm' suffix];
	scm_psd_vname    = [prefix 'psd'        '_omb_scm' suffix];
	scm_gain_vname   = [prefix 'gain'       '_omb_scm' suffix];
	scm_offset_vname = [prefix 'phaseshift' '_omb_scm' suffix];
	scm_psdrat_vname = [prefix 'psdrat'     '_omb_scm' suffix];
	flag_vname       = [prefix 'flag'                  suffix];
	flag_hist_vname  = [prefix 'flag'       '_hist'    suffix];
	
	% SCM Histogram
	scm_hamp_vname    = [prefix 'amp'        '_hist_scm'  suffix];
	scm_hphase_vname  = [prefix 'phase'      '_hist_scm'  suffix];
	scm_hpsd_vname    = [prefix 'psd'        '_hist_scm'  suffix];
	scm_hgain_vname   = [prefix 'gain'       '_hist_scm'  suffix];
	scm_hoffset_vname = [prefix 'phaseshift' '_hist_scm'  suffix];
	scm_hpsdrat_vname = [prefix 'psdrat'     '_hist_scm'  suffix];
	
	% Histogram bins
	scm_bamp_vname    = [prefix 'amp'        '_bins_scm'  suffix];
	scm_bphase_vname  = [prefix 'phase'      '_bins_scm'  suffix];
	scm_bpsd_vname    = [prefix 'psd'        '_bins_scm'  suffix];
	scm_bgain_vname   = [prefix 'gain'       '_bins_scm'  suffix];
	scm_boffset_vname = [prefix 'phaseshift' '_bins_scm'  suffix];
	scm_bpsdrat_vname = [prefix 'psdrat'     '_bins_scm'  suffix];
	
	% Noise floor
	scm_famp_vname    = [prefix 'amp'        '_floor_scm' suffix];
	scm_fphase_vname  = [prefix 'phase'      '_floor_scm' suffix];
	scm_fpsd_vname    = [prefix 'psd'        '_floor_scm' suffix];
	scm_fgain_vname   = [prefix 'gain'       '_floor_scm' suffix];
	scm_foffset_vname = [prefix 'phaseshift' '_floor_scm' suffix];
	scm_fpsdrat_vname = [prefix 'psdrat'     '_floor_scm' suffix];
	
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
	len       = length(scm_bkgd.amp_bins);
	powten    = floor( log10( len ) ) + 1;
	fmt       = sprintf('%%0%ii', powten);
	bamp_labl = strcat( {'amp'}, num2str( [1:len]', fmt) )';
	
	% Component label
	c_labl = {'X', 'Y', 'Z'};
	
	% Frequency Labels
	len    = length(scm_bkgd.f);
	powten = floor( log10( len ) ) + 1;
	fmt    = sprintf('%%0%ii', powten);
	f_labl = strcat( {'freq'}, num2str( [1:len]', fmt) )';
	
	% Hist Flag Label
	len        = length(scm_bkgd.hist_flag);
	powten     = floor( log10( len ) ) + 1;
	fmt        = sprintf('%%0%ii', powten);
	hflag_labl = strcat( {'histflag'}, num2str( [1:len]', fmt) )';
	
	% Phase bins Label
	len         = length(scm_bkgd.phase_bins);
	powten      = floor( log10( len ) ) + 1;
	fmt         = sprintf('%%0%ii', powten);
	bphase_labl = strcat( {'phase'}, num2str( [1:len]', fmt) )';
	
	% PSD bins Label
	len       = length(scm_bkgd.psd_bins);
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
	var_list = { t_vname,            scm_bkgd.t,     ...
	             f_vname,            scm_bkgd.f,     ...
	             c_index_vname,      c_index,            ...
	             scm_amp_vname,      permute( scm_bkgd.amp,         [3,1,2] ), ...
	             scm_phase_vname,    permute( scm_bkgd.phase,       [3,1,2] ), ...
	             scm_psd_vname,      permute( scm_bkgd.psd,         [3,1,2] ), ...
	             flag_vname,         scm_bkgd.flag,      ...
	             bamp_labl_vname,    bamp_labl,          ...
	             bphase_labl_vname,  bphase_labl,        ...
	             bpsd_labl_vname,    bpsd_labl,          ...
	             c_labl_vname,       c_labl,             ...
	             f_labl_vname,       f_labl,             ...
	             hflag_labl_vname,   hflag_labl,         ...
	             flag_hist_vname,    scm_bkgd.hist_flag, ...
	             scm_hamp_vname,     permute( scm_bkgd.amp_hist,         [3,2,1,4] ), ...
	             scm_hphase_vname,   permute( scm_bkgd.phase_hist,       [3,2,1,4] ), ...
	             scm_hpsd_vname,     permute( scm_bkgd.psd_hist,         [3,2,1,4] ), ...
	             scm_famp_vname,     permute( scm_bkgd.amp_floor,         [2,1,3] ),  ...
	             scm_fphase_vname,   permute( scm_bkgd.phase_floor,       [2,1,3] ),  ...
	             scm_fpsd_vname,     permute( scm_bkgd.psd_floor,         [2,1,3] ),  ...
	             scm_bamp_vname,     scm_bkgd.amp_bins,         ...
	             scm_bphase_vname,   scm_bkgd.phase_bins,       ...
	             scm_bpsd_vname,     scm_bkgd.psd_bins          ...
	           };
%	             scm_gain_vname,     permute( scm_bkgd.gain,        [3,1,2] ), ...
%	             scm_offset_vname,   permute( scm_bkgd.phase_shift, [3,1,2] ), ...
%	             scm_psdrat_vname,   permute( scm_bkgd.psd_rat,     [3,1,2] ), ...
%	             scm_hgain_vname,    permute( scm_bkgd.gain_hist,        [3,2,1,4] ), ...
%	             scm_hoffset_vname,  permute( scm_bkgd.phase_shift_hist, [3,2,1,4] ), ...
%	             scm_hpsdrat_vname,  permute( scm_bkgd.psd_rat_hist,     [3,2,1,4] ), ...
%	             scm_fgain_vname,    permute( scm_bkgd.gain_floor,        [2,1,3] ),  ...
%	             scm_foffset_vname,  permute( scm_bkgd.phase_shift_floor, [2,1,3] ),  ...
%	             scm_fpsdrat_vname,  permute( scm_bkgd.psd_rat_floor,     [2,1,3] ),  ...
%	             scm_bgain_vname,    scm_bkgd.gain_bins,        ...
%	             scm_boffset_vname,  scm_bkgd.phase_shift_bins, ...
%	             scm_bpsdrat_vname,  scm_bkgd.psd_rat_bins      ...
	clear scm_bkgd

	% Record Variance
	recbound = { t_vname,           ...
	             scm_amp_vname,     ...
	             scm_phase_vname,   ...
	             scm_psd_vname,     ...
	             flag_vname,        ...
	             scm_hamp_vname,    ...
	             scm_hphase_vname,  ...
	             scm_hpsd_vname,    ...
	             scm_famp_vname,    ...
	             scm_fphase_vname,  ...
	             scm_fpsd_vname     ...
	           };
%	             scm_gain_vname,    ...
%	             scm_offset_vname,  ...
%	             scm_psdrat_vname,  ...
%	             scm_hgain_vname,   ...
%	             scm_hoffset_vname, ...
%	             scm_hpsdrat_vname, ...
%	             scm_fgain_vname,   ...
%	             scm_foffset_vname, ...
%	             scm_fpsdrat_vname  ...

	% Data types
	vardatatypes = { t_vname,           'cdf_time_tt2000', ...
	                 f_vname,           'cdf_float', ...
	                 c_index_vname,     'cdf_uint1', ...
	                 scm_amp_vname,     'cdf_float', ...
	                 scm_phase_vname,   'cdf_float', ...
	                 scm_psd_vname,     'cdf_float', ...
	                 flag_vname,        'cdf_uint1', ...
	                 bamp_labl_vname,   'cdf_char',  ...
	                 bphase_labl_vname, 'cdf_char',  ...
	                 bpsd_labl_vname,   'cdf_char',  ...
	                 c_labl_vname,      'cdf_char',  ...
	                 f_labl_vname,      'cdf_char',  ...
	                 hflag_labl_vname,  'cdf_char',  ...
	                 flag_hist_vname,   'cdf_uint1', ...
	                 scm_hamp_vname,    'cdf_uint4', ...
	                 scm_hphase_vname,  'cdf_uint4', ...
	                 scm_hpsd_vname,    'cdf_uint4', ...
	                 scm_famp_vname,    'cdf_float', ...
	                 scm_fphase_vname,  'cdf_float', ...
	                 scm_fpsd_vname,    'cdf_float', ...
	                 scm_bamp_vname,    'cdf_float', ...
	                 scm_bphase_vname,  'cdf_float', ...
	                 scm_bpsd_vname     'cdf_float'  ...
	               };
%	                 scm_gain_vname,    'cdf_float', ...
%	                 scm_offset_vname,  'cdf_float', ...
%	                 scm_psdrat_vname,  'cdf_float', ...
%	                 scm_hgain_vname,   'cdf_uint4', ...
%	                 scm_hoffset_vname, 'cdf_uint4', ...
%	                 scm_hpsdrat_vname, 'cdf_uint4', ...
%	                 scm_fgain_vname,   'cdf_float', ...
%	                 scm_foffset_vname, 'cdf_float', ...
%	                 scm_fpsdrat_vname, 'cdf_float', ...
%	                 scm_bgain_vname,   'cdf_float', ...
%	                 scm_boffset_vname, 'cdf_float', ...
%	                 scm_bpsdrat_vname, 'cdf_float'  ...
	
	% CDF Compression Level
	varcompress = { scm_amp_vname,     'gzip.6', ...
	                scm_phase_vname,   'gzip.6', ...
	                scm_psd_vname,     'gzip.6', ...
	                flag_vname,        'gzip.6', ...
	                flag_hist_vname,   'gzip.6', ...
	                scm_hamp_vname,    'gzip.6', ...
	                scm_hphase_vname,  'gzip.6', ...
	                scm_hpsd_vname,    'gzip.6', ...
	                scm_famp_vname,    'gzip.6', ...
	                scm_fphase_vname,  'gzip.6', ...
	                scm_fpsd_vname,    'gzip.6', ...
	                scm_bamp_vname,    'gzip.6', ...
	                scm_bphase_vname,  'gzip.6', ...
	                scm_bpsd_vname     'gzip.6'  ...
	              };
%	                scm_gain_vname,    'gzip.6', ...
%	                scm_offset_vname,  'gzip.6', ...
%	                scm_psdrat_vname,  'gzip.6', ...
%	                scm_hgain_vname,   'gzip.6', ...
%	                scm_hoffset_vname, 'gzip.6', ...
%	                scm_hpsdrat_vname, 'gzip.6', ...
%	                scm_fgain_vname,   'gzip.6', ...
%	                scm_foffset_vname, 'gzip.6', ...
%	                scm_fpsdrat_vname, 'gzip.6', ...
%	                scm_bgain_vname    'gzip.6', ...
%	                scm_boffset_vname  'gzip.6', ...
%	                scm_bpsdrat_vname  'gzip.6'  ...
	
	singleton = { scm_amp_vname,    ...
	              scm_phase_vname,  ...
	              scm_psd_vname,    ...
	              scm_hamp_vname,   ...
	              scm_hphase_vname, ...
	              scm_hpsd_vname,   ...
	              scm_famp_vname,   ...
	              scm_fphase_vname, ...
	              scm_fpsd_vname    ...
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
	catdesc       = [ catdesc,       scm_amp_vname, 'Wave amplitude of the SCM magnetic field.' ];
	depend_0      = [ depend_0,      scm_amp_vname, t_vname ];
	depend_1      = [ depend_1,      scm_amp_vname, f_vname ];
	depend_2      = [ depend_2,      scm_amp_vname, c_index_vname ];
	display_type  = [ display_type,  scm_amp_vname, 'spectrogram'];
	fieldnam      = [ fieldnam,      scm_amp_vname, 'Amplitude' ];
	fillval       = [ fillval,       scm_amp_vname, single(-1e31) ];
	format        = [ format,        scm_amp_vname, 'f11.4' ];
	labl_ptr_1    = [ labl_ptr_1,    scm_amp_vname, f_labl_vname ];
	labl_ptr_2    = [ labl_ptr_2,    scm_amp_vname, c_labl_vname ];
	si_conversion = [ si_conversion, scm_amp_vname, '1e-9>T' ];
	units         = [ units,         scm_amp_vname, 'nT' ];
	validmin      = [ validmin,      scm_amp_vname, single(0.0) ];
	validmax      = [ validmax,      scm_amp_vname, single(1e5) ];
	var_type      = [ var_type,      scm_amp_vname, 'data' ];
	
	% PHASE
	catdesc       = [ catdesc,       scm_phase_vname, 'Wave phase of the SCM magnetic field.' ];
	depend_0      = [ depend_0,      scm_phase_vname, t_vname ];
	depend_1      = [ depend_1,      scm_phase_vname, f_vname ];
	depend_2      = [ depend_2,      scm_phase_vname, c_index_vname ];
	display_type  = [ display_type,  scm_phase_vname, 'spectrogram'];
	fieldnam      = [ fieldnam,      scm_phase_vname, 'Phase' ];
	fillval       = [ fillval,       scm_phase_vname, single(-1e31) ];
	format        = [ format,        scm_phase_vname, 'f11.4' ];
	labl_ptr_1    = [ labl_ptr_1,    scm_phase_vname, f_labl_vname ];
	labl_ptr_2    = [ labl_ptr_2,    scm_phase_vname, c_labl_vname ];
	si_conversion = [ si_conversion, scm_phase_vname, '1e-9>T' ];
	units         = [ units,         scm_phase_vname, 'degrees' ];
	validmin      = [ validmin,      scm_phase_vname, single(0.0) ];
	validmax      = [ validmax,      scm_phase_vname, single(360.0) ];
	var_type      = [ var_type,      scm_phase_vname, 'data' ];
	
	% PSD
	catdesc       = [ catdesc,       scm_psd_vname, 'Power spectral density of the SCM magnetic field.' ];
	depend_0      = [ depend_0,      scm_psd_vname, t_vname ];
	depend_1      = [ depend_1,      scm_psd_vname, f_vname ];
	depend_2      = [ depend_2,      scm_psd_vname, c_index_vname ];
	display_type  = [ display_type,  scm_psd_vname, 'spectrogram'];
	fieldnam      = [ fieldnam,      scm_psd_vname, 'Phase' ];
	fillval       = [ fillval,       scm_psd_vname, single(-1e31) ];
	format        = [ format,        scm_psd_vname, 'f11.4' ];
	labl_ptr_1    = [ labl_ptr_1,    scm_psd_vname, f_labl_vname ];
	labl_ptr_2    = [ labl_ptr_2,    scm_psd_vname, c_labl_vname ];
	si_conversion = [ si_conversion, scm_psd_vname, '1e-9>T' ];
	units         = [ units,         scm_psd_vname, 'nT^2/Hz' ];
	validmin      = [ validmin,      scm_psd_vname, single(0.0) ];
	validmax      = [ validmax,      scm_psd_vname, single(1e4) ];
	var_type      = [ var_type,      scm_psd_vname, 'data' ];
	
	% GAIN
%	catdesc       = [ catdesc,       scm_gain_vname, 'Ratio of DFG to SCM wave amplitude.' ];
%	depend_0      = [ depend_0,      scm_gain_vname, t_vname ];
%	depend_1      = [ depend_1,      scm_gain_vname, f_vname ];
%	display_type  = [ display_type,  scm_gain_vname, 'spectrogram'];
%	fieldnam      = [ fieldnam,      scm_gain_vname, 'Gain' ];
%	fillval       = [ fillval,       scm_gain_vname, single(-1e31) ];
%	format        = [ format,        scm_gain_vname, 'f11.4' ];
%	labl_ptr_1    = [ labl_ptr_1,    scm_gain_vname, c_labl_vname ];
%	si_conversion = [ si_conversion, scm_gain_vname, ' ' ];
%	units         = [ units,         scm_gain_vname, ' ' ];
%	validmin      = [ validmin,      scm_gain_vname, single(0.0) ];
%	validmax      = [ validmax,      scm_gain_vname, single(1e4) ];
%	var_type      = [ var_type,      scm_gain_vname, 'data' ];
	
	% PHASE OFFSET
%	catdesc       = [ catdesc,       scm_offset_vname, 'Difference of DFG and SCM phase.' ];
%	depend_0      = [ depend_0,      scm_offset_vname, t_vname ];
%	depend_1      = [ depend_1,      scm_offset_vname, f_vname ];
%	display_type  = [ display_type,  scm_offset_vname, 'spectrogram'];
%	fieldnam      = [ fieldnam,      scm_offset_vname, 'Phase Diff.' ];
%	fillval       = [ fillval,       scm_offset_vname, single(-1e31) ];
%	format        = [ format,        scm_offset_vname, 'f11.4' ];
%	labl_ptr_1    = [ labl_ptr_1,    scm_offset_vname, c_labl_vname ];
%	si_conversion = [ si_conversion, scm_offset_vname, '1e0>degrees' ];
%	units         = [ units,         scm_offset_vname, 'degrees' ];
%	validmin      = [ validmin,      scm_offset_vname, single(-360.0) ];
%	validmax      = [ validmax,      scm_offset_vname, single(360.0) ];
%	var_type      = [ var_type,      scm_offset_vname, 'data' ];
	
	% PSD RATIO
%	catdesc       = [ catdesc,       scm_psdrat_vname, 'Ratio of DFG to SCM power spectral density.' ];
%	depend_0      = [ depend_0,      scm_psdrat_vname, t_vname ];
%	depend_1      = [ depend_1,      scm_psdrat_vname, f_vname ];
%	display_type  = [ display_type,  scm_psdrat_vname, 'spectrogram'];
%	fieldnam      = [ fieldnam,      scm_psdrat_vname, 'PSD' ];
%	fillval       = [ fillval,       scm_psdrat_vname, single(-1e31) ];
%	format        = [ format,        scm_psdrat_vname, 'f11.4' ];
%	labl_ptr_1    = [ labl_ptr_1,    scm_psdrat_vname, c_labl_vname ];
%	si_conversion = [ si_conversion, scm_psdrat_vname, '1e-18>T^2/Hz' ];
%	units         = [ units,         scm_psdrat_vname, 'nT^2/Hz' ];
%	validmin      = [ validmin,      scm_psdrat_vname, single(1e4) ];
%	validmax      = [ validmax,      scm_psdrat_vname, single(0) ];
%	var_type      = [ var_type,      scm_psdrat_vname, 'data' ];
	
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
% SCM HISTOGRAM INFO /////////////////////////////////////////////////////
%-------------------------------------------------------------------------
	
	% AMPLITUDE
	catdesc       = [ catdesc,       scm_hamp_vname, 'SCM wave amplitude, Histogrammed in time at each frequency.' ];
	depend_0      = [ depend_0,      scm_hamp_vname, f_vname ];
	depend_1      = [ depend_1,      scm_hamp_vname, scm_bamp_vname ];
	depend_2      = [ depend_2,      scm_hamp_vname, flag_hist_vname ];
	depend_3      = [ depend_3,      scm_hamp_vname, c_index_vname ];
	display_type  = [ display_type,  scm_hamp_vname, 'spectrogram'];
	fieldnam      = [ fieldnam,      scm_hamp_vname, 'Histogrammed wave amplitude (SCM)' ];
	fillval       = [ fillval,       scm_hamp_vname, uint32(4294967295) ];
	format        = [ format,        scm_hamp_vname, 'I10' ];
	labl_ptr_1    = [ labl_ptr_1,    scm_hamp_vname, bamp_labl_vname ];
	labl_ptr_2    = [ labl_ptr_2,    scm_hamp_vname, hflag_labl_vname ];
	labl_ptr_3    = [ labl_ptr_3,    scm_hamp_vname, c_labl_vname ];
	si_conversion = [ si_conversion, scm_hamp_vname, '1e-9>T' ];
	units         = [ units,         scm_hamp_vname, 'degrees' ];
	validmin      = [ validmin,      scm_hamp_vname, uint32(0) ];
	validmax      = [ validmax,      scm_hamp_vname, uint32(4294967294) ];
	var_type      = [ var_type,      scm_hamp_vname, 'data' ];
	
	% PHASE
	catdesc       = [ catdesc,       scm_hphase_vname, 'SCM wave phase, Histogrammed in time at each frequency.' ];
	depend_0      = [ depend_0,      scm_hphase_vname, f_vname ];
	depend_1      = [ depend_1,      scm_hphase_vname, scm_bphase_vname ];
	depend_2      = [ depend_2,      scm_hphase_vname, flag_hist_vname ];
	depend_3      = [ depend_3,      scm_hphase_vname, c_index_vname ];
	display_type  = [ display_type,  scm_hphase_vname, 'spectrogram'];
	fieldnam      = [ fieldnam,      scm_hphase_vname, 'Histogrammed wave phase (SCM)' ];
	fillval       = [ fillval,       scm_hphase_vname, uint32(4294967295) ];
	format        = [ format,        scm_hphase_vname, 'I10' ];
	labl_ptr_1    = [ labl_ptr_1,    scm_hphase_vname, bphase_labl_vname ];
	labl_ptr_2    = [ labl_ptr_2,    scm_hphase_vname, hflag_labl_vname ];
	labl_ptr_3    = [ labl_ptr_3,    scm_hphase_vname, c_labl_vname ];
	si_conversion = [ si_conversion, scm_hphase_vname, '1e-9>T' ];
	units         = [ units,         scm_hphase_vname, 'degrees' ];
	validmin      = [ validmin,      scm_hphase_vname, uint32(0) ];
	validmax      = [ validmax,      scm_hphase_vname, uint32(4294967294) ];
	var_type      = [ var_type,      scm_hphase_vname, 'data' ];
	
	% POWER SPECTRAL DENSITY
	catdesc       = [ catdesc,       scm_hpsd_vname, 'SCM wave power, Histogrammed in time at each frequency.' ];
	depend_0      = [ depend_0,      scm_hpsd_vname, f_vname ];
	depend_1      = [ depend_1,      scm_hpsd_vname, scm_bpsd_vname ];
	depend_2      = [ depend_2,      scm_hpsd_vname, flag_hist_vname ];
	depend_3      = [ depend_3,      scm_hpsd_vname, c_index_vname ];
	display_type  = [ display_type,  scm_hpsd_vname, 'spectrogram'];
	fieldnam      = [ fieldnam,      scm_hpsd_vname, 'Histogrammed wave power (SCM)' ];
	fillval       = [ fillval,       scm_hpsd_vname, uint32(4294967295) ];
	format        = [ format,        scm_hpsd_vname, 'I10' ];
	labl_ptr_1    = [ labl_ptr_1,    scm_hpsd_vname, bpsd_labl_vname ];
	labl_ptr_2    = [ labl_ptr_2,    scm_hpsd_vname, hflag_labl_vname ];
	labl_ptr_3    = [ labl_ptr_3,    scm_hpsd_vname, c_labl_vname ];
	si_conversion = [ si_conversion, scm_hpsd_vname, '1e-9>T' ];
	units         = [ units,         scm_hpsd_vname, 'degrees' ];
	validmin      = [ validmin,      scm_hpsd_vname, uint32(0) ];
	validmax      = [ validmax,      scm_hpsd_vname, uint32(4294967294) ];
	var_type      = [ var_type,      scm_hpsd_vname, 'data' ];
	
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
%	catdesc       = [ catdesc,       scm_hgain_vname, 'Ratio of DFG to SCM wave amplitude, histogrammed in time at each frequency.' ];
%	depend_0      = [ depend_0,      scm_hgain_vname, f_vname ];
%	depend_1      = [ depend_1,      scm_hgain_vname, scm_bgain_vname ];
%	depend_2      = [ depend_2,      scm_hgain_vname, flag_hist_vname ];
%	depend_3      = [ depend_3,      scm_hgain_vname, c_index_vname ];
%	display_type  = [ display_type,  scm_hgain_vname, 'spectrogram'];
%	fieldnam      = [ fieldnam,      scm_hgain_vname, 'Histogrammed gain' ];
%	fillval       = [ fillval,       scm_hgain_vname, uint32(65535) ];
%	format        = [ format,        scm_hgain_vname, 'I5' ];
%	labl_ptr_1    = [ labl_ptr_1,    scm_hgain_vname, c_labl_vname ];
%	si_conversion = [ si_conversion, scm_hgain_vname, ' ' ];
%	units         = [ units,         scm_hgain_vname, ' ' ];
%	validmin      = [ validmin,      scm_hgain_vname, uint32(0) ];
%	validmax      = [ validmax,      scm_hgain_vname, uint32(65533) ];
%	var_type      = [ var_type,      scm_hgain_vname, 'data' ];
	
	% PHASE SHIFT
%	catdesc       = [ catdesc,       scm_hoffset_vname, 'Difference between DFG and SCM wave phase, histogrammed in time at each frequency.' ];
%	depend_0      = [ depend_0,      scm_hoffset_vname, f_vname ];
%	depend_1      = [ depend_1,      scm_hoffset_vname, scm_boffset_vname ];
%	depend_2      = [ depend_2,      scm_hoffset_vname, flag_hist_vname ];
%	depend_3      = [ depend_3,      scm_hoffset_vname, c_index_vname ];
%	display_type  = [ display_type,  scm_hoffset_vname, 'spectrogram'];
%	fieldnam      = [ fieldnam,      scm_hoffset_vname, 'Histogrammed phase shift' ];
%	fillval       = [ fillval,       scm_hoffset_vname, uint32(65535) ];
%	format        = [ format,        scm_hoffset_vname, 'I5' ];
%	labl_ptr_1    = [ labl_ptr_1,    scm_hoffset_vname, c_labl_vname ];
%	si_conversion = [ si_conversion, scm_hoffset_vname, '1e0>degrees' ];
%	units         = [ units,         scm_hoffset_vname, 'degrees' ];
%	validmin      = [ validmin,      scm_hoffset_vname, uint32(0) ];
%	validmax      = [ validmax,      scm_hoffset_vname, uint32(65533) ];
%	var_type      = [ var_type,      scm_hoffset_vname, 'data' ];
	
	% PSDRAT
%	catdesc       = [ catdesc,       scm_hpsdrat_vname, 'Ratio of DFG to SCM power spectral density, histogrammed in time at each frequency.' ];
%	depend_0      = [ depend_0,      scm_hpsdrat_vname, f_vname ];
%	depend_1      = [ depend_1,      scm_hpsdrat_vname, scm_bpsdrat_vname ];
%	depend_2      = [ depend_2,      scm_hpsdrat_vname, flag_hist_vname ];
%	depend_3      = [ depend_3,      scm_hpsdrat_vname, c_index_vname ];
%	display_type  = [ display_type,  scm_hpsdrat_vname, 'spectrogram'];
%	fieldnam      = [ fieldnam,      scm_hpsdrat_vname, 'Histogrammed phase shift' ];
%	fillval       = [ fillval,       scm_hpsdrat_vname, uint32(65535) ];
%	format        = [ format,        scm_hpsdrat_vname, 'I5' ];
%	labl_ptr_1    = [ labl_ptr_1,    scm_hpsdrat_vname, c_labl_vname ];
%	si_conversion = [ si_conversion, scm_hpsdrat_vname, '1e0>degrees' ];
%	units         = [ units,         scm_hpsdrat_vname, 'degrees' ];
%	validmin      = [ validmin,      scm_hpsdrat_vname, uint32(0) ];
%	validmax      = [ validmax,      scm_hpsdrat_vname, uint32(65533) ];
%	var_type      = [ var_type,      scm_hpsdrat_vname, 'data' ];
	
%-------------------------------------------------------------------------
% SCM NOISE FLOOR INFO ///////////////////////////////////////////////////
%-------------------------------------------------------------------------
	
	% FLOOR_AMP_SCM
	catdesc       = [ catdesc,       scm_famp_vname, ['Log base 10 of the most common wave amplitude, obtained by histogramming ' ...
	                                                  'wave amplitude at each frequency over the course of a day.'] ];
	depend_0      = [ depend_0,      scm_famp_vname, f_vname ];
	depend_1      = [ depend_1,      scm_famp_vname, flag_hist_vname ];
	depend_2      = [ depend_2,      scm_famp_vname, c_index_vname ];
	display_type  = [ display_type,  scm_famp_vname, 'time_series'];
	fieldnam      = [ fieldnam,      scm_famp_vname, 'Wave amplitude' ];
	fillval       = [ fillval,       scm_famp_vname, single(-1e31) ];
	format        = [ format,        scm_famp_vname, 'f11.4' ];
	labl_ptr_1    = [ labl_ptr_1,    scm_famp_vname, hflag_labl_vname ];
	labl_ptr_2    = [ labl_ptr_2,    scm_famp_vname, c_labl_vname ];
	si_conversion = [ si_conversion, scm_famp_vname, '1e-9>T' ];
	units         = [ units,         scm_famp_vname, 'degrees' ];
	validmin      = [ validmin,      scm_famp_vname, single(-10.0) ];
	validmax      = [ validmax,      scm_famp_vname, single(10.0) ];
	var_type      = [ var_type,      scm_famp_vname, 'data' ];
	
	% FLOOR_PHASE_SCM
	catdesc       = [ catdesc,       scm_fphase_vname, ['The most common wave phase, obtained by histogramming ' ...
	                                                    'wave phase at each frequency over the course of a day.'] ];
	depend_0      = [ depend_0,      scm_fphase_vname, f_vname ];
	depend_1      = [ depend_1,      scm_fphase_vname, flag_hist_vname ];
	depend_2      = [ depend_2,      scm_fphase_vname, c_index_vname ];
	display_type  = [ display_type,  scm_fphase_vname, 'time_series'];
	fieldnam      = [ fieldnam,      scm_fphase_vname, 'Wave phase' ];
	fillval       = [ fillval,       scm_fphase_vname, single(-1e31) ];
	format        = [ format,        scm_fphase_vname, 'f11.4' ];
	labl_ptr_1    = [ labl_ptr_1,    scm_fphase_vname, hflag_labl_vname ];
	labl_ptr_2    = [ labl_ptr_2,    scm_fphase_vname, c_labl_vname ];
	si_conversion = [ si_conversion, scm_fphase_vname, '1e0>degrees' ];
	units         = [ units,         scm_fphase_vname, 'degrees' ];
	validmin      = [ validmin,      scm_fphase_vname, single(-180.0) ];
	validmax      = [ validmax,      scm_fphase_vname, single(180.0) ];
	var_type      = [ var_type,      scm_fphase_vname, 'data' ];
	
	
	% FLOOR_PSD_SCM
	catdesc       = [ catdesc,       scm_fpsd_vname, ['Log base 10 of the most common wave power, obtained by histogramming ' ...
	                                                  'wave power at each frequency over the course of a day.'] ];
	depend_0      = [ depend_0,      scm_fpsd_vname, f_vname ];
	depend_1      = [ depend_1,      scm_fpsd_vname, flag_hist_vname ];
	depend_2      = [ depend_2,      scm_fpsd_vname, c_index_vname ];
	display_type  = [ display_type,  scm_fpsd_vname, 'time_series'];
	fieldnam      = [ fieldnam,      scm_fpsd_vname, 'Power spectral density' ];
	fillval       = [ fillval,       scm_fpsd_vname, single(-1e31) ];
	format        = [ format,        scm_fpsd_vname, 'f11.4' ];
	labl_ptr_1    = [ labl_ptr_1,    scm_fpsd_vname, hflag_labl_vname ];
	labl_ptr_2    = [ labl_ptr_2,    scm_fpsd_vname, c_labl_vname ];
	si_conversion = [ si_conversion, scm_fpsd_vname, '1e-18>T^2 Hz^-1' ];
	units         = [ units,         scm_fpsd_vname, 'nT^2 Hz^-1' ];
	validmin      = [ validmin,      scm_fpsd_vname, single(-15.0) ];
	validmax      = [ validmax,      scm_fpsd_vname, single(15.0) ];
	var_type      = [ var_type,      scm_fpsd_vname, 'data' ];
	
	% FLOOR GAIN
%	catdesc       = [ catdesc,       scm_fgain_vname, ['Log base 10 of the most common DFG to SCM gain factor, obtained by histogramming ' ...
%	                                                   'gain at each frequency over the course of a day.'] ];
%	depend_0      = [ depend_0,      scm_fgain_vname, f_vname ];
%	depend_1      = [ depend_1,      scm_fgain_vname, flag_hist_vname ];
%	depend_2      = [ depend_2,      scm_fgain_vname, c_index_vname ];
%	display_type  = [ display_type,  scm_fgain_vname, 'time_series'];
%	fieldnam      = [ fieldnam,      scm_fgain_vname, 'Gain' ];
%	fillval       = [ fillval,       scm_fgain_vname, single(-1e31) ];
%	format        = [ format,        scm_fgain_vname, 'f11.4' ];
%	labl_ptr_1    = [ labl_ptr_1,    scm_fgain_vname, c_labl_vname ];
%	si_conversion = [ si_conversion, scm_fgain_vname, ' ' ];
%	units         = [ units,         scm_fgain_vname, ' ' ];
%	validmin      = [ validmin,      scm_fgain_vname, single(-15.0) ];
%	validmax      = [ validmax,      scm_fgain_vname, single(15.0) ];
%	var_type      = [ var_type,      scm_fgain_vname, 'data' ];
	
	% FLOOR PHASE SHIFT
%	catdesc       = [ catdesc,       scm_foffset_vname, ['The most common DFG to SCM phase shift, obtained by histogramming ' ...
%	                                                     'phase shift at each frequency over the course of a day.'] ];
%	depend_0      = [ depend_0,      scm_foffset_vname, f_vname ];
%	depend_1      = [ depend_1,      scm_foffset_vname, flag_hist_vname ];
%	depend_2      = [ depend_2,      scm_foffset_vname, c_index_vname ];
%	display_type  = [ display_type,  scm_foffset_vname, 'time_series'];
%	fieldnam      = [ fieldnam,      scm_foffset_vname, 'Gain' ];
%	fillval       = [ fillval,       scm_foffset_vname, single(-1e31) ];
%	format        = [ format,        scm_foffset_vname, 'f11.4' ];
%	labl_ptr_1    = [ labl_ptr_1,    scm_foffset_vname, c_labl_vname ];
%	si_conversion = [ si_conversion, scm_foffset_vname, ' ' ];
%	units         = [ units,         scm_foffset_vname, ' ' ];
%	validmin      = [ validmin,      scm_foffset_vname, single(-15.0) ];
%	validmax      = [ validmax,      scm_foffset_vname, single(15.0) ];
%	var_type      = [ var_type,      scm_foffset_vname, 'data' ];
	
	% FLOOR PSD RATIO
%	catdesc       = [ catdesc,       scm_fpsdrat_vname, ['Log base 10 of the most common DFG to SCM PSD ratio, obtained by histogramming ' ...
%	                                                     'PSD ratio at each frequency over the course of a day.'] ];
%	depend_0      = [ depend_0,      scm_fpsdrat_vname, f_vname ];
%	depend_1      = [ depend_1,      scm_fpsdrat_vname, flag_hist_vname ];
%	depend_2      = [ depend_2,      scm_fpsdrat_vname, c_index_vname ];
%	display_type  = [ display_type,  scm_fpsdrat_vname, 'time_series'];
%	fieldnam      = [ fieldnam,      scm_fpsdrat_vname, 'Gain' ];
%	fillval       = [ fillval,       scm_fpsdrat_vname, single(-1e31) ];
%	format        = [ format,        scm_fpsdrat_vname, 'f11.4' ];
%	labl_ptr_1    = [ labl_ptr_1,    scm_fpsdrat_vname, c_labl_vname ];
%	si_conversion = [ si_conversion, scm_fpsdrat_vname, ' ' ];
%	units         = [ units,         scm_fpsdrat_vname, ' ' ];
%	validmin      = [ validmin,      scm_fpsdrat_vname, single(-15.0) ];
%	validmax      = [ validmax,      scm_fpsdrat_vname, single(15.0) ];
%	var_type      = [ var_type,      scm_fpsdrat_vname, 'data' ];
	
%-------------------------------------------------------------------------
% SCM HISTOGRAM BINS INFO ////////////////////////////////////////////////
%-------------------------------------------------------------------------
	
	% BINS_AMP_SCM
	catdesc       = [ catdesc,       scm_bamp_vname, 'Histogram bins for wave amplitude.'];
	fieldnam      = [ fieldnam,      scm_bamp_vname, 'Wave amplitude bins' ];
	fillval       = [ fillval,       scm_bamp_vname, single(-1e31) ];
	format        = [ format,        scm_bamp_vname, 'f11.4' ];
	lablaxis      = [ lablaxis,      scm_bamp_vname, 'Amp bins' ];
	si_conversion = [ si_conversion, scm_bamp_vname, '1e-9>T' ];
	units         = [ units,         scm_bamp_vname, 'nT' ];
	validmin      = [ validmin,      scm_bamp_vname, single(-10.0) ];
	validmax      = [ validmax,      scm_bamp_vname, single(10.0) ];
	var_type      = [ var_type,      scm_bamp_vname, 'support_data' ];
	
	% BINS_PHASE_SCM
	catdesc       = [ catdesc,       scm_bphase_vname, 'Histogram bins for wave phase.'];
	fieldnam      = [ fieldnam,      scm_bphase_vname, 'Wave amplitude bins' ];
	fillval       = [ fillval,       scm_bphase_vname, single(-1e31) ];
	format        = [ format,        scm_bphase_vname, 'f11.4' ];
	lablaxis      = [ lablaxis,      scm_bphase_vname, 'Phase bins' ];
	si_conversion = [ si_conversion, scm_bphase_vname, '1e0>degrees' ];
	units         = [ units,         scm_bphase_vname, 'degrees' ];
	validmin      = [ validmin,      scm_bphase_vname, single(-180.0) ];
	validmax      = [ validmax,      scm_bphase_vname, single(180.0) ];
	var_type      = [ var_type,      scm_bphase_vname, 'support_data' ];
	
	% BINS_PSD_SCM
	catdesc       = [ catdesc,       scm_bpsd_vname, 'Histogram bins for wave amplitude.' ];
	fieldnam      = [ fieldnam,      scm_bpsd_vname, 'Wave amplitude bins' ];
	fillval       = [ fillval,       scm_bpsd_vname, single(-1e31) ];
	format        = [ format,        scm_bpsd_vname, 'f11.4' ];
	lablaxis      = [ lablaxis,      scm_bpsd_vname, 'PSD bins' ];
	si_conversion = [ si_conversion, scm_bpsd_vname, '1e-18>T^2 Hz^-1' ];
	units         = [ units,         scm_bpsd_vname, 'nT^2 Hz^-1' ];
	validmin      = [ validmin,      scm_bpsd_vname, single(-10.0) ];
	validmax      = [ validmax,      scm_bpsd_vname, single(10.0) ];
	var_type      = [ var_type,      scm_bpsd_vname, 'support_data' ];
	
	% GAIN BINS
%	catdesc       = [ catdesc,       scm_bgain_vname, 'Histogram bins for gain factor.' ];
%	fieldnam      = [ fieldnam,      scm_bgain_vname, 'Gain bins' ];
%	fillval       = [ fillval,       scm_bgain_vname, single(-1e31) ];
%	format        = [ format,        scm_bgain_vname, 'f11.4' ];
%	si_conversion = [ si_conversion, scm_bgain_vname, ' ' ];
%	units         = [ units,         scm_bgain_vname, ' ' ];
%	validmin      = [ validmin,      scm_bgain_vname, single(-10.0) ];
%	validmax      = [ validmax,      scm_bgain_vname, single(10.0) ];
%	var_type      = [ var_type,      scm_bgain_vname, 'support_data' ];
	
	% PHASE SHIFT BINS
%	catdesc       = [ catdesc,       scm_boffset_vname, 'Histogram bins for phase shift.' ];
%	fieldnam      = [ fieldnam,      scm_boffset_vname, 'Phase shift bins' ];
%	fillval       = [ fillval,       scm_boffset_vname, single(-1e31) ];
%	format        = [ format,        scm_boffset_vname, 'f11.4' ];
%	si_conversion = [ si_conversion, scm_boffset_vname, '1e0>degrees' ];
%	units         = [ units,         scm_boffset_vname, 'degrees' ];
%	validmin      = [ validmin,      scm_boffset_vname, single(-180.0) ];
%	validmax      = [ validmax,      scm_boffset_vname, single(180.0) ];
%	var_type      = [ var_type,      scm_boffset_vname, 'support_data' ];
	
	% PSD RATIO BINS
%	catdesc       = [ catdesc,       scm_bpsdrat_vname, 'Histogram bins for PSD ratio.' ];
%	fieldnam      = [ fieldnam,      scm_bpsdrat_vname, 'PSD ratio bins' ];
%	fillval       = [ fillval,       scm_bpsdrat_vname, single(-1e31) ];
%	format        = [ format,        scm_bpsdrat_vname, 'f11.4' ];
%	si_conversion = [ si_conversion, scm_bpsdrat_vname, ' ' ];
%	units         = [ units,         scm_bpsdrat_vname, ' ' ];
%	validmin      = [ validmin,      scm_bpsdrat_vname, single(-10.0) ];
%	validmax      = [ validmax,      scm_bpsdrat_vname, single(10.0) ];
%	var_type      = [ var_type,      scm_bpsdrat_vname, 'support_data' ];
	
%-------------------------------------------------------------------------
% SCM METADATA INFO //////////////////////////////////////////////////////
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