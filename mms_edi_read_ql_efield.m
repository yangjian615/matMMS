%
% Name
%   mms_edi_read_ql_efield
%
% Purpose
%   Read Quick-Look (QL) electron drift instrument (EDI) data from MMS.
%
% Calling Sequence
%   EDI_QL = mms_edi_read_ql_efield(FILES)
%     Read electron drift instrument (EDI) Quick-Look (QL) data from
%     files with file names FILES.
%
%   EDI_QL = mms_edi_read_ql_efield(FILES, TSTART, TEND)
%     Read electron drift instrument (EDI) Quick-Look (QL) data from
%     files with file names FILES between the time interval beginning at
%     TSTART and ending at TEND. TSTART and TEND should be ISO-8601
%     strings: yyyy-mm-ddThh:mm:ss
%
% Parameters
%   FILES           in, required, type = char/cell
%   TSTART          in, required, type = char, default = ''
%   TEND            in, required, type = char, default = ''
%
% Returns
%   EDI_QL          out, required, type=struct
%                   Fields are:
%                       'tt2000'      - TT2000 epoch times
%                       'E_DSL'       - Electric field in DSL coordinates
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-08-13      Written by Matthew Argall
%   2015-08-22      Updated to read file version 0.2.0. - MRA
%
function edi_ql = mms_edi_read_ql_efield(files, tstart, tend)

	% Defaults
	if nargin() < 2
		tstart = '';
	end
	if nargin() < 3
		tend = '';
	end

%------------------------------------%
% Check File Names                   %
%------------------------------------%
	% Number of files given
	if iscell(files)
		nFiles = length(files);
	else
		nFiles = 1;
	end

	% Dissect the file name
	[sc, instr, mode, level, ~, ~, optdesc] = mms_dissect_filename(files);

	% Instr, Level, Mode
	assert( min( strcmp(instr, 'edi') == 1 ), 'Only EDI files are allowed.' );
	assert( min( strcmp(level, 'ql') ) == 1, 'Only QL files are allowed.' );
	if nFiles > 1
		assert( min( strcmp(instr, instr{1}) ) == 1, 'All files must be from the same instrument.' );
		assert( min( strcmp(mode, mode{1}) )   == 1, 'All files must have the same telemetry mode.' );
	end

	% We now know all the files match, so keep on the the first value.
	if nFiles > 1
		sc      = sc{1};
		instr   = instr{1};
		mode    = mode{1};
		level   = level{1};
		optdesc = optdesc{1};
	end

%------------------------------------%
% Read EDI Data                      %
%------------------------------------%

	% Create variable names
	E_vname           = mms_construct_varname(sc, instr, 'E',            'dmpa');
	v_vname           = mms_construct_varname(sc, instr, 'v_ExB',        'dmpa');
	d_vname           = mms_construct_varname(sc, instr, 'd',            'dmpa');
	E_bc_vname        = mms_construct_varname(sc, instr, 'E',            'bc_dmpa');
	v_bc_vname        = mms_construct_varname(sc, instr, 'v_ExB',        'bc_dmpa');
	d_bc_vname        = mms_construct_varname(sc, instr, 'd',            'bc_dmpa');
	q_bc_vname        = mms_construct_varname(sc, instr, 'quality',      'bc');
	b_vname           = mms_construct_varname(sc, instr, 'B',            'dmpa');
	b_std_vname       = mms_construct_varname(sc, instr, 'B',            'std_dmpa');
	recnum_vname      = mms_construct_varname(sc, instr, 'recnum');
	recnum_gd12_vname = mms_construct_varname(sc, instr, 'recnum',       'gd12');
	recnum_gd21_vname = mms_construct_varname(sc, instr, 'recnum',       'gd21');
	fv_gd12_vname     = mms_construct_varname(sc, instr, 'fv',           'gd12_dmpa');
	fv_gd21_vname     = mms_construct_varname(sc, instr, 'fv',           'gd21_dmpa');
	q_gd12_vname      = mms_construct_varname(sc, instr, 'beam_quality', 'gd12');
	q_gd21_vname      = mms_construct_varname(sc, instr, 'beam_quality', 'gd21');

	% Read the magnetometer data
	[E_dmpa, t]  = MrCDF_nRead(files, E_vname,           'sTime', tstart, 'eTime', tend);
	v_dmpa       = MrCDF_nRead(files, v_vname,           'sTime', tstart, 'eTime', tend);
	d_dmpa       = MrCDF_nRead(files, d_vname,           'sTime', tstart, 'eTime', tend);
	E_bc_dmpa    = MrCDF_nRead(files, E_bc_vname,        'sTime', tstart, 'eTime', tend);
	v_bc_dmpa    = MrCDF_nRead(files, v_bc_vname,        'sTime', tstart, 'eTime', tend);
	d_bc_dmpa    = MrCDF_nRead(files, d_bc_vname,        'sTime', tstart, 'eTime', tend);
	q_bc_dmpa    = MrCDF_nRead(files, q_bc_vname,        'sTime', tstart, 'eTime', tend);
	b_dmpa       = MrCDF_nRead(files, b_vname,           'sTime', tstart, 'eTime', tend);
	b_std_dmpa   = MrCDF_nRead(files, b_std_vname,       'sTime', tstart, 'eTime', tend);
	recnum       = MrCDF_nRead(files, recnum_vname,      'sTime', tstart, 'eTime', tend);
	recnum_gd12  = MrCDF_nRead(files, recnum_gd12_vname, 'sTime', tstart, 'eTime', tend);
	recnum_gd21  = MrCDF_nRead(files, recnum_gd21_vname, 'sTime', tstart, 'eTime', tend);
	fv_gd12_dmpa = MrCDF_nRead(files, fv_gd12_vname,     'sTime', tstart, 'eTime', tend);
	fv_gd21_dmpa = MrCDF_nRead(files, fv_gd21_vname,     'sTime', tstart, 'eTime', tend);
	q_gd12       = MrCDF_nRead(files, q_gd12_vname,      'sTime', tstart, 'eTime', tend);
	q_gd21       = MrCDF_nRead(files, q_gd21_vname,      'sTime', tstart, 'eTime', tend);

	% Transpose the data to be row vectors.
	edi_ql = struct( 'tt2000',        t,             ...
	                 'E_dmpa',        E_dmpa,        ...
	                 'v_ExB_dmpa',    v_dmpa,        ...
	                 'd_dmpa',        d_dmpa,        ...
	                 'E_bc_dmpa',     E_bc_dmpa,     ...
	                 'v_ExB_bc_dmpa', v_bc_dmpa,     ...
	                 'd_bc_dmpa',     d_bc_dmpa,     ...
	                 'q_bc_dmpa',     q_bc_dmpa,     ...
	                 'B_dmpa',        b_dmpa,        ...
	                 'B_std_dmpa',    b_std_dmpa,    ...
	                 'recnum',        recnum,        ...
	                 'recnum_gd12',   recnum_gd12,   ...
	                 'recnum_gd21',   recnum_gd21,   ...
	                 'fv_gd12_dmpa',  fv_gd12_dmpa,  ...
	                 'fv_gd21_dmpa',  fv_gd21_dmpa,  ...
	                 'q_gd12',        q_gd12,        ...
	                 'q_gd21',        q_gd21         ...
	               );
end