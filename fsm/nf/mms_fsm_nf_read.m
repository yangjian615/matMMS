%
% Name
%   mms_fsm_nf_read
%
% Purpose
%   Read FSM noise floor data.
%
% Calling Sequence
%   DATA = mms_fsm_nf_read(FILES)
%     Read FSM noise floor data from files named FILES and return at
%     a data structure DATA.
%
% Parameters
%   FILES:          in, required, type=string/cell
%
% Returns
%   DATA            out, required, type=struct
%                   Structure with the following fields:
%                     't'     -  Time (cdf_time_tt2000)
%                     'dt'    -  Time interval over which the noise floor is valid.
%                     'nf'    -  Noise floor
%                     'std'   -  Standard deviation of the noise floor
%                     'flag'  -  Operational flag
%                     'comp'  -  Component index (0=x, 1=y, 2=z)
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2016-10-23      Written by Matthew Argall
%
function data = mms_fsm_nf_read(files, tstart, tend)

	if nargin() < 3
		tend = '';
	end
	if nargin() < 2
		tstart = '';
	end

	% Make sure we have a cell array of file names
	if ischar(files)
		files = { files };
	end
	nFiles = length(files);

%------------------------------------------------------
% Check File Compatibility                            |
%------------------------------------------------------
	
	% Dissect file name
	[sc, instr, mode, level, ~, ~, optdesc] = mms_dissect_filename( files{1} );

	% Check file names
	assert( all( strcmp(instr,   'fsm'    ) ), 'Only FSM files are allowed.' );
	if nFiles > 1
		assert( all( strcmp(mode,    mode{1}    ) ), 'All files must have the same telemetry mode.' );
		assert( all( strcmp(optdesc, optdesc{1} ) ), 'All files must have the same optional descriptor.' );
	end
	assert( all( strcmp(level,   'l2plus' ) ), 'Only L2PLUS files are allowed.' );

	% We now know all the files match, so keep on the the first value.
	if nFiles > 1
		sc      = sc{1};
		instr   = instr{1};
		mode    = mode{1};
		level   = level{1};
		optdesc = optdesc{1};
	end

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
% Read Data                                           |
%------------------------------------------------------
	% Read the magnetometer data
	[nf, t, f, comp, flag] = MrCDF_nRead(files, nf_vname,  'sTime', tstart, 'eTime', tend, 'RowMajor', true);
	dt                     = MrCDF_nRead(files, dt_vname,  'sTime', tstart, 'eTime', tend, 'RowMajor', true);
	std                    = MrCDF_nRead(files, std_vname, 'sTime', tstart, 'eTime', tend, 'RowMajor', true);

	if length(t) == 1
		f    = f';
		comp = comp';
		flag = flag';
	end
	
	% Transpose the data to be row vectors.
	data = struct( 't',    t,    ...
	               'dt',   dt,   ...
	               'f',    f,    ...
	               'nf',   nf,   ...
	               'std',  std,  ...
	               'flag', flag, ...
	               'comp', comp  ...
	             );
end