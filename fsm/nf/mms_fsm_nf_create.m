%
% Name
%   mms_fsm_nf_create
%
% Purpose
%   Read noise floor results from individual orbits (srvy) or intervals (brst)
%   and compute the mean noise floor.
%
% Calling Sequence
%   NF = mms_fsm_nf_create( FILES )
%     Read data from either FGM or SCM background calibration files with names FILES
%     and determine the mean noise floor NF.
%
% Inputs
%   FILES           in, required, type=cell
%
% Returns
%   NF              out, required, type=struct
%                   Fields include:
%                     'f'     - Frequencies of noise floor
%                     'comp'  - Corresponding vector component
%                     'flag'  - Operational flags
%                     'nf'    - Noise floor determination
%                     'std'   - Standard deviation of the noise floor
%
% MATLAB release(s) MATLAB 9.0.0.341360 (R2016a)
% Required Products None
%
% History:
%   2016-10-20      Written by Matthew Argall
%
%***************************************************************************
function nf = mms_fsm_nf_create(files)
	
	if ischar(files)
		files = { files };
	end

	[sc, instr, mode, level, tstart, ~, optdesc] = mms_dissect_filename( files{1} );

%------------------------------------%
% Variable Names                     %
%------------------------------------%

	% Parse file name to get variable name prefix and suffix
	prefix    = [sc '_' instr '_'];
	suffix    = ['_' mode '_' level];
	mag_instr = optdesc(5:end);

	floor_vname = [prefix 'psd' '_floor' suffix];

%------------------------------------%
% Read Data                          %
%------------------------------------%
	nFiles     = length(files);
	noiseFloor = cell(1, nFiles);
	hFlag      = cell(1, nFiles);

	% Loop through each file.
	for ii = 1 : nFiles
		[noiseFloor{ii}, f, hFlag{ii}, comp] = MrCDF_Read(files{ii}, floor_vname, 'RowMajor', true);
	end
	
	% Concatenate data along flag dimension
	hFlag      = cat( 2, hFlag{:} );
	noiseFloor = cat( 1, noiseFloor{:} );

%------------------------------------%
% Setup Figure                       %
%------------------------------------%
	% Unique flags
	uniqFlags = unique(hFlag);
	nFlags    = length(uniqFlags);
	
	% Allocate memory
	nComp  = length(comp);
	nFreq  = length(f);
	nf     = zeros( nFlags, nComp, nFreq );
	nf_std = zeros( nFlags, nComp, nFreq );
	
	% Loop over each component
	for iComp = 1 : length(comp)
	for iFlag = 1 : nFlags
		% Elements with matching flags
		theFlag = uniqFlags(iFlag);
		idx     = find( hFlag == theFlag );

		% Mean and standard deviation
		nf(iFlag, iComp, :)     = mean( noiseFloor(idx,iComp,:), 1 );
		nf_std(iFlag, iComp, :) = std( noiseFloor(idx,iComp,:), 1, 1 );
	end
	end

%------------------------------------%
% Output Data Structure              %
%------------------------------------%
	% Reform so that
	%   [DEPEND_1, DEPEND_2, DEPEND_3, DEPEND_0]
	%   [   F,       FLAG,     COMP,     TIME  ]
	%   Matlab will remove a shallow trailing dimenion, so no need to reshape
	nf = struct( 'f',     f',         ...
	             'comp',  comp',      ...
	             'flag',  uniqFlags', ...
	             'nf',    single( permute( nf,     [3,2,1] ) ), ...
	             'std',   single( permute( nf_std, [3,2,1] ) )  ...
	           );
end
