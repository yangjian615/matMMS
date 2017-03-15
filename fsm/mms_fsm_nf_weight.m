%
% Name
%   mms_fsm_nf_create
%
% Purpose
%   Create a weight function for merging FGM and SCM datasets.
%
% Calling Sequence
%   W = mms_fsm_nf_create( NF_FGM, NF_SCM, T0 )
%     Create a weight function W using noise floor estimates for FGM and SCM
%     (NF_FGM, NF_SCM, respectively) obtained at the time T0.
%
% Inputs
%   NF_FGM          in, required, type=struct
%                   Fields include:
%                     'f'     - Frequencies of noise floor
%                     'comp'  - Corresponding vector component
%                     'flag'  - Operational flags
%                     'nf'    - Noise floor determination
%                     'std'   - Standard deviation of the noise floor
%   NF_SCM          in, required, type=struct
%                   Fields include:
%                     'f'     - Frequencies of noise floor
%                     'comp'  - Corresponding vector component
%                     'flag'  - Operational flags
%                     'nf'    - Noise floor determination
%                     'std'   - Standard deviation of the noise floor
%   T0              out, required, type=int64 (cdf_time_tt2000)
%
% Returns
%   W               out, required, type=structure
%
% MATLAB release(s) MATLAB 9.0.0.341360 (R2016a)
% Required Products None
%
% History:
%   2016-10-25      Written by Matthew Argall
%
%***************************************************************************
function weight = mms_fsm_nf_weight(nf_fgm, nf_scm, t0)

	% Low-frequency cut-off when SCM should not be used
	%   - SCM has a low-pass filter at these frequencies
	fc_srvy = 0.5;
	fc_brst = 1.0;

%------------------------------------%
% Check Inputs                       %
%------------------------------------%
	
	N_fgm = length(nf_fgm.t);
	N_scm = length(nf_scm.t);
	
	% Upper time limit
	t_plus_fgm = nf_fgm.t + nf_fgm.dt;
	t_plus_scm = nf_scm.t + nf_scm.dt;
	
	% Locate relevant noise floor
	it_fgm = find( nf_fgm.t >= t0 & t_plus_fgm > t0 );
	it_scm = find( nf_scm.t >= t0 & t_plus_scm > t0 );
	
	% More than one found
	%   - Pick the first (closest start time)
	if length(it_fgm) > 1
		mrfprintf( 'logwarn', 'FGM noise floor entries found: %i', length(it_fgm) )
		it_fgm = it_fgm(1);
	end
	if length(it_scm) > 1
		mrfprintf( 'logwarn', 'SCM noise floor entries found: %i', length(it_scm) )
		it_scm = it_scm(1);
	end

%------------------------------------%
% Separate Data                      %
%------------------------------------%
	
	% Create weight function by flag
	%   - SCM has only slow/fast
	%   - FGM has slow/fast; pre/post-perigee; hi/lo-range; DEC 32/64
	%   - Apply SCM slow to all FGM slow (regardless of other flags)
	%   - Apply SCM fast to all FGM fast (regardless of other flags)
	%   - Slow/fast is bit 1
	tf_fast_fgm = bitget( nf_fgm.flag, 1 );
	tf_fast_scm = bitget( nf_scm.flag, 1 );

	% Slow/fast survey
	islow_fgm = find( ~tf_fast_fgm );
	islow_scm = find( ~tf_fast_scm );
	ifast_fgm = find( tf_fast_fgm );
	ifast_scm = find( tf_fast_scm );
	
	% For days after about 2015-09-01, SCM operates only in fast survey
	if isempty(islow_scm)
		islow_scm = ifast_scm;
	end

%------------------------------------%
% Slow Srvy Weight Function          %
%------------------------------------%
	% Allocate memory
	nFreq = length( nf_scm.f );
	nFlag = length( nf_fgm.flag );
	w     = zeros(nFreq, 3, nFlag);
	
	% Transition to pure SCM
	ifc = find( nf_scm.f >= fc_srvy, 1, 'first' );
	ifN = find( nf_scm.f >= 4.0,     1, 'first' );

	% W: Slow survey
	for ii = 1 : length(islow_fgm)
		iFlag = islow_fgm(ii);
		
		% Weighting factors
		temp_fgm = squeeze( 1.0 ./ nf_fgm.nf(ifc:ifN, :, iFlag,     it_fgm) );
		temp_scm = squeeze( 1.0 ./ nf_scm.nf(ifc:ifN, :, islow_scm, it_scm) );
		
		% Weight function
		%   - Of the total noise floor NF, W is the % of NF made up by SCM.
		%   - This means FGM is better than SCM by that percentage
		%   - Multiply the FGM signal by W and SCM signal by (1-W)
		w(ifc:ifN,   :, iFlag) = temp_scm ./ (temp_fgm + temp_scm);
		w(1:ifc-1,   :, iFlag) = 1.0;
		w(ifN+1:end, :, iFlag) = 0.0;
	end

%------------------------------------%
% Fast/Brst Weight Function          %
%------------------------------------%
	
	% Transition to pure SCM
	if max(nf_scm.f) == 16
		ifc = find( nf_scm.f >= fc_srvy, 1, 'first' );
		ifN = find( nf_scm.f >= 8,       1, 'first' );
	else
		ifc = find( nf_scm.f >= fc_brst, 1, 'first' );
		ifN = find( nf_scm.f >= 16,      1, 'first' );
	end
	
	% W: Fast survey
	for ii = 1 : length(ifast_fgm)
		iFlag = ifast_fgm(ii);
		
		% Weighting factor
		temp_fgm = squeeze( 1.0 ./ nf_fgm.nf(ifc:ifN, :, iFlag,     it_fgm) );
		temp_scm = squeeze( 1.0 ./ nf_scm.nf(ifc:ifN, :, islow_scm, it_scm) );
		
		% Weight function
		w(ifc:ifN,   :, iFlag) = temp_scm ./ (temp_fgm + temp_scm);
		w(1:ifc-1,   :, iFlag) = 1.0;
		w(ifN+1:end, :, iFlag) = 0.0;
	end

%------------------------------------%
% Output Structure                   %
%------------------------------------%
	weight = struct( 'flag', nf_fgm.flag, ...
	                 'f',    nf_scm.f,    ...
	                 'w',    w );
end
