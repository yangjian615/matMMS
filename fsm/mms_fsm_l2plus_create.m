%
% Name
%   mms_fsm_create
%
% Purpose
%   Merge MMS fluxgate and search coil magnetometer data in the frequency
%   domain.
%
% Calling Sequence
%   [B_MERGED, T_MERGED] = mms_fsm_merge(SC. TSTART, TEND);
%     Given the MMS spacecraft number SC (e.g. 'mms1'), and a time
%     interval, [TSTART, TEND), gather all of the required search coil and
%     fluxgate magnetometer data and merge them into a single dataset
%     B_MERGED with time stamps T_MERGED.
%
%   [..., B_FG, T_FG] = mms_fsm_merge(__);
%     Also return the calibrated FGM magnetic field B_FG and its time
%     stamps T_FG.
%
%   [..., B_SC, T_SC] = mms_fsm_merge(__);
%     Also return the *UN*calibrated SCM magnetic field B_SC and its time
%     stamps T_SC.
%
% Parameters
%   SC:             in, required, type=char
%   TSTART:         in, required, type=char
%   TEND:           in, required, type=char
%   'Duration':     in, required, type=double
%                   The duration of each merging interval. Sets the
%                     frequency resolution of the final dataset.
%   'f_max':        in, required, type=double, default=Nyquist frequency
%                   The maximum of the frequency range to merge.
%   'f_min':        in, required, type=double, default=df
%                   The minimum ( > 0 ) of the frequency range to merge.
%   'fg_dir':       in, required, type=char, default=pwd();
%                   Directory in which to find FGM data.
%   'fg_cal_dir':   in, required, type=char, default=pwd();
%                   Directory in which to find FGM calibration data.
%   'sc_dir':       in, required, type=char, default=pwd();
%                   Directory in which to find SCM data.
%   'sc_cal_dir':   in, required, type=char, default=pwd();
%                   Directory in which to find SCM calibration data.
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-06      Written by Matthew Argall
%
function [t_fsm, b_fsm] = mms_fsm_l2plus_create( fgm, scm, duration, weight)

%------------------------------------%
% Check Inputs                       %
%------------------------------------%
	% Defaults
	if nargin < 3
		duration = [];
	end
	if nargin < 4
		weight = [];
	end
	nShift  = 0.25;
	nfilter = 4096;
	ndelay  = nfilter/2;
	fc      = 4.0;

%------------------------------------%
% Allocate Memory                    %
%------------------------------------%
	nPts   = length( scm.t );
	t_fill = MrCDF_Epoch_Compute( [9999, 12, 31, 0, 0, 0] );
	t_fsm  = repmat( t_fill, 1, nPts );
	b_fsm  = NaN( 3, nPts );
	
%------------------------------------%
% Find Unique Flag Values            %
%------------------------------------%
	% Unique Flags
	[flag_fgm, iFlag] = unique( fgm.flag );
	flag_scm          = unique( scm.flag );

	% Sort flags by time
	[~, isort] = sort( fgm.t(iFlag) );
	flag_fgm   = flag_fgm(isort);

	% Does scm contain slow srvy data
	%   - Where the Fast survey flag is not set.
	tf_slow = ~isempty( find( ~(bitget(flag_scm, 1)), 1, 'first' ) );
	if tf_slow
		is0_scm = 1;
		if0_scm = find( bitget(scm.flag, 1), 1, 'first' );
		is1_scm = if0_scm - 1;
		if1_scm = length( scm.t );
	else
		is0_scm = 1;
		is1_scm = length(scm.t);
		if0_scm = 1;
		if1_scm = is1_scm;
	end

%------------------------------------%
% Cycle Through Each FGM Flag        %
%------------------------------------%
	for ii = 1 : length(flag_fgm)
		% Indices associated with the current flag
		ifgm = find( fgm.flag == flag_fgm(ii) );
		
		% Weight function
%		iw = find( weight.flag == flag_fgm(ii) );
		
%		mrfprintf('logwarn', 'FIX: Hardwire index. Must regenerate cal files.');
%		iw = 1;

%		w  = weight.w(:,:,iw)';
%		f_w = weight.f;

	%------------------------------------%
	% Synchronize Data Intervals         %
	%------------------------------------%

		% Pick the appropriate SCM data
		%   - FGM bit 1 is for fast survey
		if bitget( flag_fgm(ii), 1 )
			i0 = if0_scm;
			i1 = if1_scm;
		else
			i0 = is0_scm;
			i1 = is1_scm;
		end
		
		% Syncronize intervals
		[iFGM, iSCM] = MrIntervalsXY( fgm.t(ifgm), scm.t(i0:i1), ...
		                              'Remove', true, ...
		                              'Sync',   true, ...
		                              'TolX',   2,    ...
		                              'TolY',   2 );

	%------------------------------------%
	% Loop Through Data Intervals        %
	%------------------------------------%
		% Number of data intervals found
		nIntervals = size(iFGM, 2);
		temp_fgm   = struct([]);
		temp_scm   = struct([]);
		for jj = 1 : nIntervals
			% Subinterval to be processed
			i0_fgm = ifgm(1) + iFGM(1,jj) - 1;
			i1_fgm = ifgm(1) + iFGM(2,jj) - 1;
			i0_scm = i0      + iSCM(1,jj) - 1;
			i1_scm = i0      + iSCM(2,jj) - 1;

		%------------------------------------%
		% Interpolate FGM to SCM             %
		%------------------------------------%

		%------------------------------------%
		% Apply the Model                    %
		%------------------------------------%
			%
			% TODO:
			%   1) Adjust filter model for 8, 16, 32, 64 Hz data
			%   2) Apply orthogonality correction between FGM and SCM
			%   3) Improve frequency-domain merging filter
			%   4) Create time-domain merging filter
			%   5) Apply SCM frequency compensation in time domain
			%

			% Find the model
			range = bitget(flag_fgm(ii), 2);
			mode  = bitget(flag_fgm(ii), 3);
			fs    = fgm.sr(i0_fgm);
			fsnew = scm.sr(i0_scm);
			model = mms_fsm_fgm_load_model(fgm.model_dir, fgm.sc, fgm.instr, range, mode);
			delay = mms_fsm_fgm_delay( fgm.sc, fgm.instr, fs, mode);

			% Apply frequency compensation
			[t_fgm, b_fgm] = mms_fsm_fgm_compensate( model, fgm.t(i0_fgm:i1_fgm)', fgm.b(:,i0_fgm:i1_fgm)', ...
			                                         delay, fs, fsnew );
			
			% Fix delay between FGM & SCM
			[t_fgm, b_fgm, t_scm, b_scm] = mms_fsm_align_time( model, t_fgm, b_fgm, scm.t(i0_scm:i1_scm)', scm.b(:,i0_scm:i1_scm)', fsnew );
			i1_scm = i0_scm + length(t_scm) - 1;
	
		%------------------------------------%
		% Merge Datasets                     %
		%------------------------------------%

			% Merge the data
			for kk = 1 : 3
				b_fsm_fgm = MrFilter_WinSinc( b_fgm(:,kk)', 1/fsnew, fc, nfilter+1, nfilter, 'blackman', 'low' );
				b_fsm_scm = MrFilter_WinSinc( b_scm(:,kk)', 1/fsnew, fc, nfilter+1, nfilter, 'blackman', 'high' );
				
				% Remove the delay
				b_fsm_fgm = b_fsm_fgm(ndelay+1:end);
				b_fsm_scm = b_fsm_scm(ndelay+1:end);

				% Store the signal
				b_fsm(kk,i0_scm:i1_scm-ndelay) = b_fsm_fgm + b_fsm_scm;
			end
			
			% Time array
			t_fsm( i0_scm:i1_scm-ndelay ) = t_scm(1:end-ndelay)';
			
			clear b_fsm_fgm b_fsm_scm
			
%			try
%				b_fsm(:, i0_scm:i1_scm ) = ...
%					fsm_merge_interp_nocal_shift( b_fgm', b_scm', fsnew, duration, nShift, w, f_w );
%			catch ME
%				mrfprintf('logerr', ME);
%				mrfprintf('logtext', 'Skipping interval %i of %i for flag %i.', jj, nIntervals, flag_fgm(ii));
%			end
			
			clear t_fgm b_fgm t_scm b_scm
		end
	end

	% Remove data not included in the merging process
	igood = find( t_fsm ~= t_fill );
	t_fsm = t_fsm(igood);
	b_fsm = b_fsm(:,igood);
end
