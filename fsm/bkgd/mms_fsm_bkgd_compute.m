%
% Name
%   mms_fsm_bkgd_compute
%
% Purpose
%   Compute the amplitude, phase, and power of magnetic field measurments, histogram
%   them over time, then find the noise floor based on the distribution of amplitude
%   at each frequency.
%
% Calling Sequence
%   [FGM, SCM] = mms_fsm_bkgd_compute(FGM, SCM, T)
%     Use the data structures FGM and SCM containing time and magnetic field
%     wave forms to compute the amplitude, power, and phase over a sliding
%     window of duration T seconds. Append the results to the input structures.
%
% Parameters
%   FGM             in, required, type = struct
%                   Contains fields:
%                       't' - Time (TT2000)
%                       'b' - 3-component magnetic field
%   SCM             in, required, type = struct
%                   Contains fields:
%                       't' - Time (TT2000)
%                       'b' - 3-component magnetic field
%   T               in, required, type = double
%
% Returns
%   BKGD_FGM        out, required, type=struct
%                   Fields are:
%                       'amp'         - Wave amplitude
%                       'phase'       - Phase
%                       'psd'         - Power spectral density
%                       'amp_hist'    - 
%                       'phase_hist'  - 
%                       'psd_hist'    - 
%                       'amp_floor'   - 
%                       'phase_floor' - 
%                       'psd_floor'   - 
%                       'flag'        - 
%                       'hist_flag'   - 
%   BKGD_SCM        out, required, type=struct
%                   Fields are:
%                       'amp'   - Wave amplitude
%                       'phase' - Phase
%                       'psd'   - Power spectral density
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2016-05-31      Written by Matthew Argall
%
function [bkgd_fgm, bkgd_scm, bkgd_cross] = mms_fsm_bkgd_compute(fgm, scm, T)

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
	tf_slow = ~isempty( find( ~(bitget(flag_scm, 1)) ) );
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
% Cycle Throught Each FGM Flag       %
%------------------------------------%
	for ii = 1 : length(flag_fgm)
		% Indices associated with the current flag
		ifgm = find( fgm.flag == flag_fgm(ii) );

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
		temp_x     = struct([]);
		for jj = 1 : nIntervals
			% Subinterval to be processed
			i0_fgm = ifgm(1) + iFGM(1,jj) - 1;
			i1_fgm = ifgm(1) + iFGM(2,jj) - 1;
			i0_scm = i0      + iSCM(1,jj) - 1;
			i1_scm = i0      + iSCM(2,jj) - 1;
	
		%------------------------------------%
		% FGM Freq Compensation              %
		%------------------------------------%
			%
			% TODO: Separate brst & srvy data
			%
			
			% Find the model
			range = bitget(flag(ii), 2);
			mode  = bitget(flag(ii), 3);
			fs    = fgm.sr(i0);
			fsnew = scm.sr(i0_scm);
			model = mms_fsm_fgm_load_model( fgm.sc, fgm.instr, range, mode );
			delay = mms_fsm_fgm_delay( fgm.sc, fgm.instr, fs, mode);
			
			% Apply frequency compensation
			[t_fgm, b_fgm] = mms_fsm_fgm_compensate( model, fgm.t(i0_fgm:i1_fgm)', fgm.b(:,i0_fgm:i1_fgm)', ...
			                                         delay, fs, fsnew );
			
			% Transpose back
			t_ref     = t_fgm(1);
			t_fgm_ssm = MrCDF_epoch2ssm(t_fgm', t_ref);
			b_fgm     = b_fgm';
	
		%------------------------------------%
		% SCM Data                           %
		%------------------------------------%
			t_scm_ssm = MrCDF_epoch2ssm( scm.t( i0_scm:i1_scm ), t_ref );
			b_scm     = scm.b(:, i0_scm:i1_scm);
	
		%------------------------------------%
		% Interpolate                        %
		%------------------------------------%
			
			% Convert time to seconds since midnight
			%   - Doubles needed for interpolation
%			t_ref     = scm.t( i0_scm );
%			t_fgm_ssm = MrCDF_epoch2sse( fgm.t( i0_fgm:i1_fgm ), t_ref);
%			t_scm_ssm = MrCDF_epoch2sse( scm.t( i0_scm:i1_scm ), t_ref);
	
			% Intepolate FGM to SCM
%			b_fgm_interp      = zeros(3, i1_scm-i0_scm+1);
%			b_fgm_interp(1,:) = interp1(t_fgm_ssm, fgm.b(1, i0_fgm:i1_fgm), t_scm_ssm, 'linear', 'extrap');
%			b_fgm_interp(2,:) = interp1(t_fgm_ssm, fgm.b(2, i0_fgm:i1_fgm), t_scm_ssm, 'linear', 'extrap');
%			b_fgm_interp(3,:) = interp1(t_fgm_ssm, fgm.b(3, i0_fgm:i1_fgm), t_scm_ssm, 'linear', 'extrap');
	
		%------------------------------------%
		% Compute Spectral Components        %
		%------------------------------------%
			
			% Compute amplitude, phase, power
%			temp_fgm_jj = mms_fsm_bkgd_spectra(t_fgm_ssm, fgm.b(:,i0_fgm:i1_fgm), T);
			temp_fgm_jj = mms_fsm_bkgd_spectra(t_scm_ssm, b_fgm, T);
			temp_scm_jj = mms_fsm_bkgd_spectra(t_scm_ssm, b_scm, T);
			
			% Gain, Offset, Power Ratio
			temp_x_jj.('t')           = temp_scm_jj.t;
			temp_x_jj.('f')           = temp_scm_jj.f;
			temp_x_jj.('gain')        = temp_fgm_jj.('amp') ./ temp_scm_jj.('amp');
			temp_x_jj.('phase_shift') = temp_fgm_jj.('phase') - temp_scm_jj.('phase');
			temp_x_jj.('psd_rat')     = temp_fgm_jj.('psd') ./ temp_scm_jj.('psd');
	
		%------------------------------------%
		% Histogram Spectral Components      %
		%------------------------------------%
			
			% Enough data?
			%   - FGM and SCM should have identical number of points and time stamps
			if ~isempty(temp_fgm_jj.t)
				% Histogram spectral components
				temp_fgm_jj = mms_fsm_bkgd_histogram(temp_fgm_jj);
				temp_scm_jj = mms_fsm_bkgd_histogram(temp_scm_jj);
				temp_x_jj   = mms_fsm_bkgd_histogram(temp_x_jj);
				
				% Finalize structure
				%  - Convert time back to TT2000
				%  - Include flag
				temp_fgm_jj.('t')         = MrCDF_sse2epoch(temp_fgm_jj.t, t_ref);
				temp_fgm_jj.('flag')      = repmat( flag_fgm(ii), 1, length(temp_fgm_jj.t) );
				temp_fgm_jj.('hist_flag') = flag_fgm(ii);
				
				temp_scm_jj.('t')         = MrCDF_sse2epoch(temp_scm_jj.t, t_ref);
				temp_scm_jj.('flag')      = repmat( flag_fgm(ii), 1, length(temp_scm_jj.t) );
				temp_scm_jj.('hist_flag') = flag_fgm(ii);
				
				temp_x_jj.('t')           = MrCDF_sse2epoch(temp_x_jj.t, t_ref);
				temp_x_jj.('flag')        = repmat( flag_fgm(ii), 1, length(temp_fgm_jj.t) );
				temp_x_jj.('hist_flag')   = flag_fgm(ii);
				
				% Sum Histogram Occurrence over Multiple Intervals
				if isempty(temp_fgm)
					temp_fgm   = temp_fgm_jj;
					temp_scm   = temp_scm_jj;
					temp_cross = temp_x_jj;
				else
					temp_fgm.('amp_hist')   = temp_fgm.('amp_hist')   + temp_fgm_jj.('amp_hist');
					temp_fgm.('phase_hist') = temp_fgm.('phase_hist') + temp_fgm_jj.('phase_hist');
					temp_fgm.('psd_hist')   = temp_fgm.('psd_hist')   + temp_fgm_jj.('psd_hist');
					
					temp_scm.('amp_hist')   = temp_scm.('amp_hist')   + temp_scm_jj.('amp_hist');
					temp_scm.('phase_hist') = temp_scm.('phase_hist') + temp_scm_jj.('phase_hist');
					temp_scm.('psd_hist')   = temp_scm.('psd_hist')   + temp_scm_jj.('psd_hist');
					
					temp_cross.('gain_hist')        = temp_cross.('gain_hist')        + temp_x_jj.('gain_hist');
					temp_cross.('phase_shift_hist') = temp_cross.('phase_shift_hist') + temp_x_jj.('phase_shift_hist');
					temp_cross.('psd_rat_hist')     = temp_cross.('psd_rat_hist')     + temp_x_jj.('psd_rat_hist');
				end
			end
			
			clear temp_fgm_jj temp_scm_jj
		end
		
	%------------------------------------%
	% Determine Noise Floor              %
	%------------------------------------%
		% Amplitude, Power, Phase
%		temp_fgm.('amp_floor')   = mms_fsm_bkgd_noisefloor( temp_fgm.('amp_hist'),   temp_fgm.('amp_bins')   );
%		temp_fgm.('phase_floor') = mms_fsm_bkgd_noisefloor( temp_fgm.('phase_hist'), temp_fgm.('phase_bins') );
%		temp_fgm.('psd_floor')   = mms_fsm_bkgd_noisefloor( temp_fgm.('psd_hist'),   temp_fgm.('psd_bins')   );
		
%		temp_scm.('amp_floor')   = mms_fsm_bkgd_noisefloor( temp_scm.('amp_hist'),   temp_scm.('amp_bins')   );
%		temp_scm.('phase_floor') = mms_fsm_bkgd_noisefloor( temp_scm.('phase_hist'), temp_scm.('phase_bins') );
%		temp_scm.('psd_floor')   = mms_fsm_bkgd_noisefloor( temp_scm.('psd_hist'),   temp_scm.('psd_bins')   );

%		temp_scm.('gain_floor')        = mms_fsm_bkgd_noisefloor( temp_scm.('gain_hist'),        temp_scm.('gain_bins')   );
%		temp_scm.('phase_shift_floor') = mms_fsm_bkgd_noisefloor( temp_scm.('phase_shift_hist'), temp_scm.('phase_shift_bins') );
%		temp_scm.('psd_rat_floor')     = mms_fsm_bkgd_noisefloor( temp_scm.('psd_rat_hist'),     temp_scm.('psd_rat_bins')     );
	
	%------------------------------------%
	% Save Data from Multiple Flags      %
	%------------------------------------%
		if ii == 1
			bkgd_fgm   = temp_fgm;
			bkgd_scm   = temp_scm;
			bkgd_cross = temp_cross;
		else
			bkgd_fgm   = [bkgd_fgm    temp_fgm];
			bkgd_scm   = [bkgd_scm    temp_scm];
			bkgd_cross = [bkgd_crosss temp_cross];
		end
		
		clear temp_fgm temp_scm temp_cross
	end

%------------------------------------%
% Concatenate Multiple Tags          %
%------------------------------------%
	% FGM and SCM should have the same number of structure array elements
	%   - Concatenate along the time or hist_flag dimension
	if length(bkgd_fgm) > 1
		bkgd_fgm   = MrCatStruct(bkgd_fgm,   [2,0,2,2,2,2,0,2,0,2,0,2,2]);
		bkgd_scm   = MrCatStruct(bkgd_scm,   [2,0,2,2,2,2,0,2,0,2,0,2,2]);
		bkgd_cross = MrCatStruct(bkgd_cross, [2,0,2,2,2,2,0,2,0,2,0,2,2]);
	end
end