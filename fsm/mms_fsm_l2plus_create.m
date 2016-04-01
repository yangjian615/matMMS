%
% Name
%   mms_fsm_create_l1b
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
function [t_fsm, b_fsm] = mms_fsm_l2plus_create(fgm_file, scm_file, scm_cal_file, zmpa, duration, f_min, f_max)

	if nargin < 6
		f_min = 1.5;
	end
	if nargin < 7
		f_max = 4.0;
	end

%------------------------------------%
% Prepare Data                       %
%------------------------------------%

	% Read SCM and transfer function
	[tt2000_scm, b_scm, xfr_scm, sr_scm] = mms_fsm_l2plus_scm_prep(scm_file, scm_cal_file, duration);

	% Read and calibrate FGM
	[tt2000, b_fgm_omb, b_scm_omb] = mms_fsm_l2plus_fgm_prep(fgm_file, zmpa, tt2000_scm, b_scm);

	clear tt2000_scm b_scm

%------------------------------------%
% Find Coninuous, Overlapping Data   %
%------------------------------------%
	% Continuous data intervals
	t_ssm       = MrCDF_epoch2ssm(tt2000);
	[int, nint] = MrIntervalsX(t_ssm);

%------------------------------------%
% Loop Over Intervals                %
%------------------------------------%

	%
	% TODO
	%   1) Noise floor
	%

	% Allocate memory to output
	t_fsm = zeros(1, size(b_scm_omb,2), 'int64');
	b_fsm = zeros(size(b_scm_omb), 'single');

	% Step through each interval
	for ii = 1 : nint 
		% Absolute start and end indices.
		istart = int(1,ii);
		iend   = int(2,ii);
	
		% Extract the data for the current interval
		time  = tt2000(istart:iend);
		b_fgm = b_fgm_omb(:, istart:iend);
		b_scm = b_scm_omb(:, istart:iend);

		% Merge the data
		t_fsm(istart:iend) = time;
		try
			b_fsm(:, istart:iend) = ...
				fsm_merge_vinterp(b_fgm, b_scm, xfr_scm, sr_scm, duration);
		catch ME
			mrfprintf('logerr', ME);
			mrfprintf('logtext', 'Skipping interval %i of %i.', ii, nint);
			b_merged(:, istart:iend) = NaN;
		end
	end
	clear t_fgm t_scm b_fgm b_scm sr_fgm sr_scm xfr_scm

	% Remove data not included in the merging process
	igood = find( t_fsm ~= 0 );
	t_fsm = t_fsm(igood);
	b_fsm = b_fsm(:,igood);
end



%
% Name
%   mms_fsm_ql_prep_fgm
%
% Purpose
%   Read FGM L1A data and calibrate it, creating a data product in OMB coordinates.
%
% Calling Sequence
%   [TT2000, B_BCS] = mms_fsm_ql_prep_fgm(FILES, T_SCM, TSTART, TEND);
%     Using FGM L2PRE data files, named in the structure FILES, in the
%     time interval [TSTART, TEND), read data and upsample to SCM times,
%     given by T_SCM.
%
% Parameters
%   FILES:          in, required, type=struct
%                   Required fields:
%                       'fgm'        - FGM file name
%   T_SCM:          in, required, type=int64 (cdf_time_tt2000)
%   TSTART:         in, required, type=char
%   TEND:           in, required, type=char
%
% Returns
%   TT2000          out, required, type=1xN INT64 (cdf_time_tt2000)
%   B_BCS           out, required, type=3xN float
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-11-30      Written by Matthew Argall
%
function [tt2000, b_omb, b_scm_out] = mms_fsm_l2plus_fgm_prep(file, zmpa, t_scm, b_scm)

%------------------------------------%
% Read                               %
%------------------------------------%
	% Dissect the file name
	[~, ~, mode] = mms_dissect_filename(file);

	%
	% L2Pre BCS data is already calibrated. We are going to up-sample
	% to SCM time stamps, but that must be done at a later stage,
	% after data gaps have been removed.
	%

	% Uncalibrated FG in 123
	fgm_l2pre  = mms_fg_read_l2pre(file);
	
	% Extract data
	tt2000 = fgm_l2pre.tt2000;
	b_bcs  = fgm_l2pre.b_bcs(1:3,:);
	clear fgm_l2pre

%------------------------------------%
% Separate SLOW, FAST, BRST          %
%------------------------------------%
	%
	% In order to find gaps in the FGM data, we must know
	% what the sampling rate is. In the SRVY data, SLOW
	% and FAST data are combined. Our first step is to
	% separate them. We then find data gaps and interpolate
	% over continuous, overlapping intervals.
	%

	% Convert to seconds
	t_fgm_ssm = MrCDF_epoch2ssm(tt2000);
	
	% Compute the sampling rate
	%   - integer values
	%   - assume (t_0 - t_(-1)) = (t1 - t0)
	sr = round( 1.0 ./ ( diff(t_fgm_ssm) ) );
	sr = [sr(1) sr];
	
	% Find unexpected sampling rates
	%   slow = 8
	%   fast = 16
	%   brst = 128
	ibad = find( sr ~= 8 & sr ~= 16 & sr ~= 128 );
	if ~isempty(ibad)
		nbad = length(ibad);
		mrfprintf('logwarn', '%d FGM points found with unexpected sampling rate.', nbad);
		for ii = 1 : min( [10, nbad] )
			mrfprintf('logtext', '    Sampling rate: %d', ibad(ii));
		end
	end

	% Separate FAST from SLOW
	if strcmp(mode, 'brst')
		ifast = find( sr == 128 );
		islow = [];
	else
		islow = find( sr ==  8 );
		ifast = find( sr == 16 );
	end
	
	%
	% TODO: Split HIRANGE and LORANGE
	%

%------------------------------------%
% Interpolate to SCM                 %
%------------------------------------%
	% Interpolate slow and fast data
	if ~isempty(ifast)
		[t_fast, b_fast, b_scm_fast] ...
			= mms_fsm_l2plus_fgm_interp(tt2000(ifast), b_bcs(:,ifast), t_scm, b_scm);
	end
	if ~isempty(islow)
		[t_slow, b_slow, b_scm_slow] ...
			= mms_fsm_l2plus_fgm_interp(tt2000(islow), b_bcs(:,islow), t_scm, b_scm);
	end
	
%------------------------------------%
% SRVY Data                          %
%------------------------------------%
	
	%
	% TODO: Do fast and slow srvy data overlap?
	%

	% Combine fast and slow now that they have the same sampling rates
	if strcmp(mode, 'brst')
		tt2000    = t_fast;
		b_bcs     = b_fast;
		b_scm_out = b_scm_fast;
	else
		% Combine slow and fast into srvy
		tt2000    = [t_slow, t_fast];
		b_bcs     = [b_slow, b_fast];
		b_scm_out = [b_scm_slow, b_scm_fast];
	
		% Sort in time
		[tt2000, isort] = sort( tt2000 );
		b_bcs           = b_bcs(:, isort);
		b_scm_out       = b_scm_out(:, isort);
	end
	
	% Delete data
	clear t_slow t_fast b_slow b_fast
	
%------------------------------------%
% BCS --> OMB                        %
%------------------------------------%
	% Transformation matrixes
	bcs2smpa = mms_fg_xbcs2smpa(zmpa);
	omb2smpa = mms_fg_xomb2smpa();
	smpa2omb = permute(omb2smpa, [2,1]);
	
	% Rotate to OMB
	b_smpa   = mrvector_rotate(bcs2smpa, b_bcs);
	b_omb    = mrvector_rotate(smpa2omb, b_smpa);
end


%
% Name
%   mms_fsm_l2plus_fgm_interp
%
% Purpose
%   Read FGM L1A data and calibrate it, creating a data product in OMB coordinates.
%
% Calling Sequence
%   [TT2000, B_BCS] = mms_fsm_ql_prep_fgm(FILES, T_SCM, TSTART, TEND);
%     Using FGM L2PRE data files, named in the structure FILES, in the
%     time interval [TSTART, TEND), read data and upsample to SCM times,
%     given by T_SCM.
%
% Parameters
%   FILES:          in, required, type=struct
%                   Required fields:
%                       'fgm'        - FGM file name
%   T_SCM:          in, required, type=int64 (cdf_time_tt2000)
%   TSTART:         in, required, type=char
%   TEND:           in, required, type=char
%
% Returns
%   TT2000          out, required, type=1xN INT64 (cdf_time_tt2000)
%   B_BCS           out, required, type=3xN float
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-11-30      Written by Matthew Argall
%
function [t_out, b_out, bs_out] = mms_fsm_l2plus_fgm_interp(t_fgm, b_fgm, t_scm, b_scm)

%------------------------------------%
% Convert to Seconds                 %
%------------------------------------%
	t0        = t_scm(1);
	t_fgm_ssm = MrCDF_epoch2ssm(t_fgm, t0);
	t_scm_ssm = MrCDF_epoch2ssm(t_scm,  t0);

%------------------------------------%
% Interpolate to SCM                 %
%------------------------------------%
	
	%
	% Note: The range change affects the sampling rate.
	%       To deal with it, I hvae set TOL=2, but really
	%       I should check the RANGE flag in the L2 data.
	%

	% Find SLOW segments
	[iseg, nseg] = MrIntervalsX(t_fgm_ssm, [], 2);

	% Allocate memory & prepare to loop
	n_out  = length(t_scm);
	t_out  = zeros(1, n_out, 'int64');
	b_out  = zeros(3, n_out);
	bs_out = zeros(3, n_out);
	istart = 1;
	iend   = 0;

	% Report number of segments
	if nseg > 1
		mrfprintf('logtext', 'FGM is split into %d segments.', nseg);
	end

	% Step through each segment
	%   - This loop is designed to loop over segments with
	%     different sampling rates.
	for ii = 1 : nseg
		% Current segment of data
		b_seg     = b_fgm(:,iseg(1,ii):iseg(2,ii));
		t_seg     = t_fgm(iseg(1,ii):iseg(2,ii));
		t_seg_ssm = t_fgm_ssm(iseg(1,ii):iseg(2,ii));

		% Find overlap with SCM
		[ifgm, iscm] = MrIntervalsXY(t_seg_ssm, t_scm_ssm, ...
		                             'Remove', true, 'Sync', true, ...
		                             'TolX', 2);

		% Report number of overlapping segments
		if size(ifgm,2) > 1
			mrfprintf('logtext', '%d overlapping segments.', size(iseg,2));
		end

		% Overlapping intervals
		%   - Loop over overlapping data segments in FGM and SCM.
		%     Throw outand points in FGM that are not in SCM and
		%     vice versa
		%   - Interpolate FGM up to SCM
		for jj = 1 : size(ifgm, 2)
			% Extract the interval to be interpolated
			b_in   = b_seg(:,ifgm(1,jj):ifgm(2,jj));
			t_in   = t_seg_ssm(ifgm(1,jj):ifgm(2,jj));
			t_terp = t_scm_ssm(iscm(1,jj):iscm(2,jj));

			% Interpolate
			b_out = zeros(3, length(t_terp));
			b_out(1,:) = interp1(t_in, b_in(1,:), t_terp, 'linear', 'extrap');
			b_out(2,:) = interp1(t_in, b_in(2,:), t_terp, 'linear', 'extrap');
			b_out(3,:) = interp1(t_in, b_in(3,:), t_terp, 'linear', 'extrap');
		
			% Store the data
			istart = iend + 1;
			iend   = istart + length(t_terp) - 1;
			t_out(istart:iend)     = t_scm(iscm(1,jj):iscm(2,jj));
			b_out(:,istart:iend)   = b_out;
			bs_out(:, istart:iend) = b_scm(:, iscm(1,jj):iscm(2,jj));
		end
	end
	
	% Change in total number of points
	mrfprintf('logwarn', '   Number of points reduced by %d from %d to %d.', n_out-iend, n_out, iend);
	
	% Trim the data
	t_out  = t_out(1:iend);
	b_out  = b_out(:,1:iend);
	bs_out = bs_out(:,1:iend);
end


%
% Name
%   mms_fsm_l2plus_scm_prep
%
% Purpose
%   Read FGM L1A data and calibrate it, creating a data product in OMB coordinates.
%
% Calling Sequence
%   [TT2000_FGM, B_FGM_OMB, SR_FGM] = mms_fsm_ql_prep_scm(FILES, FGM_INSTR, TSTART, TEND);
%     Using data files in the structure FILES from instrument FGM_INSTR between the
%     time interval of [TSTART, TEND), read and calibrate data, producing time stamps
%     TT2000_FGM, magnetic field in OMB B_FGM_OMB, and sample rate SR_FGM.
%
% Parameters
%   FILES:          in, required, type=struct
%                   Required fields:
%                       'scm'      - SCM file name
%                       'scm_cal'  - Sensor temperature HK10E file
%   TSTART:         in, required, type=char
%   TEND:           in, required, type=char
%
% Returns
%   TT2000_SCM      out, required, type=1xN INT64 (cdf_time_tt2000)
%   B_SCM_OMB       out, required, type=3xN float
%   XFR_SCM         out, reuqired, type=3xN 
%   SR_SCM          out, reuqired, type=1xN 
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-10-27      Written by Matthew Argall
%
function [tt2000, b_omb, xfr, sr_med] = mms_fsm_l2plus_scm_prep(file, cal_file, duration)

	% Read uncalibrated SCM in 123
	scm_l1a = mms_scm_l1a_read(file);
	
	% Read calibration data
	[transfr_fn, freqs] = mms_scm_read_caltab(cal_file);

	% Convert numbers to nano-Tesla
	%   - Call this OMB. The SCM team considers 123 to be orthogonalized already.
	%   - Technically, data in OMB is fully calibrated, but the remainder of our
	%     calibration process will be performed simultaneously as the data is merged.
	%   - SCM is inverted with respect to AFG and DFG. Negate it.
	b_omb = -mms_scm_number2nT(scm_l1a.b_123);

	% Extract the time
	tt2000 = scm_l1a.tt2000;
	
	% Determine the sample rate
	%   - Integer sampling rate
	%   - 16=slow, 32=fast, 128=brst
	sr     = round( 1.0 ./ ( double( diff(tt2000) ) * 1e-9 ) );
	sr_med = median( sr );

	% Frequency resolution
	%   - Number of samples per calibration interval
	df = 1.0 / duration;
	N  = duration * sr_med;

	% Compensation function
	%   - Interpolate the transfer function to DF
	xfr = mms_scm_tf_compensate(transfr_fn, freqs, N, df);
end


%
% Name
%   mms_fsm_ql_prep_fgm
%
% Purpose
%   Read FGM L1A data and calibrate it, creating a data product in OMB coordinates.
%
% Calling Sequence
%   [TT2000, B_BCS] = mms_fsm_ql_prep_fgm(FILES, T_SCM, TSTART, TEND);
%     Using FGM L2PRE data files, named in the structure FILES, in the
%     time interval [TSTART, TEND), read data and upsample to SCM times,
%     given by T_SCM.
%
% Parameters
%   FILES:          in, required, type=struct
%                   Required fields:
%                       'fgm'        - FGM file name
%   T_SCM:          in, required, type=int64 (cdf_time_tt2000)
%   TSTART:         in, required, type=char
%   TEND:           in, required, type=char
%
% Returns
%   TT2000          out, required, type=1xN INT64 (cdf_time_tt2000)
%   B_BCS           out, required, type=3xN float
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-11-30      Written by Matthew Argall
%
function [v_out] = mms_fsm_l2pre_interp(v, t, t_out)

	% Convert time to seconds
	t0        = t(1);
	t_ssm     = MrCDF_epoch2ssm(t,     t0);
	t_out_ssm = MrCDF_epoch2ssm(t_out, t0);
	
	% Extrapolate before
	ibefore = find(t_out >= t(1), 1, 'first');
	if ibefore > 1
		mrfprintf('logwarn', 'Extrapolating %d points before.', ibefore-1);
	end
	
	% Extrapolate after
	iafter = find(t_out <= t(end), 1, 'last');
	if iafter < length(t_out)
		mrfprintf('logwarn', 'Extrapolating %d points after.', length(t_out)-iafter);
	end

	% Interpolate
	v_out      = zeros(3, length(t_out));
	v_out(1,:) = interp1(t_ssm, v(1,:), t_out_ssm, 'linear', 'extrap');
	v_out(2,:) = interp1(t_ssm, v(2,:), t_out_ssm, 'linear', 'extrap');
	v_out(3,:) = interp1(t_ssm, v(3,:), t_out_ssm, 'linear', 'extrap');
end


%
% Name
%   mms_fsm_ql_omb2smpa
%
% Purpose
%   Transform the merged magnetic field from OMB to DMPA coordinates.
%
% Calling Sequence
%   [B_MERGED_DMPA, B_FGM_DMPA] = mms_fsm_ql_omb2smpa(FILES, TSTART, TEND, T_MERGED, B_MERGED, T_FGM, B_FGM);
%     Using data files in the structure FILES between the time interval of [TSTART, TEND),
%     transform the data from OMB to SMPA, then despin SMPA into DMPA. Data have
%     TT2000 time stamps of T_MERGED and T_FGM, respectively.
%
% Parameters
%   FILES:          in, required, type=struct
%                   Required fields:
%                       'defatt'  - Definitive attitude files
%                       'dss'     - Sunpulse times from HK101
%   TSTART:         in, required, type=char
%   TEND:           in, required, type=char
%   B_MERGED:       in, required, type=3xN float
%   B_FGM_OMB:      in, required, type=3xN float
%
% Returns
%   B_MERGED_DMPA   out, required, type=1xN INT64 (cdf_time_tt2000)
%   B_FGM_DMPA      out, required, type=3xN float
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-10-27      Written by Matthew Argall
%
function [flag] = mms_fsm_l2pre_sr_flag(tt2000)

	% Find the sampling rate
	sr     = round( 1.0 / ( diff(tt2000) * 1e-9 ) );
	sr_med = median( sr );
	
	% Create the flag
	flag = sr == sr_med;
end