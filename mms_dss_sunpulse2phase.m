%
% Name
%   mms_hk_read_sunpulse
%
% Purpose
%   Read houskeeping 0x101 sun pulse and related data.
%
% Calling Sequence
%   SUNPULSE = mms_hk_sunpulse(FILENAME)
%     Read housekeeping 0x101 packet times, 'Epoch', time of last sun
%     pulse in tt2000, 'SunPulse', sun pulse flag, 'Flag', and external sun
%     pulse period in microseconds, 'Period'. Return a structure SUNPULSE
%     with those as fields.
%
%   [__] = mms_hk_sunpulse(__, 'ParamName', ParamValue)
%     Any parameter name-value pair given below.
%
% Parameters
%   FILENAME        in, required, type=char
%   'UniqPackets'   in, optional, type=boolean, default=false
%                   If set, information from unique packets is returned.
%                     Packets can can contain overlapping data.
%   'UniqPulse'     in, optional, type=boolean, default=false
%                   Return unique sun pulse times. Sun pulse times are
%                     reported multiple times per packet, even if unchanged.
%                     Setting this parameter true automatically sets
%                     'UniqPackets' to true.
%
% Returns
%   SUNPULSE        out, required, type=struct
%                   Fields are:
%                     'Epoch'    - Packet times
%                     'SunPulse' - Sun pulse times
%                     'Period'   - Period (micro-sec) of revolution Only
%                                    returned when FLAG=0 and only on the
%                                    second and subsequent received sun
%                                    pulses from the s/c.
%                     'Flag'     - Status flag
%                                    0: s/c sun pulse
%                                    1: s/c pseudo sun pulse
%                                    2: s/c CIDP generated speudo sun pulse
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products
%   Signal Processing Toolbox  --  medfilt1.m
%
% History:
%   2015-03-25      Written by Matthew Argall
%   2015-04-12      Added the 'UniqPackets' and 'UniqPulse' parameters. - MRA
%   2015-04-14      Return a structure. - MRAs
%   2015-04-27      Corrected check of extrapolation at end of array. - MRAs
%
function phase = mms_dss_sunpulse2phase(sunpulse, times)

	% Extract for ease of access
	period   = int64( sunpulse.('Period') ) * int64(1e3);
	hk_epoch = sunpulse.('Epoch');
	pulse    = sunpulse.('SunPulse');
	flag     = sunpulse.('Flag');
	nPts     = length(period);
	
	% Check if there are sunpulse times
	assert(~isempty(pulse), 'SunPulse array is empty. Cannot determine phase.');

	% A period is valid if FLAG = 0 & 2 < PERIOD < 50 seconds.
	valid_period = flag == 0 & period < 50e9 & period > 2e9;
	
	% Do not count the first period as valid
	%   - PULSE is recorded every sun pulse, so the difference between
	%     pulse times should be the spin period. This should be comparable
	%     to PERIOD.
	dPulse      = double( diff( pulse ) );
	dPulse_flag = zeros(1, nPts, 'int8');
	
	% If the first period is a valid value, it can be used to create
	% an epoch time just prior to the end of a gap.
	if valid_period(1)
		pulse        = [ pulse(1) - period(1) pulse ];
		dPulse       = [ period(1) dPulse ];
		dPulse_flag  = [ 0 dPulse_flag ];
		period       = [ 0 period ];
		valid_period = [ 0 valid_period ];
		flag         = [ 3 flag ];
	end
	
%------------------------------------%
% Intervals and Gaps                 %
%------------------------------------%
	% Find data gaps by examining the period
	T_median        = median( double( dPulse ) );
	[iGaps, nInGap] = MrGapsX( pulse, T_median);
	
	% Fill in gaps that are > NMINGAP periods. Ignore others
	nMinGap  = 4;
	iFill    = find( nInGap > nMinGap );
	nBigGaps = length(iFill);
	
	% Warn about gaps
	if nBigGaps > 0
		msg = sprintf('%d gaps of more than %f seconds.', nBigGaps, double(nMinGap*T_median)*1e-9);
		warning('SunPulse:Phase', msg);
	end
	
	% Insert a pseudo sun pulse before each big gap
	%   - The first period is a valid value, so can be used to create
	%     an epoch time just prior to the end of a gap.
	for ii = 1 : nBigGaps
		% Index at which to add a pseudo pulse
		%   - We already took care of the first interval
		%   - Do so only for the gaps larger than NBIGGAPS periods
		%   - IGAP is the index of the gap of interest.
		%   - We want to add a point before the gap ends.
		igap = iFill(ii);
		idx  = iGaps(2, igap);
	
		% Insert pseudo sunpulse only if the reported period is valid
		if valid_period(idx)
			speudopulse  = pulse(idx) - period(idx);
			pulse        = [ pulse(1:idx-1)         speudopulse  pulse(idx:end)        ];
			period       = [ period(1:idx-1)        period(idx)  period(idx:end)       ];
			flag         = [ flag(1:idx-1)          3            flag(idx:end)         ];
			dPulse       = [ dPulse(1:idx-1)        period(idx)  dPulse(idx:end)       ];
			dPulse_flag  = [ dPulse_flag(1:idx-1)   0            dPulse_flag(idx:end)  ];
			valid_period = [ valid_period(1:idx-1)  0            valid_period(idx:end) ];
	
			% Bump the gap indices forward one
			iGaps( igap+1:end ) = iGaps( igap+1:end ) + 1;
		end
	end
	
%------------------------------------%
% Smooth Data Intervals              %
%------------------------------------%

	% medfilt1 requires double precision
	dPulse = double( dPulse );
	
	% Process each continuous intervals
	%   - There is one more data interval than there are big gaps
	nMedFilt = 7;
	halfWin  = fix( nMedFilt/2 );
	istart   = 1;
	for ii = 1 : nBigGaps + 1
		if ii == nBigGaps + 1
			igap = [];
			iend = length(dPulse);
		else
			igap = iFill(ii);
			iend = iGaps(1, igap);
		end

		% Median smooth if we can
		if (iend - istart) >= nMedFilt
			% Median smooth the period during this interval
			medT = medfilt1( dPulse(istart:iend), nMedFilt );

			% Care for the edges
			medT( 1:halfWin )       = medT( halfWin+1 );
			medT( end-halfWin:end ) = medT( end-halfWin-1 );
		else
			medT    = zeros(1, iend - istart + 1);
			medT(:) = median( dPulse(istart:iend) );
		end

		% Number of spins between each pulse
		%   - Recall, T is the time between sun pulses, not the reported period
		%   - There could still be gaps smaller than NMINGAP periods
		nSpins = medT / T_median;

		% NSPINS should be roughly an integer. If not, the sun pulse is
		% changing in an unexpected manner.
		iNotIntSpin = find( abs( nSpins - round(nSpins) ) >= 0.25 );
		if ~isempty(iNotIntSpin)
			% Send a warning
			msg = sprintf( '%d periods that are non-integer multiples of the median spin period.', ...
			               sum(iNotIntSpin));
			warning('SunPulse:Phase', msg);
			
			% Set the flag
			T_flag(iNotIntSpin) = 2;
		end
		
		%
		% Other things to consider:
		%   - Could fill gaps of length < NMINGAP by using the spin period
		%     if it is valid (i.e. = 0).
		%
		
		% Compute the average spin period between pulses.
		%   - If there NSPINS=3 spins between pulses, this will give
		%     the average of the tree periods.
		dPulse(istart:iend) = dPulse(istart:iend) / nSpins;
		
		% Move to the next interval.
		istart = iGaps(2, igap);
	end
	
%------------------------------------%
% Interpolate over Large Gaps        %
%------------------------------------%
	%
	% "Interpolation" involves using the same period for all
	% points in TIMES that fall within a data gap.
	%
	% Figure out the average period within a large data gap
	%   1. Determine the period before the gap
	%   2. Determine the period after the gap
	%   3. Compare the two and decide which period to use
	%
	
	for ii = 1 : nBigGaps
		% Index of the point just prior to a large data gap
		igap = iFill(ii);
		idx  = iGaps(1, igap);
		
		% Period just before a large data gap
		%   - Use the DSS period, if it is valid
		%   - Cannot determine the period before the data begins
		if valid_period(idx)
			T1 = period(idx);
		elseif idx > 1
			T1 = dPulse(idx);
		else
			t_warn = spdfencodett2000( pulse(idx) );
			msg    = sprintf('Cannot determine period before gap at %s', t_warn{1});
			warning('SunPulse:Phase', msg);
			T1 = nan;
		end
		
		% Period just after a large data gap
		%   - The DSS period will not be valid, because we added
		%     a pseudopulse above.
		%   - Make sure the next valid point is not also the beginning
		%     of another data gap. This would occur if there is a
		%     single point stranded between two gaps.
		if idx + 1 == iGaps(2, igap)
			t_warn = spdfencodett2000( pulse(idx) );
			msg    = sprintf('Cannot determine period after gap at %s', t_warn{1});
			warning('SunPulse:Phase', msg);
			T2 = nan;
		elseif idx > 1
			T2 = dPulse(idx+1);
		end
		
		% Determine the number of spins in the gap. Use the period
		% before and after the gap. They should be the same.
		nSpin1 = dPulse(idx+1) / T1;
		nSpin2 = dPulse(idx+1) / T2;
		
		% Check if they are the same
		if round(nSpin1) == round(nSpin2)
			nSpins           = nSpin1;
			dPulse_flag(idx) = 1;
			
			% Do they differ by more than half a spin?
			if abs(nSpin1 - nSpin2) >= 0.5
				t_warn = spdfencodett2000( pulse(idx) );
				msg    = sprintf('Spin rates that bound gap at %s differ by >= 0.5 spins.', ...
				                 t_warn{1});
				warning('SunPulse:Phase', msg);
				dPulse_flag(idx) = 3;
			end
			
		% Use the first spin period less than 10 degrees from expected
		%    - 360 * 0.25 * 0.111 = 10
		elseif mod( nSpin1, 1 ) < 0.25 * 0.111
			nSpins           = nSpin1;
			dPulse_flag(idx) = 4;
			
			% Warn
			t_warn = spdfencodett2000( pulse(idx) );
			msg    = sprintf('Using period at beginning of gap %s', t_warn{1});
			warning('SunPulse:Phase', msg);
		
		% Use second spin if it is an integer to within 10%
		elseif mod( nSpin2, 1 ) < 0.1
			nSpins           = nSpin2;
			dPulse_flag(idx) = 4;
			
			% Warn
			t_warn = spdfencodett2000( pulse(idx) );
			msg    = sprintf('Using period at end of gap %s', t_warn{1});
			warning('SunPulse:Phase', msg);
		
		% Use average spin
		else
			nSpins           = (nSpin1 + nSpin2) / 2.0;
			dPulse_flag(idx) = 5;
			
			% Warn
			t_warn = spdfencodett2000( pulse(idx) );
			msg    = sprintf('Using average period before and after gap %s', t_warn{1});
			warning('SunPulse:Phase', msg);
		end
		
		% Average period
		dPulse(idx) = dPulse(idx) / nSpins;
	end
	
%------------------------------------%
% Compuate the Spin Phase            %
%------------------------------------%
	% Determine where TIMES is located within the sun pulse array.
	sp_sse = MrCDF_epoch2sse( pulse );
	t_sse  = MrCDF_epoch2sse( times, pulse(1) );
	inds   = MrValue_Locate(sp_sse, t_sse, 'RoundUp', true);
	
	% Calculate the phase in degrees
	%   - degrees * seconds / ( seconds / spin )
	phase = 360.0 * double(times - pulse(inds)) ./ double(period(inds));
	phase = mod(phase, 360.0);
	
%------------------------------------%
% Did we Extrapolate?                %
%------------------------------------%
	
	% Were there times before the first sun pulse
	iBefore = find(times >= pulse(1), 1, 'first');
	if isempty(iBefore)
		iBefore = 0;
	else
		iBefore = iBefore - 1;
	end
	
	% Were there times after the last sun pulse?
	iAfter = find(times <= pulse(end), 1, 'last');
	if isempty(iAfter)
		iAfter = length(times) + 1;
	else
		iAfter = iAfter + 1;
	end
	
	% Warn if excessive extrapolation
	if iBefore > 0
		if 360.0 * (pulse(1) - times(1)) / period(1) > 3.0 * 360.0
			warning('SunPulse:Phase', 'Extrapolating more than 3 spins before first point.');
		end
	end

	% Warn if excessive extrapolation
	if iAfter < length(times)
		if 360.0 * (times(end) - pulse(end)) / period(1) > 3.0 * 360.0
			warning('SunPulse:Phase', 'Extrapolating more than 3 spins after last point.');
		end
	end
end