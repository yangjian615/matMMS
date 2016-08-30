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
%   BKGD = mms_fsm_bkgd_compute(DATA, T)
%     Use the data structures FGM and SCM containing time and magnetic field
%     wave forms to compute the amplitude, power, and phase over a sliding
%     window of duration T seconds. Append the results to the input structures.
%
% Parameters
%   DATA            in, required, type = struct
%                   Contains fields:
%                       't' - Time (TT2000)
%                       'b' - 3-component magnetic field
%   T               in, required, type = double
%
% Returns
%   BKGD            out, required, type=struct
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
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2016-08-14      Written by Matthew Argall
%
function bkgd = mms_fsm_bkgd_compute_one(data, T)

%------------------------------------%
% Find Unique Flag Values            %
%------------------------------------%
	% Unique Flags
	[flag, iFlag] = unique( data.flag );

	% Sort flags by time
	[~, isort] = sort( data.t(iFlag) );
	flag       = flag(isort);

%------------------------------------%
% Cycle Throught Each FGM Flag       %
%------------------------------------%
	for ii = 1 : length(flag)
		% Indices associated with the current flag
		idx = find( data.flag == flag(ii) );

	%------------------------------------%
	% Avoid Data Gaps                    %
	%------------------------------------%
		
		% Search for continuous data intervals
		[iInt, nInt] = MrIntervalsX( data.t(idx) );

	%------------------------------------%
	% Loop Through Data Intervals        %
	%------------------------------------%
		% Number of data intervals found
		temp   = struct([]);
		for jj = 1 : nInt
			% Subinterval to be processed
			i0 = idx(1) + iInt(1,jj) - 1;
			i1 = idx(1) + iInt(2,jj) - 1;
			
			% Convert time to seconds since midnight
			%   - Doubles needed for interpolation
			t_ref = data.t( i0 );
			t_ssm = MrCDF_epoch2sse( data.t( i0:i1 ), t_ref);
	
		%------------------------------------%
		% Compute Spectral Components        %
		%------------------------------------%
			
			% Compute amplitude, phase, power
			temp_jj = mms_fsm_bkgd_spectra(t_ssm, data.b(:,i0:i1), T);
	
		%------------------------------------%
		% Histogram Spectral Components      %
		%------------------------------------%
			
			% Enough data?
			%   - FGM and SCM should have identical number of points and time stamps
			if ~isempty(temp_jj.t)
				% Histogram spectral components
				temp_jj = mms_fsm_bkgd_histogram(temp_jj);

				% Finalize structure
				%  - Convert time back to TT2000
				%  - Include flag
				temp_jj.t             = MrCDF_sse2epoch(temp_jj.t, t_ref);
				temp_jj.('flag')      = repmat( flag(ii), 1, length(temp_jj.t) );
				temp_jj.('hist_flag') = flag(ii);
				
				% Sum Histogram Occurrence over Multiple Intervals
				if isempty(temp)
					temp = temp_jj;
				else
					temp.('amp_hist')   = temp.('amp_hist')   + temp_jj.('amp_hist');
					temp.('phase_hist') = temp.('phase_hist') + temp_jj.('phase_hist');
					temp.('psd_hist')   = temp.('psd_hist')   + temp_jj.('psd_hist');
				end
			end
			
			clear temp_jj
		end
		
	%------------------------------------%
	% Determine Noise Floor              %
	%------------------------------------%
		% Amplitude, Power, Phase
		temp.('amp_floor')   = mms_fsm_bkgd_noisefloor( temp.('amp_hist'),   temp.('amp_bins')   );
		temp.('phase_floor') = mms_fsm_bkgd_noisefloor( temp.('phase_hist'), temp.('phase_bins') );
		temp.('psd_floor')   = mms_fsm_bkgd_noisefloor( temp.('psd_hist'),   temp.('psd_bins')   );
	
	%------------------------------------%
	% Save Data from Multiple Flags      %
	%------------------------------------%
		if ii == 1
			bkgd = temp;
		else
			bkgd = [bkgd temp];
		end
		
		clear temp
	end

%------------------------------------%
% Concatenate Multiple Tags          %
%------------------------------------%
	% Combine different operational modes
	if length(bkgd) > 1
		% Ensure the length of the frequency dimension is the same for slow and fast survey
		nSlow = length( bkgd(1).f );
		nFast = length( bkgd(end).f );
		if nSlow ~= nFast

			% If the are not, extend the slow survey data
			names = fieldnames(bkgd);
			for ii = 1 : length(bkgd) - 1
				for jj = 1 : length( names )
					if ismember(names{jj}, {'amp', 'phase', 'psd', 'amp_hist', 'phase_hist', 'psd_hist', 'amp_floor', 'phase_floor', 'psd_floor'})
						% Array extension size (frequency is the last dimension)
						dims      = size( bkgd(ii).(names{jj}) );
						dims(end) = nFast - nSlow;
					
						% Concatenate NaNs to fill space
						bkgd(ii).(names{jj}) = cat( length(dims), bkgd(ii).(names{jj}), NaN(dims) );
					end
				end
			
				% Copy the frequency array
				bkgd(ii).f = bkgd(end).f;
			end
		end

		% Concatenate the fields
		bkgd = MrCatStruct(bkgd, [2,0,2,2,2,2,0,2,0,2,0,2,2,2,2,2]);
	end
end