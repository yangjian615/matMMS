%
% Name
%   mms_fg_interp_stemp
%
% Purpose
%   Interpolate fluxgate sensor temperatures.
%
%   Method:
%       1. Find gaps larger than T_MAXGAP.
%       2. Boxcar smooth each data interval by T_NSMOOTH points.
%       3. Check distance to nearest neighbor.
%       4. Linearly interpolate temperature to T_OUT.
%       5. Find gaps in T_OUT larger than T_OUT_MAXGAP
%       6. Boxcar smooth each data interval by T_OUT_NSMOOTH.
%
% Calling Sequence
%   TEMP = mms_fg_interp_stemp(T_STEMP, STEMP, T_OUT)
%     Interpolate sensor temperatures STEMP measured at times T_STEMP to the
%     times given by T_OUT.
%
% Parameters
%   T_STEMP         in, required, type = int64 (cdf tt2000 times)
%   STEMP           in, required, type = N-element float
%   T_OUT           in, required, type = int64 (cdf tt2000 times)
%
% Returns
%   STEMP           out, required, type=float
%                   Interpolated sensor temperature. Has the same number of
%                       elements as T_OUT.
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-07-27      Written by Matthew Argall
%
function stemp = mms_fg_interp_stemp(t_stemp, stemp, t_out)

	% Temperature correction cadence
	%   - Configuration
	t_Nsmooth    = 13;         % number of points
	t_maxgap     = 4 * 60;     % seconds
	t_out_smooth = 50;         % seconds
	t_out_maxgap = 4 * 60;     % seconds

%------------------------------------%
% Find Gaps in Temperature           %
%------------------------------------%
	% Sampling interval
	dt_temp = diff( t_stemp );
	
	% Median sampling interval
	%   - Could always be 64 seconds.
	mdt_temp = median( dt_temp );
	
	% Find data gaps greater than T_MAXGAP
	%   - IGAPS   is the index before a gap
	%   - IGAPS+1 is the index after  a gap
	igaps    = find( double(dt_temp) / 1e9 > t_maxgap );
	ngaps    = length( igaps, 1 );
	
	% Warn about gaps
	if ngaps > 0
		warning('HK_10e:Smooth', '%d gaps > %f seconds in temperature data.', ngaps, t_maxgap)
	end
	
	% Ensure the first and last intervals are captured
	igaps = [0 igaps length(t_stemp) ];

%------------------------------------%
% Smooth Temperature                 %
%------------------------------------%
	% Step through each gap
	for ii = 1 : ngaps
		% Number of consecutive points
		istart = igaps[ii+1];
		istop  = igaps{ii] + 1;
		npts   = istop - istart + 1;
	
		%
		% Smooth the sensor temperature.
		%
		% Ken's IDL program uses Smooth(/EDGE_TRUNCTATE), which handles endpoints
		% differently from MATLAB's Smooth() function
		%   - http://exelisvis.com/docs/SMOOTH.html
		%   - http://www.mathworks.com/help/curvefit/smooth.html?searchHighlight=smooth
		if npts > 1
			nsmooth             = min( [t_smooth_param, npts] );
			stemp(istart:istop) = smooth( stemp(istart:istop), nsmooth );
		end
	end

%------------------------------------%
% Distance to Nearest Neighbor       %
%------------------------------------%
	%
	% Determine how far each data point (given by T_OUT) is from a 
	% true temperature sample (given by T_STEMP).
	%
	
	% Locate T_OUT within T_STEMP
	iloc = MrValue_Locate(t_stemp, t_out);
	
	% Find how far away the lower/higher nearest neighbor is
	tf_up   = iloc == 0;
	tf_down = iloc+1 > length(t_stemp);
	nnlow   = t_stemp(iloc   + tf_up)   - t_out;
	nnhigh  = t_stemp(iloc+1 - tf_down) - t_out;
	
	% Find the smallest distance between neighbors
	%   - Convert from nanoseconds to seconds.
	t_gap = min( abs(nnlow), abs(nnhigh) ) / 1e9;
	
	% Are there large differences between T_OUT and T_STEMP?
	igap = find( t_gap > t_maxgap );
	ngap = length( igap );
	
	% Warn about large gaps
	if ngap > 0
		warning('HK_10e:Interpolate', '%d points inter/extrapolated %f seconds from nearest temperature reading.\n'
		                              'Maximum difference is %f seconds.', ngap, t_maxgap, max(t_gap(:)));
	end

%------------------------------------%
% Interpolate Temperature            %
%------------------------------------%
	% Linearly interpolate data
	stemp = interpol(t_stemp, stemp, t_out, 'linear');

%------------------------------------%
% Smooth Between Gaps in Output Data %
%------------------------------------%
	% Smooth and gap parameters
	t_smooth = 50;       % Seconds
	t_maxgap = 4 * 60;   % Seconds
	
	% Data gaps in interpolation times
	dt_out = diff(t_out);
	igap   = find( double(dt_out) / 1e9 > t_out_maxgap );
	ngap   = length( igap );
	
	% Warn about gaps
	if ngaps > 0
		warning('HK_10e:Smooth', '%d gaps > %f seconds in output data.', ngaps, t_out_maxgap)
	end
	
	% Ensure the first and last intervals are captured
	igaps = [0 igaps length(t_out) ];
	
	% Step through each data interval
	for ii = 1 : ngaps
		% Number of consecutive points
		istart = igaps[ii] + 1;
		istop  = igaps[ii+1];
		npts   = istop - istart + 1;
	
		% Smooth
		if npts > 1
			% Number of points to smooth
			%   - Make sure there are enough points in the interval
			dt_med  = median( double( dt_out(istart:istop) ) / 1e9 );
			nsmooth = min( [ round(t_out_smooth / dt_med), npts ] );
			
			% Smooth -- STEMP now has the same number of points as T_OUT
			stemp(istart:istop) = smooth( stemp(istart:istop), nsmooth );
		end
	end
end