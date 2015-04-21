%
% Name
%   mms_edi_bavg
%
% Purpose
%   Compute the average magnetic field. This field is used to define the plane
%   perpendicular to the magnetic field in which EDI beams are fired.
%
% Calling Sequence
%   B_DATA = mms_edi_bavg(T_FG, B_FG, T_GD12, T_GD21)
%     Create an average magnetic field using fluxgate magnetometer data B_FG
%     and its time tags T_FG, by first interpolating B_FG onto EDI beam times
%     from gun-detector pair 12, T_GD12, and 21, T_GD21, then by averagin
%     the field within 5 second intervals.
%
%   B_DATA = mms_edi_read_efield(__, DT)
%     Specify the time interval, in seconds, ing which the magnetic field is
%     averaged.
%
% Parameters
%   T_FG            in, required, type=char/cell
%   B_FG            in, required, type=char
%   T_GD12          in, required, type=char
%   T_GD21          in, required, type=char
%   DT              in, optional, type=char, default=5.0
%
% Returns
%   B_DATA          out, required, type=structure
%                   Fields are:
%                     't_avg'     -  Time of the averaged data.
%                     'dt_avg'    -  Time interval between points.
%                     'b_avg'     -  Averaged data.
%                     'b_std'     -  Standard deviation of the averaged data.
%                     'b_gd12'    -  B_FG interpolated onto T_GD12.
%                     'b_gd21'    -  B_FG interpolated onto T_GD21.
%                     'inds_gd12' -  Indices into B_AVG onto which T_GD12 map.
%                     'inds_gd21' -  Indices into B_AVG onto which T_GD21 map.
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-20      Written by Matthew Argall
%
function b_avg = mms_edi_bavg(t_fg, b_fg, t_gd12, t_gd21, dt)

	if nargin < 5
		dt = 5.0;
	end

	% Times will be based on DFG epoch.
	%   - If we do not have a magnetic field vector for a 5sec interval,
	%     there is no way to determine beam orientation.
	tvec = MrCDF_Epoch_Breakdown( [t_fg(1), t_fg(end)] );
	
	% Round down to the nearest 5 seconds and recompute
	tvec(:,6)     = tvec(:,6) - mod(tvec(:,6), 5);
	tvec(:,7:end) = 0.0;
	tedge         = MrCDF_Epoch_Compute(tvec);
	t0            = tedge(1);
	t1            = tedge(2) + int64(dt * 1d9);
	
	% Number of points
	N_fg   = length(t_fg);
	N_gd12 = length(t_gd12);
	N_gd21 = length(t_gd21);
	
	% Allocate memory
	nOut      = max( [N_fg N_gd12 N_gd21] );
	t_avg     = zeros(1, nOut, 'int64');
	b_avg     = zeros(3, nOut);
	b_std     = zeros(3, nOut);
	inds_gd12 = zeros(1, N_gd12, 'int32');
	inds_gd21 = zeros(1, N_gd21, 'int32');
	b_gd12    = zeros(3, N_gd12);
	b_gd21    = zeros(3, N_gd21);
	
%------------------------------------%
% Loop Through Each Interval         %
%------------------------------------%

	% Interpolation requires floats. Convert to seconds
	t_fg_sse   = MrCDF_epoch2sse(t_fg,   t0);
	t_gd12_sse = MrCDF_epoch2sse(t_gd12, t0);
	t_gd21_sse = MrCDF_epoch2sse(t_gd21, t0);

	% Initial conditions
	%    - TT:             Start at the first magnetometer data point.
	%    - DT_INTERP:      Number of seconds on either side of the averaging interval that are interpolated
	%    - TRANGE_AVG:     First time interval to be averaged.
	%    - TRANGE_INTEPR:  First time interval to be interpolated.
	%    - N_BEAM[12]_TOT: Total number of beams included throughout the averaging process.
	tt = t0;
	dt_interp      = dt;
	t_avg_range    = [0, dt];
	t_interp_range = [-dt_interp, dt + dt_interp];
	count          = 1;
	n_tot_gd12     = 0;
	n_tot_gd21     = 0;
	
	% Begin loop
	while t_avg_range(1) < MrCDF_epoch2sse(t1, t0)
	%------------------------------------%
	% Find the Interpolation Interval    %
	%------------------------------------%
		% Interpolation window
		%    it_interp  = The magnetic field is interpolated onto the beam times. Cubic
		%                 Splines interpolation will be used, so to avoid interpolaing
		%                 at edges, we will interpolate over a larger interval than the
		%                 averaging interval, then extract the averaging interval.
		it_fg_interp   = find( t_fg_sse   >= t_interp_range(1) & t_fg_sse   < t_interp_range(2) );
		it_gd12_interp = find( t_gd12_sse >= t_interp_range(1) & t_gd12_sse < t_interp_range(2) );
		it_gd21_interp = find( t_gd21_sse >= t_interp_range(1) & t_gd21_sse < t_interp_range(2) );
	
		% No data
		if isempty(it_fg_interp) || ( isempty(it_gd12_interp) & isempty(it_gd21_interp) )
			epoch_range = MrCDF_sse2epoch(t_interp_range, t0);
			epoch_range = MrCDF_Epoch_Encode(epoch_range);
			warning('EDI:Bavg', 'No beams found between %s and %s', epoch_range{1}, epoch_range{2});
			
			% Move to next interval
			t_avg_range    = t_avg_range + dt;
			t_interp_range = t_interp_range + dt;
			continue
		end
		
		% Extract the interval for ease of access
		b_fg_temp   = b_fg(:, it_fg_interp);
		t_fg_temp   = t_fg_sse(:, it_fg_interp);
		t_gd12_temp = t_gd12_sse(it_gd12_interp);
		t_gd21_temp = t_gd21_sse(it_gd21_interp);
	
	%------------------------------------%
	% Interpolate                        %
	%------------------------------------%
		% Interpolate to GD12 pair
		b_fg_gd12      = zeros(3, length(t_gd12_temp));
		b_fg_gd12(1,:) = spline( t_fg_temp, b_fg_temp(1,:), t_gd12_temp );
		b_fg_gd12(2,:) = spline( t_fg_temp, b_fg_temp(2,:), t_gd12_temp );
		b_fg_gd12(3,:) = spline( t_fg_temp, b_fg_temp(3,:), t_gd12_temp );
		
		% Interpolate to GD21 pair
		b_fg_gd21      = zeros(3, length(t_gd21_temp));
		b_fg_gd21(1,:) = spline( t_fg_temp, b_fg_temp(1,:), t_gd21_temp );
		b_fg_gd21(2,:) = spline( t_fg_temp, b_fg_temp(2,:), t_gd21_temp );
		b_fg_gd21(3,:) = spline( t_fg_temp, b_fg_temp(3,:), t_gd21_temp );
	
	%------------------------------------%
	% Find Avg Interval & Average        %
	%------------------------------------%
		% Get the beams within the averaging interval
		%   --|-----{-----}-----|--
		%     |  Interpolation  |
		%           { Avg }
		it_gd12_avg = find( t_gd12_temp >= t_avg_range(1) & t_gd12_temp < t_avg_range(2) );
		it_gd21_avg = find( t_gd21_temp >= t_avg_range(1) & t_gd21_temp < t_avg_range(2) );
		
		% Extract the data
		b_fg_gd12 = b_fg_gd12(:, it_gd12_avg);
		b_fg_gd21 = b_fg_gd21(:, it_gd21_avg);
		
		% Take the average
		t_avg(count)   = int64(t_avg_range(1) * 1d9) + t0;
		b_avg(:,count) = mean( [b_fg_gd12 b_fg_gd21], 2 );
		b_std(:,count) = std( [b_fg_gd12 b_fg_gd21], 0, 2 );
	
		% Number of beams in the average
		n_gd12 = length(it_gd12_avg);
		n_gd21 = length(it_gd21_avg);
	
	%------------------------------------%
	% Store the Interpolated Data        %
	%------------------------------------%
		% Magnetic field
		b_gd12(:, n_tot_gd12+1:n_tot_gd12+n_gd12) = b_fg_gd12;
		b_gd21(:, n_tot_gd21+1:n_tot_gd21+n_gd21) = b_fg_gd21;
		
		% Index of B associated with each beam
		inds_gd12(n_tot_gd12+1:n_tot_gd12+n_gd12) = count;
		inds_gd21(n_tot_gd21+1:n_tot_gd21+n_gd21) = count;

	%------------------------------------%
	% Next Interval                      %
	%------------------------------------%
		
		% Up the number of b-fields
		count = count + 1;
		
		% Total number of beams
		n_tot_gd12 = n_tot_gd12 + n_gd12;
		n_tot_gd21 = n_tot_gd21 + n_gd21;
		
		% Next interval
		t_avg_range    = t_avg_range + dt;
		t_interp_range = t_interp_range + dt;
	end
	
%------------------------------------%
% Trim Data                          %
%------------------------------------%

	% Trim data
	t_avg = t_avg(1:count-1);
	b_avg = b_avg(:, 1:count-1);
	b_std = b_std(:, 1:count-1);
	
	% Output structure
	b_avg = struct( 't_avg',     t_avg,     ...
	                'dt_avg',    dt,        ...
	                'b_avg',     b_avg,     ...
	                'b_std',     b_std,     ...
	                'b_gd12',    b_gd12,    ...
	                'b_gd21',    b_gd21,    ...
	                'inds_gd12', inds_gd12, ...
	                'inds_gd21', inds_gd21 );
end