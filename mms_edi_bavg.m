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
%                     't_avg'       -  Time of the averaged data.
%                     'dt_avg'      -  Time interval between points.
%                     'b_avg'       -  Averaged data.
%                     'b_std'       -  Standard deviation of the averaged data.
%                     'b_gd12'      -  B_FG interpolated onto T_GD12.
%                     'b_gd21'      -  B_FG interpolated onto T_GD21.
%                     'recnum'      -  Record number of each element of T_AVG
%                     'recnum_gd12' -  Location of each GD12 beam within T_AVG.
%                     'recnum_gd21' -  Location of each GD21 beam within T_AVG.
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-20      Written by Matthew Argall
%   2015-06-19      GD12 and GD21 are handled separately. Can have beams
%                       from only one of the guns. - MRA
%   2015-06-21      Indices are now 0-based instead of 1-based. - MRA
%   2015-06-22      Added RECNUM fields to output structure. - MRA
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
	nOut        = max( [N_fg N_gd12 N_gd21] );
	t_avg       = zeros(1, nOut, 'int64');
	b_avg       = zeros(3, nOut);
	b_std       = zeros(3, nOut);
	recnum      = zeros(1, nOut,   'uint32');
	recnum_gd12 = zeros(1, N_gd12, 'uint32');
	recnum_gd21 = zeros(1, N_gd21, 'uint32');
	b_gd12      = zeros(3, N_gd12);
	b_gd21      = zeros(3, N_gd21);
	
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
	
		% Extract the interval for ease of access
		b_fg_temp   = b_fg(:, it_fg_interp);
		t_fg_temp   = t_fg_sse(:, it_fg_interp);
		
	%------------------------------------%
	% GD12                               %
	%------------------------------------%
		if isempty(it_gd12_interp)
			n_gd12 = 0;
			trange_interp = MrCDF_Epoch_Encode( MrCDF_sse2epoch( t_interp_range, t0 ) );
			warning('EDI:Bavg', 'No GD12 beams to interpolate: %s and %s', trange_interp{1}, trange_interp{2});
		else
			% Extract interpolation interval
			t_gd12_temp = t_gd12_sse(it_gd12_interp);
			
			% Allocate memory to output array
			b_fg_gd12 = zeros(3, length(t_gd12_temp));

			% Interpolate (x,y,z) components
			b_fg_gd12(1,:) = spline( t_fg_temp, b_fg_temp(1,:), t_gd12_temp );
			b_fg_gd12(2,:) = spline( t_fg_temp, b_fg_temp(2,:), t_gd12_temp );
			b_fg_gd12(3,:) = spline( t_fg_temp, b_fg_temp(3,:), t_gd12_temp );
			
			% Get the beams within the averaging interval
			%  --|-----{-----}-----|--
			%    |  Interpolation  |
			%          { Avg }
			it_gd12_avg = find( t_gd12_temp >= t_avg_range(1) & t_gd12_temp < t_avg_range(2) );
			if isempty(it_gd12_avg)
				n_gd12 = 0;
				trange_avg = MrCDF_Epoch_Encode( MrCDF_sse2epoch( t_avg_range, t0 ) );
				warning('EDI:Bavg', 'No GD12 beams to average: %s and %s', trange_avg{1}, trange_avg{2});
			else
				% Magnetic field in averaging interval
				b_fg_gd12 = b_fg_gd12(:, it_gd12_avg);
				
				% Number of beams in the average
				n_gd12 = length(it_gd12_avg);
			end
		end
		
	%------------------------------------%
	% GD21                               %
	%------------------------------------%
		if isempty(it_gd21_interp)
			n_gd21 = 0;
			trange_interp = MrCDF_Epoch_Encode( MrCDF_sse2epoch( t_interp_range, t0 ) );
			warning('EDI:Bavg', 'No GD21 beams to interpolate: %s and %s', trange_interp{1}, trange_interp{2});
		else
			% Extract interpolation interval
			t_gd21_temp = t_gd21_sse(it_gd21_interp);
			
			% Allocate memory to output array
			b_fg_gd21 = zeros(3, length(t_gd21_temp));
			
			% Interpolate (x,y,z) components
			b_fg_gd21(1,:) = spline( t_fg_temp, b_fg_temp(1,:), t_gd21_temp );
			b_fg_gd21(2,:) = spline( t_fg_temp, b_fg_temp(2,:), t_gd21_temp );
			b_fg_gd21(3,:) = spline( t_fg_temp, b_fg_temp(3,:), t_gd21_temp );
			
			% Get the beams within the averaging interval
			%  --|-----{-----}-----|--
			%    |  Interpolation  |
			%          { Avg }
			it_gd21_avg = find( t_gd21_temp >= t_avg_range(1) & t_gd21_temp < t_avg_range(2) );
			if isempty(it_gd21_avg)
				n_gd21 = 0;
				trange_avg = MrCDF_Epoch_Encode( MrCDF_sse2epoch( t_avg_range, t0 ) );
				warning('EDI:Bavg', 'No GD21 beams to average: %s and %s', trange_avg{1}, trange_avg{2});
			else
				% Magnetic field in averaging interval
				b_fg_gd21 = b_fg_gd21(:, it_gd21_avg);
				
				% Number of beams in the average
				n_gd21 = length(it_gd21_avg);
			end
		end
	
	%------------------------------------%
	% Find Avg Interval & Average        %
	%------------------------------------%
		switch 1
			case n_gd12 > 0 && n_gd21 > 0
				% Time stamp comes at the beginning of the interval
				t_avg(count) = int64(t_avg_range(1) * 1d9) + t0;
			
				% Average field & standard deviation
				b_avg(:,count) = mean( [b_fg_gd12 b_fg_gd21], 2 );
				b_std(:,count) = std(  [b_fg_gd12 b_fg_gd21], 0, 2 );
				recnum(count)  = count - 1;
			case n_gd12 > 0
				t_avg(count)   = int64(t_avg_range(1) * 1d9) + t0;
				b_avg(:,count) = mean( b_fg_gd12, 2 );
				b_std(:,count) = std(  b_fg_gd12, 0, 2 );
				recnum(count)  = count - 1;
			case n_gd21 > 0
				t_avg(count)   = int64(t_avg_range(1) * 1d9) + t0;
				b_avg(:,count) = mean( b_fg_gd21, 2 );
				b_std(:,count) = std(  b_fg_gd21, 0, 2 );
				recnum(count)  = count - 1;
			otherwise
				% Do nothing
				%  - This ensures that the Electric field, drift velocity,
				%    and drift step calculated hereafter will have the same
				%    number of points which occur at the same times as the
				%    averaged magnetic field
		end
	
	%------------------------------------%
	% Store the Interpolated Data        %
	%------------------------------------%
		% Magnetic field
		if n_gd12 > 0
			% Interpolated field
			b_gd12(:, n_tot_gd12+1:n_tot_gd12+n_gd12) = b_fg_gd12;
			
			% Index (0-based) into B_avg
			recnum_gd12(n_tot_gd12+1:n_tot_gd12+n_gd12) = count - 1;
		end
		
		if n_gd21 > 0
			% Interpolated field
			b_gd21(:, n_tot_gd21+1:n_tot_gd21+n_gd21) = b_fg_gd21;
			
			% Index (0-based) into B_avg
			recnum_gd21(n_tot_gd21+1:n_tot_gd21+n_gd21) = count - 1;
		end

	%------------------------------------%
	% Next Interval                      %
	%------------------------------------%
		if n_gd12 > 0 || n_gd21 > 0
			count = count + 1;
		end
		
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
	t_avg  = t_avg(1:count-1);
	b_avg  = b_avg(:, 1:count-1);
	b_std  = b_std(:, 1:count-1);
	recnum = recnum(1:count-1);
	
	% Output structure
	b_avg = struct( 't_avg',       t_avg,       ...
	                'dt_avg',      dt,          ...
	                'b_avg',       b_avg,       ...
	                'b_std',       b_std,       ...
	                'b_gd12',      b_gd12,      ...
	                'b_gd21',      b_gd21,      ...
	                'recnum',      recnum,      ...
	                'recnum_gd12', recnum_gd12, ...
	                'recnum_gd21', recnum_gd21 );
end