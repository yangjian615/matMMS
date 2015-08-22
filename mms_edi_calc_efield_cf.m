%
% Name
%   mms_edi_calc_efield_cf
%
% Purpose
%   Compute the electric field from EDI beam data using the
%   beam convergence method.
%
% Calling Sequence
%   EFIELD_BC = mms_edi_calc_efield_cf(EDI, B_AVG)
%     Compute the EDI electric field EFIELD_BC using the beam
%     convertence (BC) method. Input data, EDI, should be a
%     structure of the type returned by mms_edi_create_l2pre_efield.m
%     and should be inerpolated onto intervals of magnetic
%     field data B_AVG, which is a structure of the type returned
%     by mms_edi_bavg.m.
%
% Parameters
%   EDI:            in, required, type=structure
%                   Fields used:
%                      'tt2000_gd12'
%                      'tt2000_gd21'
%                      'fv_gd12_123'
%                      'fv_gd21_123'
%                      'gun1_pos_bcs'
%                      'gun1_pos_bcs'
%   FG_L1B:         in, required, type=structure
%                   Fields used:
%                      'tt2000'
%                      'b_bcs'
%
% Returns
%   EFIELD_BC       out, required, type=string
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-08-18      Written by Matthew Argall
%
function [efield_cf, b_avg_bcs] = mms_edi_calc_efield_cf( t_fg, b_fg_bcs, ...
                                                          t_edi1, pos_vg1_edi1, fv_gd12_edi1, ...
                                                          t_edi2, pos_vg2_edi2, fv_gd21_edi2, dt )

	if nargin < 3
		dt = 5.0;
	end
	method = 1;

%------------------------------------%
% Compuate Average B                 %
%------------------------------------%
	% Time range that we have data
	if isempty(t_edi1) && isempty(t_edi2)
		error( 'No EDI E-field data available.' );
	elseif isempty(t_edi1)
		t0 = t_edi2(1);
		t1 = t_edi2(end);
	elseif isempty(t_edi2)
		t0 = t_edi1(1);
		t1 = t_edi1(end);
	else
		t0 = min( [t_edi1(1)   t_edi2(1)  ] );
		t1 = max( [t_edi1(end) t_edi2(end)] );
	end
	
	% Breakdown into time vectors
	tvec = MrCDF_Epoch_Breakdown( [t0, t1] );
	
	% Round down to the nearest 5 seconds and recompute
	tvec(:,6)     = tvec(:,6) - mod(tvec(:,6), 5);
	tvec(:,7:end) = 0.0;
	tedge         = MrCDF_Epoch_Compute(tvec);
	t0            = tedge(1);
	t1            = tedge(2) + int64(dt * 1d9);

	% Find FGM data within this time interval
	ifg      = find(t_fg >= t0 & t_fg <= t1);
	t_fg     = t_fg(ifg);
	b_fg_bcs = b_fg_bcs(:, ifg);

	% Compute the averaged magnetic field
	b_avg_bcs = mms_edi_bavg(t_fg, b_fg_bcs, t_edi1, t_edi2, dt);

	% Clear uneeded data
	clear t0 t1 tvec tedge ifg

%------------------------------------%
% Beam Width                         %
%------------------------------------%
	% Rotation from BCS and EDI2 into EDI1 CS
	bcs2edi1  = mms_instr_xxyz2instr('BCS',  'EDI1');
	edi22edi1 = mms_instr_xxyz2instr('EDI2', 'EDI1');
	
	% Number of beam hits per gun
	n12 = length( t_edi1 );
	n21 = length( t_edi2 );
	
	% Rotate things into EDI1 CS
	b_avg_edi1   = mrvector_rotate( bcs2edi1,  b_avg_bcs.b_avg );
	b_gd12_edi1  = mrvector_rotate( bcs2edi1,  b_avg_bcs.b_gd12 );
	b_gd21_edi1  = mrvector_rotate( bcs2edi1,  b_avg_bcs.b_gd21 );
	fv_gd21_edi1 = mrvector_rotate( edi22edi1, fv_gd21_edi2 );
	pos_vg2_edi1 = mrvector_rotate( edi22edi1, pos_vg2_edi2 );
	
	% Combine data
	b_gdu_edi1 = [ b_gd12_edi1 b_gd21_edi1 ];
	fv_edi1    = [ fv_gd12_edi1 fv_gd21_edi1 ];
	pos_edi1   = [ repmat( pos_vg1_edi1, 1, n12 )  repmat( pos_vg2_edi1, 1, n21 ) ];
	recnum_gdu = [ b_avg_bcs.recnum_gd12 b_avg_bcs.recnum_gd21 ];
	gun_id     = [ ones(1, n12, 'uint8') ones(1, n21, 'uint8')+1 ];

	% Beam width
	beam_width = mms_edi_beam_width( fv_edi1, b_gdu_edi1, gun_id );

	% Clear unneeded data
	clear bcs2edi1 edi22edi1 b_gd12_edi1 b_gd21_edi1 fv_gd21_edi1 vg1_edi1 vg2_edi1

%------------------------------------%
% Grid of Targets                    %
%------------------------------------%
	% Setup the grid
	r    = [0.0 0.1   5.0];
	phi  = [0.0 0.1 360.0];
	grid = mms_edi_polar_grid(r, phi);

%------------------------------------%
% Step Thru Each Analysis Interval   %
%------------------------------------%
	% Allocate memory
	N       = size( b_avg_bcs.b_avg, 2 );
	d_bpp   = zeros(3, N);
	d_delta = zeros(1, N);
	
	for ii = 1 : N
		% Magnetic field and time for current interval
		t = b_avg_bcs.t_avg(ii);
		B = b_avg_edi1(:, ii);
		
		% Indices associated with interval ii.
		inds = find( recnum_gdu == b_avg_bcs.recnum(ii) );

	%------------------------------------%
	% Rotate Into BPP                    %
	%------------------------------------%
		% AVERAGE FIELD
		if method == 1
			xyz2bpp = mms_edi_xxyz2bpp( B );
		% INSTANTANEOUS FIELD
		elseif method == 2
			xyz2bpp = mms_edi_xxyz2bpp( b_gdu_edi1(inds) );
		else
			error( 'METHOD must be 1 or 2.' )
		end
			
		% Rotate from EDI1 to BPP
		fv_bpp  = mrvector_rotate( xyz2bpp, fv_edi1(:, inds )  );
		pos_bpp = mrvector_rotate( xyz2bpp, pos_edi1(:, inds ) );
		bw      = beam_width( inds );

	%------------------------------------%
	% Cost Function & Drift Step         %
	%------------------------------------%
		costFn = mms_edi_cost_function( fv_bpp, pos_bpp, bw, grid );
	
		% The virtual source point is located at the grid point
		% where the minimum of the cost function occurs. The drift
		% step points from the virtual source point toward the
		% origin.
		[delta, imin] =  min(costFn(:));
		d_bpp(1,ii)   = -grid.x(imin);
		d_bpp(2,ii)   = -grid.y(imin);
		d_delta(ii)   =  delta;
	end
	
	% Clear unneeded data
	clear N t B inds xyz2bpp fv_bpp pos_bpp bw costFn delta imin grid

%------------------------------------%
% Rotate Out of BPP                  %
%------------------------------------%
	% BPP -> BCS transformation
	bcs2bpp = mms_edi_xxyz2bpp( b_avg_bcs.b_avg );
	bpp2bcs = permute(bcs2bpp, [2 1 3]);
	
	% Transform data
	d_bcs = mrvector_rotate( bpp2bcs, d_bpp );

%------------------------------------%
% Drift Velocity & E-Field           %
%------------------------------------%
	% Caluclate the gyro-period
	%   - Replicate to 3xN so we can do element-wise math
	me     = 9.10938291e-31;                          % kg
	q      = 1.60217657e-19;                          % C
	B_mag  = mrvector_magnitude( b_avg_edi1 );
	T_gyro = repmat( (2*pi*me/q) ./ B_mag, 3, 1 );    % s

	% Drift Velocity
	%   - Convert m/s -> km/s
	v_drift_bcs = 1e-3 * (d_bcs ./ T_gyro);           % km/s
	
	% Electric Field
	%   - Convert from V/m to mV/m with 1e-3
	E_bcs = mrvector_cross( b_avg_bcs.b_avg, d_bcs ) ./ (1e3 * T_gyro);

%------------------------------------%
% Gather Data                        %
%------------------------------------%
	efield_cf = struct( 'tt2000_cf',  b_avg_bcs.t_avg,     ...
	                    'E_cf_bcs',   single(E_bcs),       ...
	                    'v_cf_bcs',   single(v_drift_bcs), ...
	                    'd_cf_bcs',   single(d_bcs),       ...
	                    'd_delta_cf', single(d_delta)      ...
	                  );
end