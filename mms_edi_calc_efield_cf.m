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
function [E, v_ExB, d, d_delta] = ...
	mms_edi_calc_efield_cf( b_avg, pos_gdu, fv_gdu, beam_width, recnum, recnum_gdu, b_gdu )

	% Method of rotating into BPP
	%   - 1 = using averaged magnetic field
	%   - 2 = using instantaneous magnetic field
	method = 1;
	if nargin == 7
		method = 2;
	end

%------------------------------------%
% Grid of Targets                    %
%------------------------------------%
	% Setup the grid
	%   - Make the phi grid on "phi-step" less than 2*pi
	%   - Avoids duplicating points at phi = 0
	r    = [0.0 0.1   6.0];
	phi  = [0.0 0.1 359.9];
	grid = mms_edi_polar_grid(r, phi);

%------------------------------------%
% Step Thru Each Analysis Interval   %
%------------------------------------%
	% Allocate memory
	N       = size( b_avg, 2 );
	d_bpp   = zeros(3, N);
	d_delta = zeros(1, N);
	
	for ii = 1 : N
if mod( ii-1, 500 ) == 0
	fprintf('Beam %d of %d\n', ii, N)
end
		% Magnetic field for current interval
		B = b_avg(:, ii);
		
		% Indices associated with interval ii.
		inds = find( recnum_gdu == recnum(ii) );

	%------------------------------------%
	% Rotate Into BPP                    %
	%------------------------------------%
		% AVERAGE FIELD
		if method == 1
			xyz2bpp = mms_edi_xxyz2bpp( B );
		% INSTANTANEOUS FIELD
		elseif method == 2
			xyz2bpp = mms_edi_xxyz2bpp( b_gdu(inds) );
		else
			error( 'METHOD must be 1 or 2.' )
		end
			
		% Rotate from EDI1 to BPP
		fv_bpp  = mrvector_rotate( xyz2bpp, fv_gdu(:, inds )  );
		pos_bpp = mrvector_rotate( xyz2bpp, pos_gdu(:, inds ) );
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
	bcs2bpp = mms_edi_xxyz2bpp( b_avg );
	bpp2bcs = permute(bcs2bpp, [2 1 3]);
	
	% Transform data
	d = mrvector_rotate( bpp2bcs, d_bpp );

%------------------------------------%
% Drift Velocity & E-Field           %
%------------------------------------%
	% Caluclate the gyro-period
	%   - Replicate to 3xN so we can do element-wise math
	me     = 9.10938291e-31;                          % kg
	q      = 1.60217657e-19;                          % C
	B_mag  = mrvector_magnitude( b_avg );
	T_gyro = repmat( (2*pi*me/q) ./ B_mag, 3, 1 );    % s

	% Drift Velocity
	%   - Convert m/s -> km/s
	v_ExB = 1e-3 * (d ./ T_gyro);           % km/s
	
	% Electric Field
	%   - Convert from V/m to mV/m with 1e-3
	E = mrvector_cross( b_avg, d ) ./ (1e3 * T_gyro);
end