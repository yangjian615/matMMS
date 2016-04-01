%
% Name
%   mms_edi_calc_efield_bc
%
% Purpose
%   Compute the electric field from EDI beam data using the
%   beam convergence method.
%
% Calling Sequence
%   EFIELD_BC = mms_edi_calc_efield_bc(EDI, B_AVG)
%     Compute the EDI electric field EFIELD_BC using the beam
%     convertence (BC) method. Input data, EDI, should be a
%     structure of the type returned by mms_edi_create_l2pre_efield.m
%     and should be inerpolated onto intervals of magnetic
%     field data B_AVG, which is a structure of the type returned
%     by mms_edi_bavg.m.
%
% Parameters
%   EDI:            in, required, type=structure
%   B_AVG:          in, required, type=structure
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
function [E v_ExB d d_std q] = ...
	mms_edi_calc_efield_bc( t_avg, b_avg, pos_gdu, fv_gdu, recnum, recnum_gdu, gun_id )

%------------------------------------%
% Beam Convergence                   %
%------------------------------------%
	szOut = size(b_avg);

	% Allocate memory
	E     = zeros( szOut, 'double');
	v_ExB = zeros( szOut, 'double');
	d     = zeros( szOut, 'double');
	d_std = zeros( szOut, 'double');
	q     = zeros( 1, szOut(2), 'uint8');

	% Declare global variables for the drift step function
	global plot_beams;      plot_beams            = false;
	global rotation_method; rotation_method       = 2;
	global use_v10502;      use_v10502            = false;
	                        sinx_wt_Q_xovr_angles = [ 8.0 30 ];
	global sinx_wt_Q_xovr;  sinx_wt_Q_xovr        = sind (sinx_wt_Q_xovr_angles).^4.0; % breakpoints for quality ranges for sin^x weighting

	% Step through each interval
	for ii = 1 : szOut(2)
		% Find data for current interval
		inds = find( recnum_gdu == recnum(ii) );

		% Compute the electric field
		if length(inds) > 1
			[ d_temp d_std_temp v_temp E_temp q_temp ] ...
				= edi_drift_step( '',                  ...
				                  t_avg(ii),           ...
				                  b_avg(:,ii),         ...
				                  pos_gdu(:,inds),     ...
				                  fv_gdu(:,inds),      ...
				                  gun_id(inds) );

			% Replace NaN with fill value
			if isnan(E_temp(1))
				E_temp     = [ -1e31; -1e31; -1e31 ];
				v_temp     = [ -1e31; -1e31; -1e31 ];
				d_temp     = [ -1e31; -1e31; -1e31 ];
				d_std_temp = [ -1e31; -1e31; -1e31 ];
				q_temp     = 0;
			end
		else
			E_temp     = [ -1e31; -1e31; -1e31 ];
			v_temp     = [ -1e31; -1e31; -1e31 ];
			d_temp     = [ -1e31; -1e31; -1e31 ];
			d_std_temp = [ -1e31; -1e31; -1e31 ];
			q_temp     = 0;
		end

		% Store results
		E(:, ii)     = E_temp;
		v_ExB(:, ii) = v_temp;
		d(:, ii)     = d_temp;
		d_std(:, ii) = d_std_temp;
		q(ii)        = q_temp;
	end
end