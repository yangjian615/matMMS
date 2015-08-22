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
function [efield_bc, b_avg_dmpa] = mms_edi_calc_efield_bc( t_fg, b_fg_dmpa, ...
                                                           t_gd12, pos_vg1_dmpa, fv_gd12_dmpa, ...
                                                           t_gd21, pos_vg2_dmpa, fv_gd21_dmpa, dt )

	if nargin < 3
		dt = 5.0;
	end

%------------------------------------%
% Compuate Average B                 %
%------------------------------------%
	% Time range that we have data
	if isempty(t_gd12) && isempty(t_gd21)
		error( 'No EDI E-field data available.' );
	elseif isempty(t_gd12)
		t0 = t_gd21(1);
		t1 = t_gd21(end);
	elseif isempty(t_gd21)
		t0 = t_gd12(1);
		t1 = t_gd12(end);
	else
		t0 = min( [t_gd12(1)   t_gd21(1)  ] );
		t1 = max( [t_gd12(end) t_gd21(end)] );
	end
	
	% Breakdown into time vectors
	tvec = MrCDF_Epoch_Breakdown( [t0, t1] );
	
	% Round down to the nearest DT seconds and recompute
	tvec(:,6)     = tvec(:,6) - mod(tvec(:,6), dt);
	tvec(:,7:end) = 0.0;
	tedge         = MrCDF_Epoch_Compute(tvec);
	t0            = tedge(1);
	t1            = tedge(2) + int64(dt * 1d9);

	% Find FGM data within this time interval
	ifg       = find(t_fg >= t0 & t_fg <= t1);
	t_fg      = t_fg(ifg);
	b_fg_dmpa = b_fg_dmpa(:, ifg);

	% Compute the averaged magnetic field
	b_avg_dmpa = mms_edi_bavg(t_fg, b_fg_dmpa, t_gd12, t_gd21, dt);

%------------------------------------%
% Beam Convergence                   %
%------------------------------------%
	szOut = size(b_avg_dmpa.b_avg);

	% Allocate memory
	E_dmpa     = zeros( szOut, 'double');
	v_dmpa     = zeros( szOut, 'double');
	d_dmpa     = zeros( szOut, 'double');
	d_std_dmpa = zeros( szOut, 'double');
	quality    = zeros( 1, szOut(2), 'uint8');

	% Declare global variables for the drift step function
	global plot_beams;      plot_beams      = false;
	global rotation_method; rotation_method = 2;
	global use_v10502;      use_v10502      = false;

	% Step through each interval
	for ii = 1 : length(b_avg_dmpa.recnum)
		recnum = b_avg_dmpa.recnum(ii);
		
		% B field data
		B_tt2000 = b_avg_dmpa.t_avg(ii);
		B_dmpa   = b_avg_dmpa.b_avg(:, ii);

		% GDU data that corresponds to the B field data: position and firing vectors
		iigd12_b_avgIntrp = find( b_avg_dmpa.recnum_gd12 == recnum );
		iigd21_b_avgIntrp = find( b_avg_dmpa.recnum_gd21 == recnum );

		% Virtual gun positions
		gd_virtual_dmpa = [ pos_vg1_dmpa(:, iigd12_b_avgIntrp), ...
		                    pos_vg2_dmpa(:, iigd21_b_avgIntrp) ];
		
		% Firing Vectors
		gd_fv_dmpa = [ fv_gd12_dmpa(:, iigd12_b_avgIntrp), ...
		               fv_gd21_dmpa(:, iigd21_b_avgIntrp) ];
		
		% Gun ID of each data point
		n_gd12                = size( iigd12_b_avgIntrp, 2 );
		gd_ID                 = ones( 1, size(gd_virtual_dmpa, 2) );
		gd_ID( n_gd12+1:end ) = 2;

		% Compute the electric field
		if size(gd_virtual_dmpa, 2) > 2
			[ d_temp d_std_temp v_temp E_temp q_temp] ...
				= edi_drift_step( '', ...
				                  B_tt2000, ...
				                  B_dmpa, ...
				                  gd_virtual_dmpa, ...
				                  gd_fv_dmpa, ...
				                  gd_ID );

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
		E_dmpa(:, ii)     = E_temp;
		v_dmpa(:, ii)     = v_temp;
		d_dmpa(:, ii)     = d_temp;
		d_std_dmpa(:, ii) = d_std_temp;
		quality(ii)       = q_temp;
	end

%------------------------------------%
% Gather Data                        %
%------------------------------------%
	efield_bc = struct( 'tt2000_bc',     b_avg_dmpa.t_avg,    ...
	                    'E_bc_dmpa',     single(E_dmpa),      ...
	                    'v_bc_dmpa',     single(v_dmpa),      ...
	                    'd_bc_dmpa',     single(d_dmpa),      ...
	                    'd_std_bc_dmpa', single(d_std_dmpa),  ...
	                    'quality_bc',    quality              ...
	                  );
end