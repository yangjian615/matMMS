% edi_drift_step using beam convergence, v1.02 (show also on plots)
%
% v0102: apply plot style during save
%        include obsID on plots
%        include local copies of myLibAppConstants, myLibScienceConstants
%        changed references to Cluster data
% v0101: update plot: sample datetime, selected confidence, legend
%        implement weighted mean: sin^2
%        tested with 'mms2_edi_slow_ql_efield_20150509_v0.0.1.cdf'
% v0100: draft, proof of concept

clc             % clear the command window
clf
clear variables
close all       % close all figures
format compact
format short g  % +, bank, hex, long, rat, short, short g, short eng
myLibAppConstants % custom colors; set default axis colors

global dotVersion; dotVersion = 'v1.02';

% Get all the data from (ex.) 
% 'mms2_edi_slow_ql_efield_20150509_v0.0.1.cdf'
% 'mms2_edp_comm_ql_dce2d_20150509000000_v0.0.0.cdf'
edi_drift_step_read_ql__EDI__B__EDP_data

nB_recs        = length (B_datenum);
driftStep      = zeros (3, nB_recs);
driftStepSigma = zeros (3, nB_recs);
driftVelocity  = zeros (3, nB_recs);
drift_E        = zeros (3, nB_recs);

% for each B interval, there are EDI records indexed to B ~ many EDI:1 B
% Here we look at each B record, find the corresponding EDI records, and filter.
disp 'looping over edi_B_dmpa'
for igdxx = 1: length (edi_B_dmpa)
	B_dmpa   = edi_B_dmpa (1:3,igdxx);
	if ~isnan (B_dmpa (1))

		% B field data
		B_tt2000 = edi_BdvE_tt2000 (igdxx);

		% GDU data that corresponds to the B field data: position and firing vectors
		% Includes both GDUs
		iigd12_b_avgIntrp = find (edi_gd_B_index == igdxx);

		gd_virtual_dmpa   = edi_gd_virtual_dmpa (:, iigd12_b_avgIntrp);
		gd_fv_dmpa        = edi_gd_fv_dmpa      (:, iigd12_b_avgIntrp);
		gd_ID             = edi_gd_ID           (:, iigd12_b_avgIntrp);

		% keyboard
		if (length (iigd12_b_avgIntrp) > 2) % More than 2 beams
			[ driftStep(:, igdxx), ...
			  driftStepSigma(:, igdxx), ...
			  driftVelocity(:, igdxx), ...
			  drift_E(:, igdxx) ] = edi_drift_step ( ...
				B_tt2000, ...
				B_dmpa, ...
				gd_virtual_dmpa, ...
				gd_fv_dmpa, ...
				obsID, ...
				gd_ID );
		else
			drift_E (:, igdxx) = [ NaN; NaN; NaN ];
		end

	end
end

[ hAxis hEDP_dce_xyz_dsl hEDI_B ] = plotyy ( ...
	edp_datenum, edp_dce_xyz_dsl(:, 1:2), ...
	B_datenum, edi_B_dmpa (1:3, :)', @plot, @plot );
% [ hAxis ] = plot ( ...
% 	edp_datenum, edp_dce_xyz_dsl(:, 1:2));

plotDateMin = double (min (min(edp_datenum), min(B_datenum)));
plotDateMax = double (max (max(edp_datenum), max(B_datenum)));

datenumOneMin = 0.0006944444496185;
datenumOneHr  = 60.0 * datenumOneMin;
% xlim ( [ B_datenum(1) B_datenum(end) ] )
% set (gca, 'XTick', [B_datenum(1): datenumOneHr: B_datenum(end)])
xlim ( [ plotDateMin plotDateMax ] )
set (gca, 'XTick', [ plotDateMin: datenumOneHr: plotDateMax]) % debug?

datetick ('x', 'HH:MM', 'keeplimits', 'keepticks') % debug?

set (hAxis(2), 'XTick', [])

hold on
% set (hEDI_B, 'Color', [ MMS_plotColorx, MMS_plotColory, MMS_plotColorz ]);
set (hEDI_B(1), 'Color', MMS_plotColorx);
set (hEDI_B(2), 'Color', MMS_plotColory);
set (hEDI_B(3), 'Color', MMS_plotColorz);

hDrift_Ex = plot (B_datenum, drift_E (1,:), ...
	'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', MMS_plotColorx, 'MarkerSize', 5.0);
% datetick ()

hDrift_Ey = plot (B_datenum, drift_E (2,:), ...
	'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', myDarkGreen, 'MarkerEdgeColor', MMS_plotColory, 'MarkerSize', 5.0);
% datetick ()
title 'SDP 2D E-fields, EDI E-fields and B avg @ 5 s intervals, DMPA'
xlabel (sprintf ( '%s', datestr( spdftt2000todatenum (edp_tt2000(1)), 'yyyy-mm-dd') ))
ylabel (hAxis (1),'mV-m^-^1')
ylabel (hAxis (2),'nT')
% legend ('E_x SDP', 'E_y SDP', 'E_x EDI', 'E_y EDI', 'B_x', 'B_y', 'B_z' );
legend ('SDP E_x', 'SDP E_y', 'EDI E_x', 'EDI E_y', 'B_x', 'B_y', 'B_z' );

hold off


% b_avg_tt2000 = b_avg.t_avg; % MATLAB does not save structure.variables
% save 'mms4_edi_slow_l1a_efield_20150506_SDP_and_EDI_driftstep_E_field_2D.mat' edp_tt2000 edp_dce_xyz_dsl b_avg_tt2000 drift_E
% load 'mms4_edi_slow_l1a_efield_20150506_SDP_and_EDI_driftstep_E_field_2D.mat'

%{

Some of the variables in the mat file
% Restore the data
mat_file = 'mms4_edi_slow_l1a_efield_20150506_v0.1.0.mat';
load (mat_file)

% Average Quantities
tt2000_avg = b_avg.t_avg;        % time tags of 5-second avg B-field
dt_avg     = b_avg.dt_avg;       % delta plus of t_avg
b_average  = b_avg.b_avg;        % Avg B-field
b_std      = b_avg.b_std;        % Standard deviation of b_avg
b_gd12     = b_avg.b_gd12;       % B interpolated to Gun1 times
b_gd21     = b_avg.b_gd21;       % B interpolated to Gun2 times
inds_gd12  = b_avg.inds_gd12;    % Array of indices indicating to which b_avg value a firing vector maps.
inds_gd21  = b_avg.inds_gd21;    % Array of indices indicating to which b_avg value a firing vector maps.

% EDI data (the structure has more)
tt2000_gd12       = edi.epoch_gd12;              % Epoch times of gun1
tt2000_gd21       = edi.epoch_gd21;              % Epoch times of gun2
virtual_gun1_dmpa = edi.virtual_gun1_dmpa;       % Location of gun1 on virtual spacecraft
virtual_gun2_dmpa = edi.virtual_gun2_dmpa;       % Location of gun2 on virtual spacecraft
fv_gd12_dmpa      = edi.fv_gd12_dmpa;            % Firing vectors from gun1
fv_gd21_dmpa      = edi.fv_gd21_dmpa;            % Firing vectors from gun2

% FG Data
fg_tt2000 = fg_ql.tt2000;       % Epoch times
b_dmpa    = fg_ql.b_dmpa;       % Magnetic field

% EDP Data
edp_tt2000 = edp_ql.tt2000;     % Epoch times
edp_dce_xyz_dsl      = edp_ql.edp_dce_xyz_dsl;      % Electric field in DSL coordinates (DSL ~= DMPA)

%}