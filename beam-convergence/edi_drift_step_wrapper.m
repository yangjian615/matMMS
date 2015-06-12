clc             % clear the command window
clf
clear variables
close all       % close all figures
format compact
format short g  % +, bank, hex, long, rat, short, short g, short eng
myColors        % custom colors; set default axis colors

% Restore data from a mat file created by mms_sdc_bestarg.m
load ('mms4_edi_slow_l1a_efield_20150506_v0.1.0.mat')
E_dmpa = zeros (3, length (b_avg.b_avg), 'double');

n = 4;
% 	for igdxx = n: n % length (b_average)
for igdxx = 1: length (b_avg.b_avg)
	% B field data
	b_tt2000   = b_avg.t_avg (    igdxx);
	b_avg_dmpa = b_avg.b_avg (1:3,igdxx);

	% GDU data that corresponds to the B field data: position and firing vectors
	iigd12_b_avgIntrp     = find (b_avg.inds_gd12 == igdxx);
	edi_gun1_quality      = edi.quality_gd12 (iigd12_b_avgIntrp);
	edi_gun1_virtual_dmpa = edi.virtual_gun1_dmpa (:, iigd12_b_avgIntrp);
	edi_gd12_fv_dmpa      = edi.fv_gd12_dmpa      (:, iigd12_b_avgIntrp);
	edi_gd12_tof          = edi.tof_gd12 (iigd12_b_avgIntrp);

	iigd21_b_avgIntrp     = find (b_avg.inds_gd21 == igdxx);
	edi_gun2_quality      = edi.q_gd21 (iigd21_b_avgIntrp);
	edi_gun2_virtual_dmpa = edi.virtual_gun2_dmpa (:, iigd21_b_avgIntrp);
	edi_gd21_fv_dmpa      = edi.fv_gd21_dmpa      (:, iigd21_b_avgIntrp);
	edi_gd21_tof          = edi.tof_gd21 (iigd21_b_avgIntrp);

	iq1 = find (edi_gun1_quality < 3);
	iq2 = find (edi_gun2_quality < 3);

% keyboard

	edi_gun1_virtual_dmpa (:, iq1) = [];
	edi_gd12_fv_dmpa      (:, iq1) = [];
	edi_gd12_tof          (:, iq1) = [];
	edi_gun2_virtual_dmpa (:, iq2) = [];
	edi_gd21_fv_dmpa      (:, iq2) = [];
	edi_gd21_tof          (:, iq2) = [];

% keyboard
	if (length (edi_gun1_virtual_dmpa) > 2) & (length (edi_gun2_virtual_dmpa) > 2)
		E_dmpa (:, igdxx) = edi_drift_step ( ...
			b_tt2000, ...
			b_avg_dmpa, ...
			edi_gun1_virtual_dmpa, ...
			edi_gd12_fv_dmpa, ...
			edi_gun2_virtual_dmpa, ...
			edi_gd21_fv_dmpa, ...
			edi_gd12_tof, ...
			edi_gd21_tof );
	else
		E_dmpa (:, igdxx) = [ NaN; NaN; NaN ];
	end
end

edp_tt2000 = edp_ql.tt2000; % Epoch times
e_dsl      = edp_ql.e_dsl;  % Electric field in DSL coordinates (DSL ~= DMPA)

edp_epoch  = double (edp_tt2000) * 1e-9; % ~> sec
bAvg_epoch = double (b_avg.t_avg) * 1e-9;
[ datestr(spdftt2000todatenum(edp_tt2000(1)), 'yyyy-mm-dd HH:MM:ss'), ' ',...
  datestr(spdftt2000todatenum(edp_tt2000(end)), 'yyyy-mm-dd HH:MM:ss') ]
edp_datenum  = spdftt2000todatenum(edp_tt2000');
bAvg_datenum = spdftt2000todatenum(b_avg.t_avg');

[ hAxis hE_dsl hB_avg ] = plotyy ( ...
	edp_datenum, e_dsl(1:2,:)', ...
	bAvg_datenum, b_avg.b_avg (1:3,:)' );

datenumOneMin = 0.0006944444496185;
xlim ( [ datenum('2015-05-06 15:30:00') datenum('2015-05-06 15:35:00') ] )
 set (gca, 'XTick', [datenum('2015-05-06 15:30:00'): datenumOneMin: datenum('2015-05-06 15:35:00')])

datetick ('x', 'HH:MM', 'keeplimits', 'keepticks')

set (hAxis(2), 'XTick', [])

hold on
% set (hB_avg, 'Color', [ MMS_plotColorx, MMS_plotColory, MMS_plotColorz ]);
set (hB_avg(1), 'Color', MMS_plotColorx);
set (hB_avg(2), 'Color', MMS_plotColory);
set (hB_avg(3), 'Color', MMS_plotColorz);

plot (bAvg_datenum, E_dmpa (1,:), ...
	'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', MMS_plotColorx, 'MarkerEdgeColor', MMS_plotColorx, 'MarkerSize', 5.0);
% datetick ()

plot (bAvg_datenum, E_dmpa (2,:), ...
	'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', MMS_plotColory, 'MarkerEdgeColor', MMS_plotColory, 'MarkerSize', 5.0);
% datetick ()
title 'SDP 2D E-fields, EDI E-fields and B avg @ 5 s intervals, DMPA'
xlabel (sprintf ( '%s', datestr( spdftt2000todatenum (edp_tt2000(1)), 'yyyy-mm-dd') ))
ylabel (hAxis (1),'mV-m^-^1')
ylabel (hAxis (2),'nT')
legend ('E_x SDP', 'E_y SDP', 'E_x EDI', 'E_y EDI', 'B_x', 'B_y', 'B_z' );

hold off

% b_avg_tt2000 = b_avg.t_avg; % MATLAB does not save structure.variables
% save 'mms4_edi_slow_l1a_efield_20150506_SDP_and_EDI_driftstep_E_field_2D.mat' edp_tt2000 e_dsl b_avg_tt2000 E_dmpa
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
e_dsl      = edp_ql.e_dsl;      % Electric field in DSL coordinates (DSL ~= DMPA)

%}