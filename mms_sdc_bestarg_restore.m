%
% Name
%   mms_sdc_bestarg_restore
%
% Purpose
%   Restore data from a mat file created by mms_sdc_bestarg.m and translate
%   structure elements to variables.
%
% Calling Sequence
%   [] = mms_sdc_bestarg_restore()
%     Restore data from a mat file created by mms_sdc_bestarg.m.
%
% Parameters
%
% Returns
%   SAVE_FILE       out, required, type=string
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-06-03      Written by Matthew Argall
%
function [] = mms_sdc_bestarg_restore()

	% Restore the data
	mat_file = '/home/argall/data/mms/matfiles/mms4_edi_slow_l1a_efield_20150506_v0.1.0.mat';
	load(mat_file);
	
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
end