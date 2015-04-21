%
% Name
%   mms_mag_view
%
% Purpose
%   Take l1a searchcoil magnetometer data, apply calibration
%   parameters, despin, and rotate into GSE. Then, plot the 
%   results of each stage.
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-14      Written by Matthew Argall
%

get_data = true;

if get_data
%------------------------------------%
% Inputs                             %
%------------------------------------%
	afg_data_dir = '/Users/argall/Documents/Work/Data/MMS/AFG/';
	dfg_data_dir = '/Users/argall/Documents/Work/Data/MMS/DFG/';
	fg_cal_dir   = '/Users/argall/Documents/Work/Data/MMS/FG_Cal/';
	sc_data_dir  = '/Users/argall/Documents/Work/Data/MMS/SCM/';
	sc_cal_dir   = '/Users/argall/Documents/Work/Data/MMS/SCM_Cal/';
	att_dir      = '/Users/argall/Documents/Work/Data/MMS/Attitude/';
	sc           = 'mms2';
	tstart       = '2015-03-20T15:00:00';
	tend         = '2015-03-20T16:00:00';
	duration     = 600.0;

%------------------------------------%
% AFG                                %
%------------------------------------%
	% Instrument and mode
	instr = 'afg';
	mode  = 'fast';

	% Get data
	[t_afg, b_afg_gse, b_afg_dsl, b_afg_smpa, b_afg_bcs, b_afg_omb] ...
		= mms_fg_gse(sc, instr, mode, tstart, tend, ...
		             'AttDir',   att_dir, ...
		             'DataDir',  afg_data_dir, ...
		             'CalDir',   fg_cal_dir);

%------------------------------------%
% DFG                                %
%------------------------------------%
	% Instrument and mode
	instr = 'dfg';
	mode  = 'f128';

	% Get data
	[t_dfg, b_dfg_gse, b_dfg_dsl, b_dfg_smpa, b_dfg_bcs, b_dfg_omb] ...
		= mms_fg_gse(sc, instr, mode, tstart, tend, ...
		             'AttDir',   att_dir, ...
		             'DataDir',  dfg_data_dir, ...
		             'CalDir',   fg_cal_dir);

%------------------------------------%
% SCM                                %
%------------------------------------%
	% Instrument and mode
	instr = 'scm';
	mode  = 'comm';

	% Get the data
	[t_scm, b_scm_gse, b_scm_dsl, b_scm_smpa, b_scm_bcs] ...
		= mms_sc_gse(sc, instr, mode, tstart, tend, ...
		             'AttDir',   att_dir, ...
		             'DataDir',  sc_data_dir, ...
		             'CalDir',   sc_cal_dir, ...
		             'Duration', duration);

%------------------------------------%
% Plot                               %
%------------------------------------%
	
	% Convert time to datenum
	t_dfg_datnum = MrCDF_epoch2datenum(t_dfg);
	t_afg_datnum = MrCDF_epoch2datenum(t_afg);
	t_scm_datnum = MrCDF_epoch2datenum(t_scm);

	% Convert to seconds for time range
	t_dfg_sse = MrCDF_epoch2sse(t_dfg);
	t_afg_sse = MrCDF_epoch2sse(t_afg, t_dfg(1));
	t_scm_sse = MrCDF_epoch2sse(t_scm, t_dfg(1));
end

it_dfg = find( t_dfg_sse < 300 );
it_afg = find( t_afg_sse < 300 );
it_scm = find( t_scm_sse < 300 );

%------------------------------------%
% Plot BCS                           %
%------------------------------------%

%
% BCS
%
f_bcs = figure();

% X-component
subplot(3,1,1)
plot( t_dfg_datnum(it_dfg), b_dfg_bcs(1,it_dfg), ...
      t_afg_datnum(it_afg), b_afg_bcs(1,it_afg), ...
      t_scm_datnum(it_scm), b_scm_bcs(1,it_scm) );
title([ 'Magnetometers in BCS '] );
xlabel( 'Time UTC' );
ylabel( {'B_{X}', '(nT)'} );
datetick();
legend('B_{DFG}', 'B_{AFG}', 'B_{SCM}');

% Y-component
subplot(3,1,2)
plot( t_dfg_datnum(it_dfg), b_dfg_bcs(2,it_dfg), ...
      t_afg_datnum(it_afg), b_afg_bcs(2,it_afg), ...
      t_scm_datnum(it_scm), b_scm_bcs(2,it_scm));
title([ 'Magnetometers in BCS '] );
xlabel( 'Time UTC' );
ylabel( {'B_{Y}', '(nT)'} );
datetick();

% Z-component
subplot(3,1,3)
plot( t_dfg_datnum(it_dfg), b_dfg_bcs(3,it_dfg), ...
      t_afg_datnum(it_afg), b_afg_bcs(3,it_afg), ...
      t_scm_datnum(it_scm), b_scm_bcs(3,it_scm) );
title([ 'Magnetometers in BCS '] );
xlabel( 'Time UTC' );
ylabel( {'B_{Z}', '(nT)'} );
datetick();

%------------------------------------%
% Plot DSL                           %
%------------------------------------%

%
% DSL
%
f_dsl = figure();

% X-component
subplot(3,1,1)
plot( t_dfg_datnum(it_dfg), b_dfg_dsl(1,it_dfg), ...
      t_afg_datnum(it_afg), b_afg_dsl(1,it_afg), ...
      t_scm_datnum(it_scm), b_scm_dsl(1,it_scm) );
title([ 'Magnetometers in DSL '] );
xlabel( 'Time UTC' );
ylabel( {'B_{X}', '(nT)'} );
datetick();
legend('B_{DFG}', 'B_{AFG}', 'B_{SCM}');

% Y-component
subplot(3,1,2)
plot( t_dfg_datnum(it_dfg), b_dfg_dsl(2,it_dfg), ...
      t_afg_datnum(it_afg), b_afg_dsl(2,it_afg), ...
      t_scm_datnum(it_scm), b_scm_dsl(2,it_scm) );
title([ 'Magnetometers in DSL '] );
xlabel( 'Time UTC' );
ylabel( {'B_{Y}', '(nT)'} );
datetick();

% Z-component
subplot(3,1,3)
plot( t_dfg_datnum(it_dfg), b_dfg_dsl(3,it_dfg), ...
      t_afg_datnum(it_afg), b_afg_dsl(3,it_afg), ...
      t_scm_datnum(it_scm), b_scm_dsl(3,it_scm) );
title([ 'Magnetometers in DSL '] );
xlabel( 'Time UTC' );
ylabel( {'B_{Z}', '(nT)'} );
datetick();

%------------------------------------%
% Plot GSE                           %
%------------------------------------%

%
% GSE
%
f_gse = figure();

% X-component
subplot(3,1,1)
plot( t_dfg_datnum(it_dfg), b_dfg_gse(1,it_dfg), ...
      t_afg_datnum(it_afg), b_afg_gse(1,it_afg), ...
      t_scm_datnum(it_scm), b_scm_gse(1,it_scm) );
title([ 'Magnetometers in GSE '] );
xlabel( 'Time UTC' );
ylabel( {'B_{X}', '(nT)'} );
datetick();
legend('B_{DFG}', 'B_{AFG}', 'B_{SCM}');

% Y-component
subplot(3,1,2)
plot( t_dfg_datnum(it_dfg), b_dfg_gse(2,it_dfg), ...
      t_afg_datnum(it_afg), b_afg_gse(2,it_afg), ...
      t_scm_datnum(it_scm), b_scm_gse(2,it_scm) );
title([ 'Magnetometers in GSE '] );
xlabel( 'Time UTC' );
ylabel( {'B_{Y}', '(nT)'} );
datetick();

% Z-component
subplot(3,1,3)
plot( t_dfg_datnum(it_dfg), b_dfg_gse(3,it_dfg), ...
      t_afg_datnum(it_afg), b_afg_gse(3,it_afg), ...
      t_scm_datnum(it_scm), b_scm_gse(3,it_scm) );
title([ 'Magnetometers in GSE '] );
xlabel( 'Time UTC' );
ylabel( {'B_{Z}', '(nT)'} );
datetick();

