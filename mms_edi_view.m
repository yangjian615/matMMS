%
% Name
%   mms_edi_view
%
% Purpose
%   Take l1a EDI  data, despin, and rotate into GSE. Then, plot the 
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
	% Inputs
	edi_data_dir = '/Users/argall/Documents/Work/Data/MMS/EDI/';
	dfg_data_dir = '/Users/argall/Documents/Work/Data/MMS/DFG/';
	fg_cal_dir   = '/Users/argall/Documents/Work/Data/MMS/FG_Cal/';
	att_dir      = '/Users/argall/Documents/Work/Data/MMS/Attitude/';
	sunpulse_dir = '/Users/argall/Documents/Work/Data/MMS/HK/';
	sc           = 'mms3';
	instr        = 'edi';
	mode         = 'slow';
	level        = 'l1a_efield';
	tstart       = '2015-04-18T00:00:00';
	tend         = '2015-04-18T24:00:00';

%------------------------------------%
% AFG                                %
%------------------------------------%
	[gd12_dsl, gd21_dsl] = mms_edi_gse(sc, instr, mode, level, tstart, tend, ...
	                                   'SunPulseDir', sunpulse_dir, ...
	                                   'DataDir',     edi_data_dir );

	
	% Convert time to datenum
	t_gd12_datnum = MrCDF_epoch2datenum(gd12_dsl.t_gd12);
	t_gd21_datnum = MrCDF_epoch2datenum(gd21_dsl.t_gd21);

	% Convert to seconds for time range
	t_gd12_sse = MrCDF_epoch2sse(gd12_dsl.t_gd12);
	t_gd21_sse = MrCDF_epoch2sse(gd21_dsl.t_gd21, gd12_dsl.t_gd12(1));
end

it_gd12 = find( t_gd12_sse < 300 );
it_gd21 = find( t_gd21_sse < 300 );

%------------------------------------%
% Plot BCS                           %
%------------------------------------%

%
% BCS
%
f_bcs = figure();

% X-component
subplot(3,1,1)
plot( t_gd12_datnum(it_gd12), gd12_dsl.gun_gd12_dsl(1, it_gd12), ...
      t_gd21_datnum(it_gd21), gd21_dsl.gun_gd21_dsl(1, it_gd21) );
title([ 'Gun Positions in DMPA '] );
xlabel( 'Time UTC' );
ylabel( {'X', '(m)'} );
datetick();
legend('GD12', 'GD21');

% Y-component
subplot(3,1,2)
plot( t_gd12_datnum(it_gd12), gd12_dsl.gun_gd12_dsl(2, it_gd12), ...
      t_gd21_datnum(it_gd21), gd21_dsl.gun_gd21_dsl(2, it_gd21) );
title([ 'Gun Positions in DMPA '] );
xlabel( 'Time UTC' );
ylabel( {'Y', '(m)'} );
datetick();

% Z-component
subplot(3,1,3)
plot( t_gd12_datnum(it_gd12), gd12_dsl.gun_gd12_dsl(3, it_gd12), ...
      t_gd21_datnum(it_gd21), gd21_dsl.gun_gd21_dsl(3, it_gd21) );
title([ 'Gun Positions in DMPA '] );
xlabel( 'Time UTC' );
ylabel( {'Z', '(m)'} );
datetick();

