%
% Name
%   mms_fg_l1b_compare
%
% Purpose
%   Compare my level 1b data results to those of the mag team to
%   see if I am applying the calibration parameters correctly.
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-21      Written by Matthew Argall
%
%***************************************************************************

%------------------------------------%
% Inputs                             %
%------------------------------------%
dfg_data_dir = '/Users/argall/Documents/Work/Data/MMS/DFG/';
afg_data_dir = '/Users/argall/Documents/Work/Data/MMS/AFG/';
fg_cal_dir   = '/Users/argall/Documents/Work/Data/MMS/FG_Cal/';
sc           = 'mms3';
instr        = 'dfg';
tstart       = '2015-04-18T05:00:00';
tend         = '2015-04-18T05:03:00';

if strcmp(instr, 'dfg')
	fg_data_dir = dfg_data_dir;
	mode        = 'f128';
else
	fg_data_dir = afg_data_dir;
	mode        = 'fast';
end

%------------------------------------%
% Find Files                         %
%------------------------------------%
	% FG L1A Data File
	[l1a_fname, count, str] = mms_file_search(sc, instr, mode, 'l1a', ...
	                                          'TStart',    tstart, ...
	                                          'TEnd',      tend, ...
	                                          'Directory', fg_data_dir);
	assert(count > 0, ['FG file not found: "' str '".']);
	
	% FG L1B Data File
	[l1b_fname, count, str] = mms_file_search(sc, instr, 'srvy', 'l1b', ...
	                                          'TStart',    tstart, ...
	                                          'TEnd',      tend, ...
	                                          'Directory', fg_data_dir);
	assert(count > 0, ['FG file not found: "' str '".']);
	
	% FG Hi-Cal File
	[hical_fname, count, str] = mms_file_search(sc, instr, 'hirangecal', 'l2pre', ...
	                                            'TStart',    tstart, ...
	                                            'TEnd',      tend, ...
	                                            'Directory', fg_cal_dir);
	assert(count > 0, ['FG Hi-Cal file not found: "' str '".']);
	
	% FG Lo-Cal File
	[local_fname, count, str] = mms_file_search(sc, instr, 'lorangecal', 'l2pre', ...
	                                            'TStart',    tstart, ...
	                                            'TEnd',      tend, ...
	                                            'Directory', fg_cal_dir);
	assert(count > 0, ['FG Hi-Cal file not found: "' str '".']);

%------------------------------------%
% Calibrated Mag in BCS              %
%------------------------------------%
[t, b_bcs] = mms_fg_bcs(l1a_fname, hical_fname, local_fname, tstart, tend);

%------------------------------------%
% Get the Official l1b Data          %
%------------------------------------%

% Get the data
fg_l1b    = mms_fg_read_l1b(l1b_fname, tstart, tend);
t_l1b     = fg_l1b.tt2000;
b_l1b_bcs = fg_l1b.b_bcs;

% Clear structure
clear fg_l1b

%------------------------------------%
% Plot the Results                   %
%------------------------------------%
% Convert time to datenumber to use the datetick function.
t_dn     = MrCDF_epoch2datenum(t');
t_l1b_dn = MrCDF_epoch2datenum(t_l1b');


f_bcs = figure();

% Magnitude
subplot(4,1,1)
plot( t_l1b_dn, b_l1b_bcs(4,:), t_dn, mrvector_magnitude(b_bcs) );
title([ upper(sc) ' ' upper(instr) ' BCS ' tstart(1:10)] );
xlabel( 'Time UTC' );
ylabel( {'|B|', '(nT)'} );
datetick();
legend('MagTeam', 'Mine');

% X-component
subplot(4,1,2)
plot( t_l1b_dn, b_l1b_bcs(1,:), t_dn, b_bcs(1,:) );
xlabel( 'Time UTC' );
ylabel( {'B_{X}', '(nT)'} );
datetick();

% Y-component
subplot(4,1,3)
plot( t_l1b_dn, b_l1b_bcs(2,:), t_dn, b_bcs(2,:) );
xlabel( 'Time UTC' );
ylabel( {'B_{Y}', '(nT)'} );
datetick();

% Z-component
subplot(4,1,4)
plot( t_l1b_dn, b_l1b_bcs(3,:), t_dn, b_bcs(3,:) );
xlabel( 'Time UTC' );
ylabel( {'B_{Z}', '(nT)'} );
datetick();