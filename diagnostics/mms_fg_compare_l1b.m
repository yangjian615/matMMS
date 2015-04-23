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

% Clear large variables before starting
clear b_l1b_bcs t_l1b b_bcs t_dn t_l1b_dn

%------------------------------------%
% Inputs                             %
%------------------------------------%
dfg_data_dir = '/Users/argall/Documents/Work/Data/MMS/DFG/';
afg_data_dir = '/Users/argall/Documents/Work/Data/MMS/DFG/';
fg_cal_dir   = '/Users/argall/Documents/Work/Data/MMS/FG_Cal/';
sc           = 'mms3';
instr        = 'dfg';
tstart       = '2015-04-18T00:00:00';
tend         = '2015-04-18T01:00:00';

if strcmp(instr, 'dfg')
	fg_data_dir = dfg_data_dir;
else
	fg_data_dir = afg_data_dir;
end

%------------------------------------%
% Calibrated Mag in BCS              %
%------------------------------------%
[b_bcs, t] = mms_fg_bcs(sc, instr, 'f128', tstart, tend, ...
                        'CalDir',  fg_cal_dir, ...
                        'DataDir', dfg_data_dir);

%------------------------------------%
% Get the Official l1b Data          %
%------------------------------------%
mode  = 'srvy';
level = 'l1b';

% Create the file name
l1b_fname = mms_construct_filename(sc, instr, mode, level, ...
                                   'Directory', fg_data_dir, ...
                                   'Tokens',    true);

% Search for the file
l1b_fname = MrFile_Search(l1b_fname, ...
                          'TStart',    tstart,   ...
                          'TEnd',      tend,     ...
                          'TimeOrder', '%Y%M%d', ...
                          'Closest',   true);

% Create variable names
b_vname = mms_construct_varname(sc, instr, mode, 'bcs');

% Get the data
[b_l1b_bcs, t_l1b] = MrCDF_Read(l1b_fname, b_vname, ...
                                'sTime', tstart, ...
                                'eTime', tend, ...
                                'ColumnMajor', true);

%------------------------------------%
% Plot the Results                   %
%------------------------------------%
% Convert time to datenumber to use the datetick function.
t_dn     = MrCDF_epoch2datenum(t);
t_l1b_dn = MrCDF_epoch2datenum(t_l1b);


f_bcs = figure();

% Magnitude
subplot(4,1,1)
plot( t_l1b_dn, b_l1b_bcs(4,:), t_dn, mrvector_magnitude(b_bcs) );
title([ 'Mag L1B vs. My L1B (' upper(instr) ', BCS)'] );
xlabel( 'Time UTC' );
ylabel( {'|B|', '(nT)'} );
datetick();
legend('MagTeam', 'Mine');

% X-component
subplot(4,1,2)
plot( t_l1b_dn, b_l1b_bcs(1,:), t_dn, b_bcs(1,:) );
title([ 'Mag L1B vs. My L1B (' upper(instr) ', BCS)'] );
xlabel( 'Time UTC' );
ylabel( {'B_{X}', '(nT)'} );
datetick();

% Y-component
subplot(4,1,3)
plot( t_l1b_dn, b_l1b_bcs(2,:), t_dn, b_bcs(2,:) );
title([ 'Mag L1B vs. My L1B (' upper(instr) ', BCS)'] );
xlabel( 'Time UTC' );
ylabel( {'B_{Z}', '(nT)'} );
datetick();

% Z-component
subplot(4,1,4)
plot( t_l1b_dn, b_l1b_bcs(3,:), t_dn, b_bcs(3,:) );
title([ 'Mag L1B vs. My L1B (' upper(instr) ', BCS)'] );
xlabel( 'Time UTC' );
ylabel( {'B_{Z}', '(nT)'} );
datetick();