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
clear t b_dmpa b_ql_dmpa t_ql t_dn t_ql_dn

%------------------------------------%
% Inputs                             %
%------------------------------------%
dfg_data_dir = '/Users/argall/Documents/Work/Data/MMS/DFG/';
afg_data_dir = '/Users/argall/Documents/Work/Data/MMS/AFG/';
fg_cal_dir   = '/Users/argall/Documents/Work/Data/MMS/FG_Cal/';
att_dir      = ''; % '/Users/argall/Documents/Work/Data/MMS/Attitude/';
sunpulse_dir = '/Users/argall/Documents/Work/Data/MMS/HK/';
sc           = 'mms3';
instr        = 'afg';
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
% Calibrated Mag in BCS              %
%------------------------------------%
[t, ~, b_dmpa] = mms_fg_gse(sc, instr, mode, tstart, tend, ...
                            'CalDir',      fg_cal_dir,     ...
                            'DataDir',     fg_data_dir,    ...
                            'AttDir',      att_dir ,       ...
                            'SunPulseDir', sunpulse_dir);

%------------------------------------%
% Get the Official QL Data           %
%------------------------------------%
mode  = 'srvy';
level = 'ql';

% Create the file name
ql_fname = mms_construct_filename(sc, instr, mode, level,   ...
                                  'Directory', fg_data_dir, ...
                                  'Tokens',    true);

% Search for the file
ql_fname = MrFile_Search(ql_fname, ...
                         'TStart',    tstart,   ...
                         'TEnd',      tend,     ...
                         'TimeOrder', '%Y%M%d', ...
                         'Closest',   true);

% Create variable names
b_vname = mms_construct_varname(sc, instr, mode, 'dmpa');

% Get the data
[b_ql_dmpa, t_ql] = MrCDF_Read(ql_fname, b_vname, ...
                               'sTime', tstart, ...
                               'eTime', tend, ...
                               'ColumnMajor', true);

%------------------------------------%
% Plot the Results                   %
%------------------------------------%
% Convert time to datenumber to use the datetick function.
t_dn    = MrCDF_epoch2datenum(t);
t_ql_dn = MrCDF_epoch2datenum(t_ql);

f_ql = figure();

% Magnitude
subplot(4,1,1)
plot( t_ql_dn, b_ql_dmpa(4,:), t_dn, mrvector_magnitude(b_dmpa) );
title([ upper(sc) ' ' upper(instr) ' DMPA ' tstart(1:10)] );
xlabel( 'Time UTC' );
ylabel( {'|B|', '(nT)'} );
datetick();
legend('MagTeam', 'Mine');

% X-component
subplot(4,1,2)
plot( t_ql_dn, b_ql_dmpa(1,:), t_dn, b_dmpa(1,:) );
xlabel( 'Time UTC' );
ylabel( {'B_{X}', '(nT)'} );
datetick();

% Y-component
subplot(4,1,3)
plot( t_ql_dn, b_ql_dmpa(2,:), t_dn, b_dmpa(2,:) );
xlabel( 'Time UTC' );
ylabel( {'B_{Y}', '(nT)'} );
datetick();

% Z-component
subplot(4,1,4)
plot( t_ql_dn, b_ql_dmpa(3,:), t_dn, b_dmpa(3,:) );
xlabel( 'Time UTC' );
ylabel( {'B_{Z}', '(nT)'} );
datetick();