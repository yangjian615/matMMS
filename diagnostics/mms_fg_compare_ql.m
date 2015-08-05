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
clear t b_dmpa fg_ql t_unh t_mag sunpulse attitude

%------------------------------------%
% Inputs                             %
%------------------------------------%
sc         = 'mms3';
instr      = 'dfg';
mode       = 'f128';
tstart     = '2015-06-30T05:00:00';
tend       = '2015-06-30T05:03:00';
fg_cal_dir = '/home/argall/data/mms/fg_cal/';
att_dir    = fullfile('/nfs', 'ancillary', sc, 'defatt');

%------------------------------------%
% Find Files                         %
%------------------------------------%
	% FGM L1A
	[fname_l1a, count, fsrch] = mms_file_search(sc, instr, mode, 'l1a', ...
	                                            'TStart', tstart, ...
	                                            'TEnd',   tend);
	assert(count > 0, ['No DFG file found: "' fsrch '".']);
	
	% FGM QL
	[fname_ql, count, fsrch] = mms_file_search(sc, instr, 'srvy', 'ql', ...
	                                           'TStart', tstart, ...
	                                           'TEnd',   tend);
	assert(count > 0, ['No DFG file found: "' fsrch '".']);
	
	% HI-CAL
	[hiCal_file, count, fsrch] = mms_file_search(sc, instr, 'hirangecal', 'l2pre', ...
	                                             'SDCroot', fg_cal_dir, ...
	                                             'SubDirs', '', ...
	                                             'TStart',  tstart, ...
	                                             'TEnd',    tend);
	assert(count > 0, ['No HiCal file found: "' fsrch '".']);
	
	% LO-CAL
	[loCal_file, count, fsrch] = mms_file_search(sc, instr, 'lorangecal', 'l2pre', ...
	                                             'SDCroot', fg_cal_dir, ...
	                                             'SubDirs', '', ...
	                                             'TStart',  tstart, ...
	                                             'TEnd',    tend);
	assert(count > 0, ['No LoCal file found: "' fsrch '".']);
	
	% Attitude files
	att_ftest = fullfile( att_dir, [upper(sc) '_DEFATT_%Y%D_%Y%D.V*'] );
	[defatt_files, count] = MrFile_Search(att_ftest, ...
	                                      'Closest',      true, ...
	                                      'TimeOrder',    '%Y%D', ...
	                                      'TStart',       tstart, ...
	                                      'TEnd',         tend, ...
	                                      'VersionRegex', 'V[0-9]{2}');
	assert(count > 0, ['No definitive attitude file found: "' att_ftest '".']);
	
	% Sunpulse files
	[dss_files, count, fsrch] = mms_file_search(sc, 'fields', 'hk', 'l1b', ...
	                                            'OptDesc',      '101',  ...
	                                            'SDCroot',      '/nfs/hk/', ...
	                                            'TStart',       tstart, ...
	                                            'TEnd',         tend);
	assert(count > 0, ['No sun sensor file found: "' fsrch '".']);

%------------------------------------%
% Get Data                           %
%------------------------------------%

% Read attitude and sunpulse data
attitude = mms_fdoa_read_defatt(defatt_files, tstart, tend);
sunpulse = mms_dss_read_sunpulse(dss_files, tstart, tend, 'UniquePulse', true);

% Create QL data from L1A
[t, ~, ~, b_dmpa] = mms_fg_create_l2(fname_l1a, hiCal_file, loCal_file, tstart, tend, ...
                                     'Attitude', attitude, ...
                                     'SunPulse', sunpulse);

% Read official QL data
fg_ql = mms_fg_read_ql(fname_ql, tstart, tend);

%------------------------------------%
% Plot the Results                   %
%------------------------------------%
% Convert time to datenumber to use the datetick function.
t_unh = MrCDF_epoch2datenum(t);
t_mag = MrCDF_epoch2datenum(fg_ql.tt2000);

f_ql = figure();

% Magnitude
subplot(4,1,1)
plot( t_mag, fg_ql.b_dmpa(4,:), t_unh, mrvector_magnitude(b_dmpa) );
title([ upper(sc) ' ' upper(instr) ' DMPA ' tstart(1:10)] );
xlabel( 'Time UTC' );
ylabel( {'|B|', '(nT)'} );
datetick();
legend('FGM', 'UNH');

% X-component
subplot(4,1,2)
plot( t_mag, fg_ql.b_dmpa(1,:), t_unh, b_dmpa(1,:) );
xlabel( 'Time UTC' );
ylabel( {'B_{X}', '(nT)'} );
datetick();

% Y-component
subplot(4,1,3)
plot( t_mag, fg_ql.b_dmpa(2,:), t_unh, b_dmpa(2,:) );
xlabel( 'Time UTC' );
ylabel( {'B_{Y}', '(nT)'} );
datetick();

% Z-component
subplot(4,1,4)
plot( t_mag, fg_ql.b_dmpa(3,:), t_unh, b_dmpa(3,:) );
xlabel( 'Time UTC' );
ylabel( {'B_{Z}', '(nT)'} );
datetick();