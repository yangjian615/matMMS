%
% Name
%   mms_fg_view
%
% Purpose
%   View calibrated, orthogonalized fluxgate data in DMPA coordinates.
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
attitude_dir = '/Users/argall/Documents/Work/Data/MMS/Attitude/';
sunpulse_dir = '/Users/argall/Documents/Work/Data/MMS/HK/';
sc           = 'mms2';
instr        = 'dfg';
tstart       = '2015-03-17T15:00:00';
tend         = '2015-03-17T16:00:00';

if strcmp(instr, 'dfg')
	fg_data_dir = dfg_data_dir;
	mode        = 'f128';
else
	fg_data_dir = afg_data_dir;
	mode        = 'fast';
end

if strcmp(instr, 'dfg')
	fg_data_dir = dfg_data_dir;
else
	fg_data_dir = afg_data_dir;
end

%------------------------------------%
% Find Files                         %
%------------------------------------%

	%
	% Data files
	%

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
	
	%
	% Cal Files
	%
	
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
	
	%
	% Definitive Attitude
	%
	
	% MMS#_DEFATT_%Y%D_%Y%D.V00
	att_tmp = fullfile(att_dir, [upper(sc) '_DEFATT_%Y%D_%Y%D.V*']);
	[att_fname, nFiles] = MrFile_Search( att_tmp, ...
	                                 'VersionRegex', 'V([0-9]{2})', ...
	                                 'TStart',       tstart, ...
	                                 'TEnd',         tend, ...
	                                 'TimeOrder',    '%Y%D' );
	% Make sure the file exists
	assert( nFiles > 0, ['Definitive attitude file not found: "' att_tmp '".']);
	
	%
	% Sunpulse
	%
	
	% HK 101
	[dss_fname, count, str] = mms_file_search(sc, 'fields', 'hk', 'l1b', ...
	                                            'OptDesc',   '101', ...
	                                            'TStart',    tstart, ...
	                                            'TEnd',      tend, ...
	                                            'Directory', sunpulse_dir);
	assert(count > 0, ['Sunpulse file not found: "' str '".']);

%------------------------------------%
% Read Attitude and Sunpulse Files   %
%------------------------------------%
	                                 
	% Read Attitude Data
	attitude = mms_fdoa_read_defatt(att_fname, tstart, tend);

	% Read sun pulse data
	sunpulse = mms_dss_read_sunpulse(dss_fname, tstart, tend, ...
	                                 'UniquePulse', true);

%------------------------------------%
% Calibrated Mag in GSE              %
%------------------------------------%
	[t, ~, b_dmpa] = mms_fg_gse(l1a_fname, hical_fname, local_fname, tstart, tend, ...
	                            'Attitude', attitude, ...
	                            'SunPulse', sunpulse);

%------------------------------------%
% Plot the Results                   %
%------------------------------------%
% Convert time to datenumber to use the datetick function.
t_dn = MrCDF_epoch2datenum(t');


fig = figure();

% Magnitude
subplot(4,1,1)
plot( t_dn, mrvector_magnitude(b_dmpa) );
title([ upper(sc) ' ' upper(instr) ' DMPA ' tstart(1:10)] );
xlabel( 'Time UTC' );
ylabel( {'|B|', '(nT)'} );
datetick();

% X-component
subplot(4,1,2)
plot( t_dn, b_dmpa(1,:) );
xlabel( 'Time UTC' );
ylabel( {'B_{X}', '(nT)'} );
datetick();

% Y-component
subplot(4,1,3)
plot( t_dn, b_dmpa(2,:) );
xlabel( 'Time UTC' );
ylabel( {'B_{Y}', '(nT)'} );
datetick();

% Z-component
subplot(4,1,4)
plot( t_dn, b_dmpa(3,:) );
xlabel( 'Time UTC' );
ylabel( {'B_{Z}', '(nT)'} );
datetick();