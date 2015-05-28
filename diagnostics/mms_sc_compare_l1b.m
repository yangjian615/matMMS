%
% Name
%   mms_sc_view
%
% Purpose
%   Compare my level 1b data results to those of the mag team to
%   see if I am applying the calibration parameters correctly.
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-05-07      Written by Matthew Argall
%
%***************************************************************************

get_data = true;

if get_data

%------------------------------------%
% Inputs                             %
%------------------------------------%
	sc_data_dir  = '/Users/argall/Documents/Work/Data/MMS/SCM/';
	sc_cal_dir   = '/Users/argall/Documents/Work/Data/MMS/SCM_Cal/';
	sunpulse_dir = '/Users/argall/Documents/Work/Data/MMS/HK/';
	att_dir      = '/Users/argall/Documents/Work/Data/MMS/Attitude/';
	sc           = 'mms2';
	instr        = 'scm';
	mode         = 'comm';
	optdesc      = 'sc256';
	duration     = 32.0;
	tstart       = '2015-03-17T14:40:00';
	tend         = '2015-03-17T15:00:00';

%------------------------------------%
% Find Files                         %
%------------------------------------%
	% SCM L1A Data File
	[l1a_fname, count, str] = mms_file_search(sc, instr, mode, 'l1a', ...
	                                          'TStart',    tstart, ...
	                                          'TEnd',      tend, ...
	                                          'OptDesc',   optdesc, ...
	                                          'Directory', sc_data_dir);
	assert(count > 0, ['SCM file not found: "' str '".']);
	
	% SCM L1B Data File
	[l1b_fname, count, str] = mms_file_search(sc, instr, 'comm', 'l1b', ...
	                                          'TStart',    tstart, ...
	                                          'TEnd',      tend, ...
	                                          'OptDesc',   optdesc, ...
	                                          'Directory', sc_data_dir);
	assert(count > 0, ['SCM file not found: "' str '".']);
	
	% SCM Cal File
	cal_fname = fullfile(sc_cal_dir, [sc '_' instr sc(4) '_caltab_%Y%M%d%H%m%S_v*.txt']);
	[cal_fname, nFiles] = MrFile_Search( cal_fname, 'VersionRegex', '([0-9])');
	assert(nFiles > 0, ['SCM cal file not found: "' str '".']);

%------------------------------------%
% Calibrated Mag in BCS              %
%------------------------------------%
	[t, ~, b_omb] = mms_sc_bcs(l1a_fname, cal_fname, tstart, tend, duration);
	
%------------------------------------%
% Official L1B Data                  %
%------------------------------------%

	%
	% Read Data
	%
	sc_l1b = mms_sc_read_l1b(l1b_fname, tstart, tend);
	


b_omb_orig = b_omb;
b_l1b_orig = sc_l1b.b_omb;
end

% Return to original data
b_123     = b_omb_orig;
b_l1b_123 = b_l1b_orig;

%------------------------------------%
% Plot the Results                   %
%------------------------------------%
% Convert time to datenumber to use the datetick function.
t_dn_orig  = MrCDF_epoch2datenum(t);
t_l1b_orig = MrCDF_epoch2datenum(sc_l1b.tt2000);

% Smaller time interval
t0 = datenum([2015 03 17 14 40 00]);
t1 = datenum([2015 03 17 14 45 00]);

% Find the time itnerval
it  = find(t_dn_orig  >= t0 & t_dn_orig  <= t1);
itb = find(t_l1b_orig >= t0 & t_l1b_orig <= t1);

% Pick out the data
t_dn = t_dn_orig(it);
b_123 = b_123(:, it);
t_l1b = t_l1b_orig(itb);
b_l1b_123 = b_l1b_123(:, itb);


yrange = [-150, 150];

f_dmpa = figure();

% Magnitude
subplot(4,1,1)
plot( t_l1b, mrvector_magnitude(b_l1b_123), t_dn, mrvector_magnitude(b_123) );
title([ upper(sc) ' ' upper(instr) ' 123 ' tstart(1:10) ]);
xlabel( 'Time UTC' );
ylabel( {'|B|', '(nT)'} );
ylim(yrange);
datetick();
legend('MagTeam', 'Mine');

% X-component
subplot(4,1,2)
plot( t_l1b, b_l1b_123(1,:), t_dn, b_123(1,:) );
xlabel( 'Time UTC' );
ylabel( {'B_{X}', '(nT)'} );
ylim(yrange);
datetick();

% Y-component
subplot(4,1,3)
plot( t_l1b, b_l1b_123(2,:), t_dn, b_123(2,:) );
xlabel( 'Time UTC' );
ylabel( {'B_{Y}', '(nT)'} );
ylim(yrange);
datetick();

% Z-component
subplot(4,1,4)
plot( t_l1b, b_l1b_123(3,:), t_dn, b_123(3,:) );
title('Calibrated SCM Data in 123');
xlabel( 'Time UTC' );
ylabel( {'B_{Z}', '(nT)'} );
ylim(yrange);
datetick();