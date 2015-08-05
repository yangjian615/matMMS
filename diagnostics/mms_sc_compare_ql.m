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
	sc         = 'mms2';
	instr      = 'scm';
	mode       = 'comm';
	optdesc    = 'sc256';
	duration   = 64.0;
	tstart     = '2015-06-30T14:40:00';
	tend       = '2015-06-30T15:00:00';
	sc_cal_dir = fullfile('/home', 'argall', 'data', 'mms', 'scm_cal');
	hk_dir     = '/nfs/hk/';
	att_dir    = fullfile('/nfs', 'ancillary', sc, 'defatt');

%------------------------------------%
% Find Files                         %
%------------------------------------%
	% SCM L1A Data File
	[l1a_fname, count, str] = mms_file_search(sc, instr, mode, 'l1a', ...
	                                          'TStart',    tstart, ...
	                                          'TEnd',      tend, ...
	                                          'OptDesc',   optdesc);
	assert(count > 0, ['SCM L1A file not found: "' str '".']);
	
	% SCM L1B Data File
	[l1b_fname, count, str] = mms_file_search(sc, instr, mode, 'l1b', ...
	                                          'TStart',    tstart, ...
	                                          'TEnd',      tend, ...
	                                          'OptDesc',   optdesc);
	assert(count > 0, ['SCM L1B file not found: "' str '".']);
	
	% SCM Cal File
	cal_fname = fullfile(sc_cal_dir, [sc '_' instr sc(4) '_caltab_%Y%M%d%H%m%S_v*.txt']);
	[cal_fname, nFiles] = MrFile_Search( cal_fname, 'VersionRegex', '([0-9])');
	assert(nFiles > 0, ['SCM cal file not found: "' str '".']);
	
	% Attitude files
	att_ftest = fullfile( att_dir, [upper(sc) '_DEFATT_%Y%D_%Y%D.V*'] );
	[defatt_files, count] = MrFile_Search(att_ftest, ...
	                                      'Closest',      true, ...
	                                      'TimeOrder',    '%Y%D', ...
	                                      'TStart',       tstart, ...
	                                      'TEnd',         tend, ...
	                                      'VersionRegex', 'V[0-9]{2}');
	assert(count > 0, ['No definitive attitude file found: "' att_ftest '".']);

%------------------------------------%
% Calibrated Mag in BCS              %
%------------------------------------%
	% Read attitude data
	[attitude, att_hdr]  = mms_fdoa_read_defatt(defatt_files, tstart, tend);

	% Calibrate SCM data
	[t, ~, b_dmpa] = mms_sc_create_l2(l1a_fname, cal_fname, tstart, tend, ...
	                                  'Attitude', attitude, ...
	                                  'Duration', duration, ...
	                                  'zMPA',     att_hdr.zMPA(:,end));
	
%------------------------------------%
% Official L1B Data                  %
%------------------------------------%
	% Read L1B Data
	sc_l1b = mms_sc_read_l1b(l1b_fname, tstart, tend);
	
	% Rotate to SMPA
	omb2smpa   = mms_fg_xomb2smpa();
	b_scm_smpa = mrvector_rotate( omb2smpa, sc_l1b.b_omb );
	
	% Despin
	smpa2dmpa  = mms_fdoa_xdespin(attitude, sc_l1b.tt2000, 'L');
	b_scm_dmpa = mrvector_rotate( smpa2dmpa, b_scm_smpa );
end

%------------------------------------%
% Plot the Results                   %
%------------------------------------%
% Convert time to datenumber to use the datetick function.
t_dn     = MrCDF_epoch2datenum(t);
t_scm_dn = MrCDF_epoch2datenum(sc_l1b.tt2000);


f_dmpa = figure();

% Magnitude
subplot(4,1,1)
plot( t_scm_dn, mrvector_magnitude(b_scm_dmpa), t_dn, mrvector_magnitude(b_dmpa) );
title([ upper(sc) ' ' upper(instr) ' 123 ' tstart(1:10) ]);
xlabel( 'Time UTC' );
ylabel( {'|B|', '(nT)'} );
ylim(yrange);
datetick();
legend('MagTeam', 'Mine');

% X-component
subplot(4,1,2)
plot( t_scm_dn, b_scm_dmpa(1,:), t_dn, b_dmpa(1,:) );
xlabel( 'Time UTC' );
ylabel( {'B_{X}', '(nT)'} );
ylim(yrange);
datetick();

% Y-component
subplot(4,1,3)
plot( t_scm_dn, b_scm_dmpa(2,:), t_dn, b_dmpa(2,:) );
xlabel( 'Time UTC' );
ylabel( {'B_{Y}', '(nT)'} );
ylim(yrange);
datetick();

% Z-component
subplot(4,1,4)
plot( t_scm_dn, b_scm_dmpa(3,:), t_dn, b_dmpa(3,:) );
title('Calibrated SCM Data in 123');
xlabel( 'Time UTC' );
ylabel( {'B_{Z}', '(nT)'} );
ylim(yrange);
datetick();