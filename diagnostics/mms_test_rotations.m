%
% Name
%   mms_test_rotations
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
%   2015-07-22      Updated to use local SDC mirror.
%
%***************************************************************************

%------------------------------------%
% Inputs                             %
%------------------------------------%
sc          = 'mms3';
instr       = 'dfg';
tstart      = '2015-08-15T13:00:00';
tend        = '2015-08-15T14:00:00';

sdc_root = '/nfs';
att_dir  = fullfile(sdc_root, 'ancillary', sc, 'defatt');
hk_dir   = fullfile(sdc_root, 'hk');


%------------------------------------%
% Find Files                         %
%------------------------------------%
	% FGM L2Pre Data File
	[l2pre_fname, count, str] = mms_file_search(sc, instr, 'srvy', 'l2pre', ...
	                                            'TStart',    tstart, ...
	                                            'TEnd',      tend);
	assert(count > 0, ['FG file not found: "' str '".']);
	
	% DSS Data File
	[dss_fname, count, str] = mms_file_search(sc, 'fields', 'hk', 'l1b', ...
	                                          'SDCroot', hk_dir, ...
	                                          'OptDesc', '101', ...
	                                          'TStart',  tstart, ...
	                                          'TEnd',    tend);
	assert(count > 0, ['FGM file not found: "' str '".']);
	
	% Attitude
	att_ftest = fullfile( att_dir, [upper(sc) '_DEFATT_%Y%D_%Y%D.V*'] );
	[defatt_files, count] = MrFile_Search(att_ftest, ...
	                                      'Closest',      true, ...
	                                      'TimeOrder',    '%Y%D', ...
	                                      'TStart',       tstart, ...
	                                      'TEnd',         tend, ...
	                                      'VersionRegex', 'V[0-9]{2}');
	assert(count > 0, ['No definitive attitude file found: "' att_ftest '".']);

%------------------------------------%
% Read and Despin                    %
%------------------------------------%

% Read the data
fgm_l2pre     = mms_fg_read_l2pre(l2pre_fname, tstart, tend);
dss           = mms_dss_read_sunpulse(dss_fname, tstart, tend, 'UniquePulse', true);
[defatt, hdr] = mms_fdoa_read_defatt(defatt_files, tstart, tend);
zmpa          = hdr.zMPA(:,1);

% Get the magnetic field
t      = fgm_l2pre.tt2000;
b_bcs  = fgm_l2pre.b_bcs(1:3,:);
b_dmpa = fgm_l2pre.b_dmpa(1:3,:);
b_gse  = fgm_l2pre.b_gse(1:3,:);
b_gsm  = fgm_l2pre.b_gsm(1:3,:);
clear fgm_l2pre

%------------------------------------%
% B_BCS ---> B_SMPA                  %
%------------------------------------%
dmpa2bcs   = mms_fg_xbcs2smpa(zmpa);
b_unh_smpa = mrvector_rotate(dmpa2bcs, b_bcs);

%------------------------------------%
% B_SMPA ---> B_DMPA                 %
%------------------------------------%
smpa2dmpa  = mms_dss_xdespin(dss, t);
b_unh_dmpa = mrvector_rotate(smpa2dmpa, b_unh_smpa);

%------------------------------------%
% B_DMPA ---> B_GSE                  %
%------------------------------------%

% Rotate to GEI
gei2dmpa  = mms_fdoa_xgei2despun(defatt, t, 'P');
dmpa2gei  = permute(gei2dmpa, [2,1,3]);
b_unh_gei = mrvector_rotate(dmpa2gei, b_unh_dmpa);

% Rotate to GSE
[mjd, utc] = tt2000toMJDutc(t');
xgei2gse   = gei2gse(mjd, utc);
b_unh_gse  = mrvector_rotate(xgei2gse, b_unh_gei);

%------------------------------------%
% B_GSE ---> B_GSM                   %
%------------------------------------%

% Get the IGRF coefficients
file      = '/home/argall/MATLAB/HapgoodRotations/igrf_coeffs.txt';
coef      = read_igrf_coeffs(file, 2015);
g10       = coef(1);
g11       = coef(2);
h11       = coef(3);
xgse2gsm  = gse2gsm(g10, g11, h11, mjd, utc);
b_unh_gsm = mrvector_rotate(xgse2gsm, b_unh_gse);

%------------------------------------%
% DMPA Results                       %
%------------------------------------%
% Convert time to datenumber to use the datetick function.
t_dn     = MrCDF_epoch2datenum(t');
f_dmpa = figure();

% X-component
subplot(3,1,1)
plot( t_dn, b_dmpa(1,:), t_dn, b_unh_dmpa(1,:) );
title([ upper(sc) ' ' upper(instr) ' DMPA ' tstart(1:10)] );
xlabel( 'Time UTC' );
ylabel( {'Bx', '(nT)'} );
datetick();
legend('FGM', 'UNH');

% Y-component
subplot(3,1,2)
plot( t_dn, b_dmpa(2,:), t_dn, b_unh_dmpa(2,:) );
xlabel( 'Time UTC' );
ylabel( {'By', '(nT)'} );
datetick();

% Z-component
subplot(3,1,3)
plot( t_dn, b_dmpa(3,:), t_dn, b_unh_dmpa(3,:) );
xlabel( 'Time UTC' );
ylabel( {'Bz', '(nT)'} );
datetick();

%------------------------------------%
% GSE Results                        %
%------------------------------------%
% Convert time to datenumber to use the datetick function.
t_dn     = MrCDF_epoch2datenum(t');
f_gse = figure();

% X-component
subplot(3,1,1)
plot( t_dn, b_gse(1,:), t_dn, b_unh_gse(1,:) );
title([ upper(sc) ' ' upper(instr) ' GSE ' tstart(1:10)] );
xlabel( 'Time UTC' );
ylabel( {'Bx', '(nT)'} );
datetick();
legend('FGM', 'UNH');

% Y-component
subplot(3,1,2)
plot( t_dn, b_gse(2,:), t_dn, b_unh_gse(2,:) );
xlabel( 'Time UTC' );
ylabel( {'By', '(nT)'} );
datetick();

% Z-component
subplot(3,1,3)
plot( t_dn, b_gse(3,:), t_dn, b_unh_gse(3,:) );
xlabel( 'Time UTC' );
ylabel( {'Bz', '(nT)'} );
datetick();

%------------------------------------%
% GSM Results                        %
%------------------------------------%
% Convert time to datenumber to use the datetick function.
t_dn     = MrCDF_epoch2datenum(t');
f_gsm = figure();

% X-component
subplot(3,1,1)
plot( t_dn, b_gsm(1,:), t_dn, b_unh_gsm(1,:) );
title([ upper(sc) ' ' upper(instr) ' GSM ' tstart(1:10)] );
xlabel( 'Time UTC' );
ylabel( {'Bx', '(nT)'} );
datetick();
legend('FGM', 'UNH');

% Y-component
subplot(3,1,2)
plot( t_dn, b_gsm(2,:), t_dn, b_unh_gsm(2,:) );
xlabel( 'Time UTC' );
ylabel( {'By', '(nT)'} );
datetick();

% Z-component
subplot(3,1,3)
plot( t_dn, b_gsm(3,:), t_dn, b_unh_gsm(3,:) );
xlabel( 'Time UTC' );
ylabel( {'Bz', '(nT)'} );
datetick();