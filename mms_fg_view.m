%
% Name
%   mms_fg_view
%
% Purpose
%   Take l1a magnetometer data, apply calibration parameters,
%   despin, and rotate into GSE. Then, plot the results of
%   each stage.
%
% Calling Sequence
%   ATTITUDE = mms_fdoa_read_defatt(SC, TIME)
%     Return structure of attitude data ATTITUDE, from MMS spacecraft SC
%     (e.g., 'mms2') between the times TSTART and TEND. Attitude data
%     files will be searched for in ATT_DIR.
%
%   [ATTITUDE, ATT_HDR] = mms_fdoa_read_defatt(__)
%     Also return a structure of header information from all files read.
%
% Parameters
%   SC              in, required, type=char
%   TIME            in, required, type=char
%   'AttDir'        in, optional, type=char, default=pwd()
%                   Directory in which to find attitude data.
%   'Type'          in, optional, type=char, default='L'
%                   Which type of phase - L, P, w, z.
%   'Offset'        in, optional, type=char, default='L'
%                   Constant offset added to phase.
%   
%
% Returns
%   ATTITUDE         out, required, type=struct
%   ATT_HDR          out, required, type=struct
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-06      Written by Matthew Argall
%
%*****************************************************************************************


%------------------------------------%
% Inputs                             %
%------------------------------------%
fg_data_dir = '/Users/argall/Documents/Work/Data/MMS/DFG/';
fg_cal_dir  = '/Users/argall/Documents/Work/Data/MMS/FG_Cal/';
att_dir     = '/Users/argall/Documents/Work/Data/MMS/Attitude/';
sc          = 'mms2';
instr       = 'dfg';
mode        = 'f128';
tstart      = '2015-03-20T04:00:00';
tend        = '2015-03-20T05:00:00';


%------------------------------------%
% Calibrated Mag in BCS              %
%------------------------------------%
disp( 'Getting mag data...' )
[B_bcs, t] = mms_fg_bcs(sc, instr, mode, tstart, tend, ...
                        'DataDir', fg_data_dir, ...
                        'CalDir',  fg_cal_dir);

%------------------------------------%
% Despin                             %
%------------------------------------%
disp( 'Rotating to DSL...' )
[bcs2dsl, ~, attitude] = mms_fdoa_xdespin(sc, tstart, tend, 'L', t, att_dir);

% Rotate to DSL
B_dsl = mrvector_rotate(bcs2dsl, B_bcs);
clear bcs2dsl

%------------------------------------%
% Rotate to GSE                      %
%------------------------------------%
disp( 'Rotating to GSE...' )

% Transform into GEI
despun2gei = mms_fdoa_xgei2despun(attitude, t, 'L');
B_gei      = mrvector_rotate(despun2gei, B_dsl);

% Transform matrix GEI -> GSE
%   - Modified Julian Date (mjd).
%   - UTC seconds since midnight (ssm).
timevec = MrCDF_Epoch_Breakdown(t);
mjd     = date2mjd(timevec(1,:), timevec(2,:), timevec(3,:));
ssm = timevec(4,:) * 3600.0 + ...
      timevec(5,:) * 60.0   + ...
      timevec(6,:)          + ...
      timevec(7,:) * 1e-3   + ...
      timevec(8,:) * 1e-6   + ...
      timevec(9,:) * 1e-9;
GEI2GSE = gei2gse(mjd, ssm);

% Transform to GSE
B_gse = mrvector_rotate(GEI2GSE, B_gei);
clear despun2gei gei2gse B_gei

%------------------------------------%
% Plot Results                       %
%------------------------------------%
	
% Convert time to datenum
t_datnum = MrCDF_epoch2datenum(t);

%
% BCS
%
figure()

% X-component
subplot(3,1,1)
plot(t_datnum, B_bcs(1,:));
title([ upper(instr) ' in BCS '] );
xlabel( 'Time UTC' );
ylabel( {'B_{X}', '(nT)'} );
datetick();

% Y-component
subplot(3,1,2)
plot(t_datnum, B_bcs(2,:));
title([ upper(instr) ' in BCS '] );
xlabel( 'Time UTC' );
ylabel( {'B_{Y}', '(nT)'} );
datetick();

% Z-component
subplot(3,1,3)
plot(t_datnum, B_bcs(3,:));
title([ upper(instr) ' in BCS '] );
xlabel( 'Time UTC' );
ylabel( {'B_{Z}', '(nT)'} );
datetick();

%
% DSL
%
figure()

% X-component
subplot(3,1,1)
plot(t_datnum, B_dsl(1,:));
title([ upper(instr) ' in DSL '] );
xlabel( 'Time UTC' );
ylabel( {'B_{X}', '(nT)'} );
datetick();

% Y-component
subplot(3,1,2)
plot(t_datnum, B_dsl(2,:));
title([ upper(instr) ' in DSL '] );
xlabel( 'Time UTC' );
ylabel( {'B_{Y}', '(nT)'} );
datetick();

% Z-component
subplot(3,1,3)
plot(t_datnum, B_dsl(3,:));
title([ upper(instr) ' in DSL '] );
xlabel( 'Time UTC' );
ylabel( {'B_{Z}', '(nT)'} );
datetick();

%
% GSE
%
figure()

% X-component
subplot(3,1,1)
plot(t_datnum, B_gse(1,:));
title([ upper(instr) ' in GSE '] );
xlabel( 'Time UTC' );
ylabel( {'B_{X}', '(nT)'} );
datetick();

% Y-component
subplot(3,1,2)
plot(t_datnum, B_gse(2,:));
title([ upper(instr) ' in GSE '] );
xlabel( 'Time UTC' );
ylabel( {'B_{Y}', '(nT)'} );
datetick();

% Z-component
subplot(3,1,3)
plot(t_datnum, B_gse(3,:));
title([ upper(instr) ' in GSE '] );
xlabel( 'Time UTC' );
ylabel( {'B_{Z}', '(nT)'} );
datetick();

