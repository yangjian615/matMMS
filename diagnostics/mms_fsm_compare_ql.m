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
sc      = 'mms1';
instr   = 'dfg-scm';
mode    = 'comm';
level   = 'ql';
tstart  = '2015-06-17T00:00:00';
tend    = '2015-06-17T24:00:00';
fsm_dir = '/nfs/edi/fsm/';

%------------------------------------%
% Find the FSM File                  %
%------------------------------------%
[fname_fsm, count, srchstr] = mms_file_search(sc, instr, mode, level, ...
                                              'Directory', fsm_dir, ...
                                              'TStart',    tstart, ...
                                              'TEnd',      tend);
assert(count > 0, ['Unable to find FSM file: "' srchstr '".']);

% Get the data
fsm_ql = mms_fsm_read_ql(fname_fsm);

%------------------------------------%
% Plot the Results                   %
%------------------------------------%
% Convert time to datenumber to use the datetick function.
t_fsm = MrCDF_epoch2datenum(fsm_ql.tt2000);
t_fgm = MrCDF_epoch2datenum(fsm_ql.tt2000_fgm);

ql_fig = figure();

% Magnitude
subplot(4,1,1)
plot( t_fsm, mrvector_magnitude(fsm_ql.b_dmpa), t_fgm, mrvector_magnitude(fsm_ql.b_fgm_dmpa) );
title([ upper(sc) ' ' upper(instr) ' DMPA ' tstart(1:10)] );
xlabel( 'Time UTC' );
ylabel( {'|B|', '(nT)'} );
datetick();
legend('FSM', 'DFG');

% X-component
subplot(4,1,2)
plot( t_fsm, fsm_ql.b_dmpa(1,:), t_fgm, fsm_ql.b_fgm_dmpa(1,:) );
xlabel( 'Time UTC' );
ylabel( {'B_{X}', '(nT)'} );
datetick();

% Y-component
subplot(4,1,3)
plot( t_fsm, fsm_ql.b_dmpa(2,:), t_fgm, fsm_ql.b_fgm_dmpa(2,:) );
xlabel( 'Time UTC' );
ylabel( {'B_{Y}', '(nT)'} );
datetick();

% Z-component
subplot(4,1,4)
plot( t_fsm, fsm_ql.b_dmpa(3,:), t_fgm, fsm_ql.b_fgm_dmpa(3,:) );
xlabel( 'Time UTC' );
ylabel( {'B_{Z}', '(nT)'} );
datetick();