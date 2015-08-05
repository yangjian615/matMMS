%
% Name
%   mms_adsm_compare_l1b
%
% Purpose
%   Compare the AFG, DFG, and SCM L1B data.
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-07-31      Written by Matthew Argall
%

%------------------------------------%
% Inputs                             %
%------------------------------------%
sc         = 'mms2';
tstart     = '2015-06-30T14:40:00';
tend       = '2015-06-30T15:00:00';

% Instrument parameters
dfg_instr   = 'dfg';
dfg_mode    = 'srvy';
dfg_optdesc = '';

afg_instr   = 'afg';
afg_mode    = 'srvy';
afg_optdesc = '';

scm_instr   = 'scm';
scm_mode    = 'comm';
scm_optdesc = 'sc256';

%------------------------------------%
% Find Files                         %
%------------------------------------%
% DFG L1A Data File
[dfg_fname, count, str] = mms_file_search(sc, dfg_instr, dfg_mode, 'l1b', ...
                                          'TStart',    tstart, ...
                                          'TEnd',      tend);
assert(count > 0, ['DFG file not found: "' str '".']);
	
% AFG L1A Data File
[afg_fname, count, str] = mms_file_search(sc, afg_instr, afg_mode, 'l1b', ...
                                          'TStart',    tstart, ...
                                          'TEnd',      tend);
assert(count > 0, ['DFG file not found: "' str '".']);

% SCM L1B Data File
[scm_fname, count, str] = mms_file_search(sc, scm_instr, scm_mode, 'l1b', ...
                                          'TStart',    tstart, ...
                                          'TEnd',      tend, ...
                                          'OptDesc',   scm_optdesc);
assert(count > 0, ['SCM L1B file not found: "' str '".']);

%------------------------------------%
% Read Data                          %
%------------------------------------%
% Get data
afg_l1b = mms_fg_read_l1b(afg_fname, tstart, tend);
dfg_l1b = mms_fg_read_l1b(dfg_fname, tstart, tend);
scm_l1b = mms_sc_read_l1b(scm_fname, tstart, tend);

%------------------------------------%
% Plot                               %
%------------------------------------%

% Convert to seconds for time range
t_dfg_sse = MrCDF_epoch2sse(afg_l1b.tt2000);
t_afg_sse = MrCDF_epoch2sse(dfg_l1b.tt2000, afg_l1b.tt2000(1));
t_scm_sse = MrCDF_epoch2sse(scm_l1b.tt2000, afg_l1b.tt2000(1));
keyboard

% Create the figure
f_omb = figure();

% X-component
subplot(3,1,1)
plot( t_dfg_sse, afg_l1b.b_omb(1,:), ...
      t_afg_sse, dfg_l1b.b_omb(1,:), ...
      t_scm_sse, scm_l1b.b_omb(1,:) );
title([ 'Magnetometers in OMB'] );
xlabel( 'Time UTC' );
ylabel( {'B_{X}', '(nT)'} );
legend('B_{AFG}', 'B_{DFG}', 'B_{SCM}');

% Y-component
subplot(3,1,2)
plot( t_dfg_sse, afg_l1b.b_omb(2,:), ...
      t_afg_sse, dfg_l1b.b_omb(2,:), ...
      t_scm_sse, scm_l1b.b_omb(2,:));
title([ 'Magnetometers in OMB'] );
xlabel( 'Time UTC' );
ylabel( {'B_{Y}', '(nT)'} );

% Z-component
subplot(3,1,3)
plot( t_dfg_sse, afg_l1b.b_omb(3,:), ...
      t_afg_sse, dfg_l1b.b_omb(3,:), ...
      t_scm_sse, scm_l1b.b_omb(3,:) );
title([ 'Magnetometers in OMB'] );
xlabel( 'Time UTC' );
ylabel( {'B_{Z}', '(nT)'} );