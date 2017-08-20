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

%------------------------------------%
% Inputs                             %
%------------------------------------%
tf_read = true;
sc      = 'mms4';
instr   = 'dfg-scm';
mode    = 'brst';
level   = 'l2plus';
optdesc = 'fsm-split';
tstart  = '2015-10-16T13:05:00';
tend    = '2015-10-16T13:10:00';
fsm_dir = '/nfs/fsm/temp';

%------------------------------------%
% Find the FSM & FGM File            %
%------------------------------------%
if tf_read
	[fname_fgm, count, srchstr] = mms_file_search(sc, 'dfg', 'srvy', 'l2pre', ...
	                                             'TStart', tstart,           ...
	                                             'TEnd',   tend);
	assert(count > 0, ['Unable to find FGM file: "' srchstr '".']);
	
	
	% Get the data
	fsm_l2plus = mms_fsm_l2plus_read(sc, 'srvy', tstart, tend);
	fgm_l2pre  = mms_fg_read_l2pre(fname_fgm,  tstart, tend);
	
	%------------------------------------%
	% DMPA Results                       %
	%------------------------------------%
	% Convert time to datenumber to use the datetick function.
	t_fsm = MrCDF_epoch2datenum(fsm_l2plus.tt2000);
	t_fgm = MrCDF_epoch2datenum(fgm_l2pre.tt2000);
end
	
fig_dmpa = figure();

% Magnitude
subplot(4,1,1)
plot( t_fsm, mrvector_magnitude(fsm_l2plus.b_dmpa), ...
      t_fgm, fgm_l2pre.b_dmpa(4,:) );
title([ upper(sc) ' ' upper(instr) ' DMPA ' tstart(1:10)] );
xlabel( 'Time UTC' );
ylabel( {'|B|', '(nT)'} );
datetick();
legend('FSM', 'DFG');

% X-component
subplot(4,1,2)
plot( t_fsm, fsm_l2plus.b_dmpa(1,:),     ...
      t_fgm, fgm_l2pre.b_dmpa(1,:) );
xlabel( 'Time UTC' );
ylabel( {'B_{X}', '(nT)'} );
datetick();

% Y-component
subplot(4,1,3)
plot( t_fsm, fsm_l2plus.b_dmpa(2,:), ...
      t_fgm, fgm_l2pre.b_dmpa(2,:) );
xlabel( 'Time UTC' );
ylabel( {'B_{Y}', '(nT)'} );
datetick();

% Z-component
subplot(4,1,4)
plot( t_fsm, fsm_l2plus.b_dmpa(3,:),     ...
      t_fgm, fgm_l2pre.b_dmpa(3,:) );
xlabel( 'Time UTC' );
ylabel( {'B_{Z}', '(nT)'} );
datetick();



%------------------------------------%
% GSE Results                        %
%------------------------------------%

fig_gse = figure();

% Magnitude
subplot(4,1,1)
plot( t_fsm, mrvector_magnitude(fsm_l2plus.b_gse), ...
      t_fgm, fgm_l2pre.b_gse(4,:) );
title([ upper(sc) ' ' upper(instr) ' GSE ' tstart(1:10)] );
xlabel( 'Time UTC' );
ylabel( {'|B|', '(nT)'} );
datetick();
legend('FSM', 'DFG');

% X-component
subplot(4,1,2)
plot( t_fsm, fsm_l2plus.b_gse(1,:),     ...
      t_fgm, fgm_l2pre.b_gse(1,:) );
xlabel( 'Time UTC' );
ylabel( {'B_{X}', '(nT)'} );
datetick();

% Y-component
subplot(4,1,3)
plot( t_fsm, fsm_l2plus.b_gse(2,:), ...
      t_fgm, fgm_l2pre.b_gse(2,:) );
xlabel( 'Time UTC' );
ylabel( {'B_{Y}', '(nT)'} );
datetick();

% Z-component
subplot(4,1,4)
plot( t_fsm, fsm_l2plus.b_gse(3,:),     ...
      t_fgm, fgm_l2pre.b_gse(3,:) );
xlabel( 'Time UTC' );
ylabel( {'B_{Z}', '(nT)'} );
datetick();