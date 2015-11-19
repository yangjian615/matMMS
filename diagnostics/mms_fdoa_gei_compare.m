%
% Name
%   mms_fdoa_gei_compare
%
% Purpose
%   Compare different methods for rotating from spacecraft coordinate
%   systems (e.g. BCS, DMPA, DSL) to GEI.
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-07-13      Written by Matthew Argall
%
%***************************************************************************
get_data = true;

if get_data
	%------------------------------------%
	% Inputs                             %
	%------------------------------------%
	sc      = 'mms1';
	instr   = 'dfg';
	mode    = 'srvy';
	optdesc = '';
	tstart  = '2015-08-16T00:00:00';
	tend    = '2015-08-16T24:00:00';
	att_dir = fullfile('/nfs', 'ancillary', sc, 'defatt');
	
	%------------------------------------%
	% Find the Data Files                %
	%------------------------------------%
	
	% L1B
	[fname_fg_l1b, count, srchstr] = mms_file_search(sc, instr, mode, 'l1b',   ...
	                                                 'OptDesc',   optdesc, ...
	                                                 'TStart',    tstart,  ...
	                                                 'TEnd',      tend);
	assert(count > 0, ['Unable to find FSM file: "' srchstr '".']);
	
	% QL
	[fname_fg_ql, count, srchstr] = mms_file_search(sc, instr, 'srvy', 'ql',   ...
	                                                'TStart',    tstart,  ...
	                                                'TEnd',      tend);
	assert(count > 0, ['Unable to find FSM file: "' srchstr '".']);
	
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
	% Read Data                          %
	%------------------------------------%
	
	% DFG
	fg_l1b = mms_fg_read_l1b(fname_fg_l1b, tstart, tend);
	fg_ql  = mms_fg_read_ql(fname_fg_ql,   tstart, tend);
	
	% Attitude
	attitude = mms_fdoa_read_defatt(defatt_files, tstart, tend);
end
	
%------------------------------------%
% Rotate to GEI                      %
%------------------------------------%

% Rotate from DMPA to GEI with RA and Dec
gei2dmpa   = mms_fdoa_xgei2despun( attitude, fg_ql.tt2000, 'L' );

keyboard

b_dmpa2gei = mrvector_rotate( permute(gei2dmpa, [2,1,3]), fg_ql.b_dmpa(1:3, :) );

% Rotate from BCS to GEI with quaternions
Q         = attitude.q;
% Q(1:3,:)  = mrvector_normalize( Q(1:3,:) );
% Q         = mrquaternion_invert(Q);
% Q(4,:)    = Q(4,:) * pi/180.0;
qbcs2gei  = mms_fdoa_xgei2bcs( attitude.tt2000, Q, fg_l1b.tt2000 );
b_bcs2gei = mrquaternion_rotate( fg_l1b.b_bcs(1:3,:), qbcs2gei );



%------------------------------------%
% Plot the Results                   %
%------------------------------------%
% Convert time to datenumber to use the datetick function.
t0    = min( [fg_l1b.tt2000(1) fg_ql.tt2000(1)] );
t_l1b = MrCDF_epoch2sse(fg_l1b.tt2000, t0);
t_ql  = MrCDF_epoch2sse(fg_ql.tt2000, t0);
t_att = MrCDF_epoch2sse(attitude.tt2000, t0);

% Interpolate
b_bcs = fg_l1b.b_bcs(1:3,:);

ql_fig = figure();

% Magnitude
subplot(4,1,1)
plot( t_l1b, mrvector_magnitude(b_bcs), t_l1b, mrvector_magnitude(b_bcs2gei), t_ql, mrvector_magnitude(b_dmpa2gei) );
title([ upper(sc) ' ' upper(instr) ' DMPA ' tstart(1:10)] );
xlabel( 'Time UTC' );
ylabel( {'|B|', '(nT)'} );
legend('b_{bcs}', 'bcs2gei', 'dmpa2gei');

% X-component
subplot(4,1,2)
plot( t_l1b, b_bcs(1,:), t_l1b, b_bcs2gei(1,:), t_ql, b_dmpa2gei(1,:) );
xlabel( 'Time UTC' );
ylabel( {'B_{X}', '(nT)'} );

% Y-component
subplot(4,1,3)
plot( t_l1b, b_bcs(2,:), t_l1b, b_bcs2gei(2,:), t_ql, b_dmpa2gei(2,:) );
xlabel( 'Time UTC' );
ylabel( {'B_{Y}', '(nT)'} );

% Z-component
subplot(4,1,4)
plot( t_l1b, b_bcs(3,:), t_l1b, b_bcs2gei(3,:), t_ql, b_dmpa2gei(3,:) );
xlabel( 'Time UTC' );
ylabel( {'B_{Z}', '(nT)'} );
