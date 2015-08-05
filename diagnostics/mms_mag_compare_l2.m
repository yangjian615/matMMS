%
% Name
%   mms_mag_view
%
% Purpose
%   Take l1a searchcoil magnetometer data, apply calibration
%   parameters, despin, and rotate into GSE. Then, plot the 
%   results of each stage.
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-14      Written by Matthew Argall
%

get_data = true;

if get_data
%------------------------------------%
% Inputs                             %
%------------------------------------%
	sc         = 'mms2';
	tstart     = '2015-06-30T14:40:00';
	tend       = '2015-06-30T15:00:00';
	
	% Directories
	sc_cal_dir = fullfile('/home', 'argall', 'data', 'mms', 'scm_cal');
	fg_cal_dir = fullfile('/home', 'argall', 'data', 'mms', 'fg_cal');
	hk_dir     = '/nfs/hk/';
	att_dir    = fullfile('/nfs', 'ancillary', sc, 'defatt');

	% Instrument parameters
	dfg_instr   = 'dfg';
	dfg_mode    = 'f128';
	dfg_optdesc = '';

	afg_instr   = 'afg';
	afg_mode    = 'fast';
	afg_optdesc = '';
	
	scm_instr   = 'scm';
	scm_mode    = 'comm';
	scm_optdesc = 'sc256';
	duration    = 64.0;

%------------------------------------%
% Find Files                         %
%------------------------------------%
	% DFG L1A Data File
	[dfg_fname, count, str] = mms_file_search(sc, dfg_instr, dfg_mode, 'l1a', ...
	                                          'TStart',    tstart, ...
	                                          'TEnd',      tend);
	assert(count > 0, ['DFG file not found: "' str '".']);
		
	% DFG Hi-Cal File
	[dfg_hical_fname, count, str] = mms_file_search(sc, dfg_instr, 'hirangecal', 'l2pre', ...
	                                                'SDCroot', fg_cal_dir, ...
	                                                'SubDirs', '', ...
	                                                'TStart',  tstart, ...
	                                                'TEnd',    tend);
	assert(count > 0, ['DFG Hi-Cal file not found: "' str '".']);
	
	% FG Lo-Cal File
	[dfg_local_fname, count, str] = mms_file_search(sc, dfg_instr, 'lorangecal', 'l2pre', ...
	                                                'SDCroot', fg_cal_dir, ...
	                                                'SubDirs', '', ...
	                                                'TStart',  tstart, ...
	                                                'TEnd',    tend);
	assert(count > 0, ['DFG Lo-Cal file not found: "' str '".']);
	
	% AFG L1A Data File
	[afg_fname, count, str] = mms_file_search(sc, afg_instr, afg_mode, 'l1a', ...
	                                          'TStart',    tstart, ...
	                                          'TEnd',      tend);
	assert(count > 0, ['DFG file not found: "' str '".']);
		
	% AFG Hi-Cal File
	[afg_hical_fname, count, str] = mms_file_search(sc, afg_instr, 'hirangecal', 'l2pre', ...
	                                                'SDCroot', fg_cal_dir, ...
	                                                'SubDirs', '', ...
	                                                'TStart',  tstart, ...
	                                                'TEnd',    tend);
	assert(count > 0, ['AFG Hi-Cal file not found: "' str '".']);
	
	% FG Lo-Cal File
	[afg_local_fname, count, str] = mms_file_search(sc, afg_instr, 'lorangecal', 'l2pre', ...
	                                                'SDCroot', fg_cal_dir, ...
	                                                'SubDirs', '', ...
	                                                'TStart',  tstart, ...
	                                                'TEnd',    tend);
	assert(count > 0, ['AFG Lo-Cal file not found: "' str '".']);
	
	% HK 0x10e files
	hk_fname = '';
	
	% SCM L1A Data File
	[scm_fname, count, str] = mms_file_search(sc, scm_instr, scm_mode, 'l1a', ...
	                                          'TStart',    tstart, ...
	                                          'TEnd',      tend, ...
	                                          'OptDesc',   scm_optdesc);
	assert(count > 0, ['SCM L1A file not found: "' str '".']);
	
	% SCM Cal File
	scm_cal_fname = fullfile(sc_cal_dir, [sc '_' scm_instr sc(4) '_caltab_%Y%M%d%H%m%S_v*.txt']);
	[scm_cal_fname, nFiles] = MrFile_Search( scm_cal_fname, 'VersionRegex', '([0-9])');
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
% Attitude                           %
%------------------------------------%
	[attitude, att_hdr] = mms_fdoa_read_defatt(defatt_files, tstart, tend);

%------------------------------------%
% AFG                                %
%------------------------------------%
	% Get data
	[t_afg, ~, ~, b_afg_dmpa] ...
		= mms_fg_create_l2(afg_fname, afg_hical_fname, afg_local_fname, tstart, tend, ...
		                   'Attitude', attitude, ...
		                   'hk_file',  hk_fname);

%------------------------------------%
% DFG                                %
%------------------------------------%
	% Get data
	[t_dfg, ~, ~, b_dfg_dmpa] ...
		= mms_fg_create_l2(dfg_fname, dfg_hical_fname, dfg_local_fname, tstart, tend, ...
		                   'Attitude', attitude, ...
		                   'hk_file',  hk_fname);

%------------------------------------%
% SCM                                %
%------------------------------------%
	% Instrument and modea
	instr = 'scm';
	mode  = 'comm';

	% Get the data
	[t_scm, ~, ~, b_scm_dmpa] ...
		= mms_sc_create_l2(scm_fname, scm_cal_fname, tstart, tend, ...
		                   'Attitude', attitude, ...
		                   'Duration', duration, ...
		                   'zMPA',     att_hdr.zMPA(:,end));

%------------------------------------%
% Plot                               %
%------------------------------------%
	
	% Convert time to datenum
	t_dfg_datnum = MrCDF_epoch2datenum(t_dfg);
	t_afg_datnum = MrCDF_epoch2datenum(t_afg);
	t_scm_datnum = MrCDF_epoch2datenum(t_scm);

	% Convert to seconds for time range
	t_dfg_sse = MrCDF_epoch2sse(t_dfg);
	t_afg_sse = MrCDF_epoch2sse(t_afg, t_dfg(1));
	t_scm_sse = MrCDF_epoch2sse(t_scm, t_dfg(1));
end

%------------------------------------%
% Plot DMPA                          %
%------------------------------------%

%
% DMPA
%
f_bcs = figure();

% X-component
subplot(3,1,1)
plot( t_dfg_sse, b_dfg_dmpa(1,:), ...
      t_afg_sse, b_afg_dmpa(1,:), ...
      t_scm_sse, b_scm_dmpa(1,:) );
title([ 'Magnetometers in DMPA '] );
xlabel( 'Time UTC' );
ylabel( {'B_{X}', '(nT)'} );
legend('B_{DFG}', 'B_{AFG}', 'B_{SCM}');

% Y-component
subplot(3,1,2)
plot( t_dfg_sse, b_dfg_dmpa(2,:), ...
      t_afg_sse, b_afg_dmpa(2,:), ...
      t_scm_sse, b_scm_dmpa(2,:));
title([ 'Magnetometers in DMPA '] );
xlabel( 'Time UTC' );
ylabel( {'B_{Y}', '(nT)'} );

% Z-component
subplot(3,1,3)
plot( t_dfg_sse, b_dfg_dmpa(3,:), ...
      t_afg_sse, b_afg_dmpa(3,:), ...
      t_scm_sse, b_scm_dmpa(3,:) );
title([ 'Magnetometers in DMPA '] );
xlabel( 'Time UTC' );
ylabel( {'B_{Z}', '(nT)'} );
