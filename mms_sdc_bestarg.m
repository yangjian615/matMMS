%
% Name
%   mms_sdc_bestarg
%
% Purpose
%   Create a MATLAB save file of inputs needed for Bestarg.
%
% Calling Sequence
%   SAVE_FILE = mms_sdc_bestarg()
%     Read and process all data required by Bestarg and save them to a
%     mat file named SAVE_FILE.
%
% Parameters
%
% Returns
%   SAVE_FILE       out, required, type=string
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-06-03      Written by Matthew Argall
%
function save_file = mms_sdc_bestarg()

% MMS2: May 9, 2015  16:08 - 16:13
% MMS4: May 6, 2015  15:30 - 15:35

	% Inputs
	sc           = 'mms4';
	tstart       = '2015-05-06T15:30:00';
	tend         = '2015-05-06T15:35:00';
	sdc_root     = '/nfs/mmsa/sdc/';
	save_dir     = '/home/argall/data/mms/matfiles/';
	attitude_dir = fullfile(sdc_root, 'ancillary', sc, 'defatt');
	hk_root      = fullfile(sdc_root, 'hk');

%------------------------------------%
% Find Files                         %
%------------------------------------%

	% FG L1B Data File
	instr   = 'dfg';
	mode    = 'srvy';
	level   = 'ql';
	optdesc = '';
	[fg_file, count, str] = mms_file_search(sc, instr, mode, level, ...
	                                        'TStart',    tstart, ...
	                                        'TEnd',      tend, ...
	                                        'OptDesc',   optdesc, ...
	                                        'SDCroot',   sdc_root);
	assert(count > 0, ['DFG file not found: "' str '".']);

	% Digital Sun Sensor L1B Data File
	instr   = 'fields';
	mode    = 'hk';
	level   = 'l1b';
	optdesc = '101';
	[dss_file, count, str] = mms_file_search(sc, instr, mode, level, ...
	                                         'TStart',    tstart, ...
	                                         'TEnd',      tend, ...
	                                         'OptDesc',   optdesc, ...
	                                         'SDCroot',   hk_root);
	assert(count > 0, ['DFG file not found: "' str '".']);

	% Attitude file
	str = fullfile( attitude_dir, [upper(sc) '_DEFATT_%Y%D_%Y%D.V*'] );
	[att_file, count] = MrFile_Search( str, ...
	                                   'Closest',      true, ...
	                                   'TStart',       tstart, ...
	                                   'TEnd',         tend, ...
	                                   'TimeOrder',    '%Y%D', ...
	                                   'VersionRegex', 'V([0-9]{2})' );
	assert(count > 0, ['Attitude file not found: "' str '".']);
	
	
	% EDP QL Data File
	%    - Find last, so file descriptors are saved.
	%    - mms4_edp_comm_ql_dce2d_20150506120000_v0.0.0.cdf
	instr   = 'edp';
	mode    = 'comm';
	level   = 'ql';
	optdesc = 'dce2d';
	[edp_file, count, str] = mms_file_search(sc, instr, mode, level, ...
	                                         'TStart',    tstart, ...
	                                         'TEnd',      tend, ...
	                                         'TimeOrder', '%Y%M%d%H%m%S', ...
	                                         'OptDesc',   optdesc, ...
	                                         'SDCroot',   sdc_root);
	assert(count > 0, ['EDP file not found: "' str '".']);
	
	
	% EDI L1A E-Field Data File
	%    - Find last, so file descriptors are saved.
	instr   = 'edi';
	mode    = 'slow';
	level   = 'l1a';
	optdesc = 'efield';
	[edi_file, count, str] = mms_file_search(sc, instr, mode, level, ...
	                                         'TStart',    tstart, ...
	                                         'TEnd',      tend, ...
	                                         'OptDesc',   optdesc, ...
	                                         'SDCroot',   sdc_root);
	assert(count > 0, ['EDI file not found: "' str '".']);

%------------------------------------%
% Read Data                          %
%------------------------------------%
	% Attitude data
	[defatt, att_hdr] = mms_fdoa_read_defatt(att_file, tstart, tend);
	
	% Sunpulse data
	sunpulse = mms_dss_read_sunpulse(dss_file, tstart, tend, 'UniquePulse', true);

	% EDI data
	edi = mms_edi_gse( edi_file, tstart, tend, ...
	                   'CS_GSE',   false, ...
	                   'CS_DMPA',  true, ...
	                   'Attitude', defatt, ...
	                   'Sunpulse', sunpulse, ...
	                   'zMPA',     att_hdr.zMPA(:,1)' );
	
	% FGM data
	fg_ql = mms_fg_read_ql(fg_file, tstart, tend);
	
	% EDP data
	edp_ql = mms_edp_read_ql(edp_file, tstart, tend);

%------------------------------------%
% Prepare Data                       %
%------------------------------------%
	% Compute the averaged magnetic field
	b_avg = mms_edi_bavg(fg_ql.tt2000, fg_ql.b_dmpa, edi.epoch_gd12, edi.epoch_gd21);
	
	
	% Create a mat file
	save_file = mms_construct_filename(sc, instr, mode, level, ....
	                                   'Directory', save_dir, ...
	                                   'OptDesc',   optdesc, ...
	                                   'TStart',    [tstart(1:4) tstart(6:7) tstart(9:10)], ...
	                                   'Version',   'v0.1.0');
	save_file(end-2:end) = 'mat';
	save(save_file, 'b_avg', 'edi', 'fg_ql', 'edp_ql');
end