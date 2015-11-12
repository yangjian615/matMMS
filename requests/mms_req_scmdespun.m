%
% Name
%   mms_req_scmdespun
%
% Purpose
%   Request for Kristoff Paulson
%
%      Burst mode data  2015-08-04T15:40-16:35
%      Survey mode data 2015-08-15T08:30-09:30
%
%   Create AFG or DFG data in GEI so that it can be compared with RBSP mag data
%   to determine spin-axis offsets. Reads in L1A data and turns it into a L2
%   product.
%
%   This will work only on a server to which the MMS SDC is mounted to the /nfs/
%   directory.
%
% Calling Sequence
%   mms_req_fggei(SC, INSTR, MODE, TSTART, TEND)
%     Read in L1A data from spacecraft SC, instrument INSTR, and telemetry mode
%     MODE during the time interval [TSTART, TEND], turn the data into a L2
%     product in GEI coordinates, then write the time and field to an ASCII file.
%
%   FNAME = mms_req_fggei(__)
%     Return the name of the ASCII file.
%
% Parameters
%   SC              in, required, type = string
%   INSTR           in, required, type = string
%   MODE            in, required, type = string
%   TSTART          in, required, type = string
%   TEND            in, required, type = string
%   OUTDIR          in, optional, type = string, default = pwd()
%
% Returns
%   FNAME           out, optional, type=string
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-05-15      Written by Matthew Argall
%
function fname = mms_req_scmdespun(sc, instr, mode, tstart, tend, outdir)

	% Default inputs
	sc      = 'mms4';
	instr   = 'scm';
	mode    = 'srvy';
	level   = 'l1a';
	optdesc = 'scm';
	tstart  = '2015-08-15T08:30:00';
	tend    = '2015-08-15T09:30:00';
	
	% Output directory
	if nargin() < 6
		outdir = '/home/argall/data/';
	end
	
	% Constants
	attdir      = fullfile('/nfs', 'ancillary', sc, 'defatt');
	hk_dir      = fullfile('/nfs', 'hk');
	scm_cal_dir = '/home/argall/data/mms/scm_cal/';
	
%------------------------------------%
% Find Files                         %
%------------------------------------%
	% DATA
	[fname_scm, count, fsrch] = mms_file_search(sc, instr, mode, level, ...
	                                            'OptDesc', optdesc,     ...
	                                            'TStart',  tstart,      ...
	                                            'TEnd',    tend);
	assert(count > 0, ['No SCM file found: "' fsrch '".']);
	
	% SCM Cal File
	switch sc
		case 'mms1'
			optdesc = 'scm1';
		case 'mms2'
			optdesc = 'scm1';
		case 'mms3'
			optdesc = 'scm4';
		case 'mms4'
			optdesc = 'scm3';
	end
	scm_ftest = fullfile(scm_cal_dir, [sc '_' optdesc '_caltab_%Y%M%d%H%M%S_v*.txt']);
	[scm_cal_file, count] = MrFile_Search(scm_ftest, ...
	                                      'Closest',      true, ...
	                                      'TimeOrder',    '%Y%M%d%H%M%S', ...
	                                      'TStart',       tstart, ...
	                                      'TEnd',         tend, ...
	                                      'VersionRegex', 'v[0-9]');
	assert(count > 0, ['No SCM calibration file found: "' scm_ftest '".']);
	
	% HK
	[fname_hk, count, fsrch] = mms_file_search(sc, 'fields', 'hk', 'l1b', ...
	                                            'OptDesc', '101',  ...
	                                            'SDCroot', hk_dir, ...
	                                            'TStart',  tstart,  ...
	                                            'TEnd',    tend);
	assert(count > 0, ['No sunpulse file found: "' fsrch '".']);
	
	% Attitude
	ftest = fullfile(attdir, [upper(sc) '_DEFATT_%Y%D_%Y%D.V*']);
	[fname_att, count] = MrFile_Search(ftest, ...
	                                   'Closest',      true, ...
	                                   'TStart',       tstart, ...
	                                   'TEnd',         tend, ...
	                                   'TimeOrder',    '%Y%D', ...
	                                   'VersionRegex', 'V[0-9]+');
	assert(count > 0, ['No attitude file found: "' ftest '".']);
	
%------------------------------------%
% Create L2                          %
%------------------------------------%
	% Get attitude data
	[attitude, att_hdr] = mms_fdoa_read_defatt(fname_att, tstart, tend);
	zMPA = att_hdr.zMPA(:,1);
	
	% Get sunpulse data
	sunpulse = mms_dss_read_sunpulse(fname_hk, tstart, tend, 'UniquePulse', true);

	% Create GEI data
	[tt2000, b_gse, ~, b_dmpa, ~, ~, b_bcs]                       ...
		= mms_sc_create_l2(fname_scm, scm_cal_file, tstart, tend, ...
		                   'Attitude', attitude,                  ...
		                   'SunPulse', sunpulse,                  ...
		                   'zMPA',     zMPA);
	
%------------------------------------%
% Save a mat file                    %
%------------------------------------%
	% Get the start time and version of the file
	[~, ~, ~, ~, fstart, version] = mms_dissect_filename(fname_scm);

	% Create a file name
	fname = mms_construct_filename(sc, instr, mode, 'l2', ...
	                               'OptDesc',   'unh', ...
	                               'TStart',    fstart, ...
	                               'Directory', outdir, ...
	                               'Version',   version);
	fname(end-2:end) = 'mat';
	
	% Save variables
	save(fname, 'tt2000', 'b_gse', 'b_dmpa', 'b_bcs');
	
%------------------------------------%
% Write the Data                     %
%------------------------------------%
	% Create a file name
	fname = mms_construct_filename(sc, instr, mode, 'l2', ...
	                               'OptDesc',   'unh-gse', ...
	                               'TStart',    fstart, ...
	                               'Directory', outdir, ...
	                               'Version',   version);
	fname(end-2:end) = 'dat';

	% Convert time to UTC
	utc = spdfencodett2000(tt2000);

	% Opend the file
	fid = fopen(fname, 'w');
	
	% Write to the file
	fprintf(fid, '                          %s              %s              %s              %s\n', 'UTC', 'Bx', 'By', 'Bz');
	for ii = 1 : length(tt2000)
		fprintf( fid, '%s   %13.6f   %13.6f   %13.6f\n', utc{ii}, b_gse(:,ii) );
	end
	
	% Close the file
	fclose(fid);
	
	if nargout() == 0
		fprintf('File written to: %s.\n', fname)
		clear fname
	end
end