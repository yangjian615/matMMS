%
% Name
%   mms_req_fggei
%
% Purpose
%   Request for Kristoff Paulson
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
function fname = mms_req_fggei(sc, instr, mode, tstart, tend, outdir)

	% Default inputs
	sc     = 'mms1';
	instr  = 'dfg';
	mode   = 'fast';
	level  = 'l1a';
	tstart = '2015-06-12T14:00:00';
	tend   = '2015-06-12T19:00:00';
	
	% Output directory
	if nargin() < 6
		outdir = '/home/argall/data/';
	end
	
	% Constants
	attdir  = fullfile('/nfs', 'ancillary', sc, 'defatt');
%	cal_dir = fullfile('/nfs', 'mag_cal');
	cal_dir = fullfile('/home', 'argall', 'data', 'mms', 'fg_cal');
	
%------------------------------------%
% Find Files                         %
%------------------------------------%
	% DATA
	[fname_fg, count, fsrch] = mms_file_search(sc, instr, mode, level, ...
	                                           'TStart', tstart, ...
	                                           'TEnd',   tend);
	assert(count > 0, ['No DFG file found: "' fsrch '".']);
	
	% HI-CAL
	[fname_hi, count, fsrch] = mms_file_search(sc, instr, 'hirangecal', 'l2pre', ...
	                                            'SDCroot', cal_dir, ...
	                                            'SubDirs', '', ...
	                                            'TStart', tstart, ...
	                                            'TEnd',   tend);
	assert(count > 0, ['No HiCal file found: "' fsrch '".']);
	
	% LO-CAL
	[fname_lo, count, fsrch] = mms_file_search(sc, instr, 'lorangecal', 'l2pre', ...
	                                            'SDCroot', cal_dir, ...
	                                            'SubDirs', '', ...
	                                            'TStart', tstart, ...
	                                            'TEnd',   tend);
	assert(count > 0, ['No LoCal file found: "' fsrch '".']);
	
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
	attitude = mms_fdoa_read_defatt(fname_att, tstart, tend);

	% Create GEI data
	[tt2000, ~, b_gei] = mms_fg_create_l2(fname_fg, fname_hi, fname_lo, tstart, tend, ...
	                                      'Attitude', attitude);
	
%------------------------------------%
% Write the Data                     %
%------------------------------------%

	[~, ~, ~, ~, fstart, version] = mms_dissect_filename(fname_fg);

	% Create a file name
	fname = mms_construct_filename(sc, instr, mode, 'l2', ...
	                               'OptDesc',   'gei', ...
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
		fprintf( fid, '%s   %13.6f   %13.6f   %13.6f\n', utc{ii}, b_gei(:,ii) );
	end
	
	% Close the file
	fclose(fid);
	
	if nargout() == 0
		clear fname
	end
end