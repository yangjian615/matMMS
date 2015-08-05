%
% Name
%   mms_req_fgdespinl1b
%
% Purpose
%   Request for Li-Jen Chen
%
%   Despin mms3 dfg l1b burst data on 2015-06-22 and 2015-06-23. Use the sunpulse
%   times to despin so that x-DMPA passes through the sun sensor.
%
%   This will work only on a server to which the MMS SDC is mounted to the /nfs/
%   directory.
%
% Calling Sequence
%   mms_req_fgdespinl1b(SC, INSTR, MODE, TSTART, TEND)
%     Read in L1A data from spacecraft SC, instrument INSTR, and telemetry mode
%     MODE during the time interval [TSTART, TEND], turn the data into a L2
%     product in GEI coordinates, then write the time and field to an ASCII file.
%
%   FNAME = mms_req_fgdespinl1b(__)
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
%   2015-05-17      Written by Matthew Argall
%
function fname = mms_req_fgdespinl1b(sc, instr, mode, tstart, tend, outdir)

	% Default inputs
	sc     = 'mms3';
	instr  = 'dfg';
	mode   = 'brst';
	level  = 'l1b';
	tstart = '2015-06-22T00:00:00';
	tend   = '2015-06-22T24:00:00';
	
	% Output directory
	if nargin() < 6
		outdir = pwd();
	end
	
	% Constants
	dssdir = fullfile('/nfs', 'hk');
	attdir = fullfile('/nfs', 'ancillary', sc, 'defatt');
	
%------------------------------------%
% Find Files                         %
%------------------------------------%
	% DATA
	[fname_fg, count, fsrch] = mms_file_search(sc, instr, mode, level, ...
	                                           'DaySubdir', true, ...
	                                           'TimeOrder', '%Y%M%d%H%m%S', ...
	                                           'TStart',    tstart, ...
	                                           'TEnd',      tend);
	assert(count > 0, ['No DFG file found: "' fsrch '".']);
	
	% Sunpulse
	[fname_dss, count, fsrch] = mms_file_search(sc, 'fields', 'hk', 'l1b', ...
	                                            'OptDesc', '101', ...
	                                            'SDCroot', dssdir, ...
	                                            'TStart',  tstart, ...
	                                            'TEnd',    tend);
	assert(count > 0, ['No DSS file found: "' fsrch '".']);
	
%------------------------------------%
% Create L2                          %
%------------------------------------%

	% Get attitude data
	sunpulse = mms_dss_read_sunpulse(fname_dss, tstart, tend, 'UniquePulse', true);

	% Read l1b burst data
	fg_l1b = mms_fg_read_l1b(fname_fg, tstart, tend);

	% Rotate OMB to SMPA
	omb2smpa    = mms_fg_xomb2smpa();
	b_smpa      = omb2smpa * fg_l1b.b_omb(1:3, :);
	smpa2despun = mms_dss_xdespin( sunpulse, fg_l1b.tt2000 );
	b_dmpa      = mrvector_rotate( smpa2despun, b_smpa );

	% Create strings
	t_utc = spdfencodett2000( fg_l1b.tt2000 );
	
	% Clear data that will no longer be used.
	clear fg_l1b b_smpa omb2smpa smpa2despun
%------------------------------------%
% Write the Data                     %
%------------------------------------%

	% Dissect the first file name to get its start time
	[~, ~, ~, ~, fstart, version] = mms_dissect_filename(fname_fg{1});
	fstart = fstart(1:8);

	% Create a file name
	fname = mms_construct_filename(sc, instr, mode, 'ql', ...
	                               'TStart',    fstart, ...
	                               'Directory', outdir, ...
	                               'Version',   version);
	fname(end-2:end) = 'dat';

	% Opend the file
	fid = fopen(fname, 'w');

	% Write to the file
	fprintf(fid, '                           %s              %s              %s              %s\n', 'UTC', 'Bx', 'By', 'Bz');
	for ii = 1 : length(t_utc)
		fprintf( fid, '%30s   %13.6f   %13.6f   %13.6f\n', t_utc{ii}, b_dmpa(:,ii) );
	end
	
	% Close the file
	fclose(fid);
	
	if nargout() == 0
		disp(['File saved to: ' fname]);
		clear fname
	end
end