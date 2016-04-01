%
% Name
%   mms_req_edpgse
%
% Purpose
%   Request for Charlie Farrugia
%
%      GSE Data on 2015-08-15 from 13:20:00 to 13:30:00
%
%   Transform EDP DCE data from DSL to GSE during the interval in which
%   the extremely large FTE is observed.
%
% Calling Sequence
%   mms_req_edpgse(SC, INSTR, MODE, TSTART, TEND)
%     Read in L1A data from spacecraft SC, instrument INSTR, and telemetry mode
%     MODE during the time interval [TSTART, TEND], turn the data into a L2
%     product in GEI coordinates, then write the time and field to an ASCII file.
%
%   FNAME = mms_req_edpgse(__)
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
function fname = mms_req_edpgse(sc, instr, mode, tstart, tend, outdir)

	% Default inputs
	sc      = 'mms1';
	tstart  = '2015-10-16T13:00:00';
	tend    = '2015-10-16T13:15:00';

	% Output directory
	if nargin() < 6
		outdir = '/home/argall/data/';
	end
	
	% Constants
	attdir = fullfile('/nfs', 'ancillary', sc, 'defatt');
	hk_dir = fullfile('/nfs', 'hk');
	
%------------------------------------%
% Find Files                         %
%------------------------------------%
	% EDP Fast
	[fedp_fast, count, fsrch] = mms_file_search(sc, 'edp', 'fast', 'ql',     ...
	                                            'OptDesc',   'dce2d',          ...
	                                            'TimeOrder', '%Y%M%d%H%m%S', ...
	                                            'TStart',    tstart,         ...
	                                            'TEnd',      tend);
	assert(count > 0, ['No EDP Fast file found: "' fsrch '".']);
	
	% EDP Slow
	[fedp_slow, count, fsrch] = mms_file_search(sc, 'edp', 'fast', 'ql',     ...
	                                            'OptDesc',   'dce2d',          ...
	                                            'TimeOrder', '%Y%M%d%H%m%S', ...
	                                            'TStart',    tstart,         ...
	                                            'TEnd',      tend);
	assert(count > 0, ['No EDP Fast file found: "' fsrch '".']);
	
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
	% Attitude 
	[attitude, att_hdr] = mms_fdoa_read_defatt(fname_att, tstart, tend);
	zMPA                = att_hdr.zMPA(:,1);

	% EDP Fast
	edp_ql_fast = mms_edp_read_ql(fedp_fast, tstart, tend);
	edp_ql_slow = mms_edp_read_ql(fedp_slow, tstart, tend);
	
	% Combine the data
	t_edp          = [edp_ql_slow.tt2000, edp_ql_fast.tt2000];
	[t_edp, isort] = sort(t_edp);
	e_dsl          = [edp_ql_slow.E_dsl, edp_ql_fast.E_dsl];
	e_dsl          = e_dsl(:,isort);
keyboard
%------------------------------------%
% Rotate to GSE                      %
%------------------------------------%
	% Transformation matrix
	gei2dsl = mms_fdoa_xgei2despun(attitude, t_edp, 'L');
	dsl2gei = permute(gei2dsl, [2,1,3]);
	
	% Transform to GEI
	e_gei = mrvector_rotate(dsl2gei, e_dsl);
	
	% Transform matrix GEI -> GSE
	%   - Modified Julian Date (mjd).
	%   - UTC seconds since midnight (ssm).
	timevec = MrCDF_Epoch_Breakdown(t_edp)';
	mjd     = date2mjd(timevec(1,:), timevec(2,:), timevec(3,:));
	ssm     = MrCDF_epoch2ssm(t_edp);
	GEI2GSE = gei2gse(mjd, ssm);

	% Transform to GSE
	e_gse = mrvector_rotate(GEI2GSE, e_gei);

%------------------------------------%
% Write the Data                     %
%------------------------------------%
	[sc, instr, mode, level, tstart, version, optdesc] = mms_dissect_filename(fedp_fast);

	% Create a file name
	fname = mms_construct_filename(sc, 'edp', 'srvy', 'l2',       ...
	                               'OptDesc',   'dce2d-unh-gse',    ...
	                               'TStart',    '20151016130000', ...
	                               'Directory', outdir,           ...
	                               'Version',   version);
	fname(end-2:end) = 'dat';

	% Convert time to UTC
	utc = spdfencodett2000(t_edp);

	% Opend the file
	fid = fopen(fname, 'w');
	
	% Write to the file
	fprintf(fid, '                          %s              %s              %s              %s\n', 'UTC', 'Ex', 'Ey', 'Ez');
	for ii = 1 : length(t_edp)
		fprintf( fid, '%s   %13.6f   %13.6f   %13.6f\n', utc{ii}, e_gse(:,ii) );
	end
	
	% Close the file
	fclose(fid);
	
	if nargout() == 0
		fprintf('File written to: %s.\n', fname)
		clear fname
	end
end