%
% Name
%   mms_req_scpos
%
% Purpose
%   Request for Project SMART students.
%
%   Create AFG or DFG data in GEI so that it can be compared with RBSP mag data
%   to determine spin-axis offsets. Reads in L1A data and turns it into a L2
%   product.
%
%   This will work only on a server to which the MMS SDC is mounted to the /nfs/
%   directory.
%
% Calling Sequence
%   POS = mms_req_scpos(SC, TSTART, TEND)
%     Return the position in ECI coordinates of spacecraft SC during the time
%     interval [TSTART, TEND]. TSTART and TEND must be given in ISO-8601 format:
%     yyyy-mm-ddTHH:MM:SS.
%
%   POS = mms_req_scpos(..., TF_MEAN)
%     Indicate you want the mean position within the given time interval.
%
%   [POS, T] = mms_req_scpos(__)
%     Also return the time stamps of the position values.
%
% Parameters
%   SC              in, required, type = string
%   TSTART          in, required, type = string
%
% Returns
%   POS             out, required, type=3xN double
%   T               out, optional, type=1xN int64
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-05-16      Written by Matthew Argall
%
function [pos, t] = mms_req_scpos(sc, tstart)
	
	% Constants
	ephdir = fullfile('/nfs', 'ancillary', sc, 'defeph');

	% Make the end time 10 minutes later than the start time.
	tstart_tt2000 = spdfparsett2000([tstart, '.000000000']);
	tstart_vec    = spdfbreakdowntt2000(tstart_tt2000);
	tend_vec      = tstart_vec;
	tend_vec(5)   = tend_vec(5) + 10;
	tend_tt2000   = spdfcomputett2000(tend_vec);
	tend          = spdfencodett2000(tend_tt2000);
	tend          = tend{1}(1:19);

%------------------------------------%
% Find Files                         %
%------------------------------------%

	% Ephemeris
	ftest = fullfile(ephdir, [upper(sc) '_DEFEPH_%Y%D_%Y%D.V*']);
	[fname_eph, count] = MrFile_Search(ftest, ...
	                                   'Closest',      true, ...
	                                   'TStart',       tstart, ...
	                                   'TEnd',         tend, ...
	                                   'TimeOrder',    '%Y%D', ...
	                                   'VersionRegex', 'V[0-9]+');
	assert(count > 0, ['No ephemeris file found: "' ftest '".']);
	
%------------------------------------%
% Read Data                          %
%------------------------------------%

	% Get ephemeris data
	ephemeris = mms_fdoa_read_defeph(fname_eph, tstart, tend);
	
	% Extract position
	pos = ephemeris.Position(:, 1);
	t   = ephemeris.tt2000(1);
	
	% Units of RE
	Re  = 6371.0;
	pos = pos / Re;
end