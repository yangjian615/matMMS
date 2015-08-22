%
% Name
%   mms_fdoa_test_interp_ephem
%
% Purpose
%   Test interpolation of FDOA ephemeris position and velocity.
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-08-21      Written by Matthew Argall
%

get_data = true;

if get_data
%------------------------------------%
% Inputs                             %
%------------------------------------%
	sc         = 'mms2';
	tstart     = '2015-07-22T09:00:00';
	tend       = '2015-07-22T17:00:00';
	eph_dir    = fullfile('/nfs', 'ancillary', sc, 'defeph');
	
%------------------------------------%
% Find Files                         %
%------------------------------------%
	% DFG QL File
	[fg_fname, count, str] = mms_file_search(sc, 'dfg', 'srvy', 'ql', ...
	                                         'TStart',  tstart, ...
	                                         'TEnd',    tend);
	assert(count > 0, ['DFG file not found: "' str '".']);
	
	% Ephemeris files
	eph_ftest = fullfile( eph_dir, [upper(sc) '_DEFEPH_%Y%D_%Y%D.V*'] );
	[defeph_fname, count] = MrFile_Search(eph_ftest, ...
	                                      'Closest',      true, ...
	                                      'TimeOrder',    '%Y%D', ...
	                                      'TStart',       tstart, ...
	                                      'TEnd',         tend, ...
	                                      'VersionRegex', 'V[0-9]{2}');
	assert(count > 0, ['No definitive ephemeris file found: "' eph_ftest '".']);

%------------------------------------%
% Read & Interpolate Data            %
%------------------------------------%
	% Get data
	fg_ql    = mms_fg_read_ql(fg_fname,           tstart, tend);
	ephem    = mms_fdoa_read_defeph(defeph_fname, tstart, tend);
	
	% Interpolate ephemeris position and velocity to DFG times
	[r_gei, v_gei] = mms_fdoa_interp_ephem(ephem.tt2000, ephem.Position, fg_ql.tt2000, ephem.Velocity);
end

%------------------------------------%
% Plot                               %
%------------------------------------%
	
% Reference time
t0 = min( [fg_ql.tt2000(1), ephem.tt2000(1)] );
t_fg  = MrCDF_epoch2sse( fg_ql.tt2000, t0 );
t_eph = MrCDF_epoch2sse( ephem.tt2000, t0 );

%
% SMPA
%
fig_ephem = figure();

% Position
subplot(2,2,1)
plot( t_eph, ephem.Position );
title([ 'Ephemeris Position'] );
ylabel( {'r', '(km)'} );

% Interpolated Position
subplot(2,2,2)
plot( t_fg, r_gei );
title([ 'Interpolated Position'] );
ylabel( {'r', '(km)'} );

% Velocity
subplot(2,2,3)
plot( t_eph, ephem.Velocity );
title([ 'Ephemeris Velocity'] );
ylabel( {'v', '(km/s)'} );

% Interpolated Velocity
subplot(2,2,4)
plot( t_fg, v_gei );
title([ 'Interpolated Velocity'] );
ylabel( {'v', '(km/s)'} );

