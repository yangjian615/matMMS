%
% Name
%   mms_fdoa_xbcs2despun
%
% Purpose
%   Create a coordinate transformation from the spinning spacecraft
%   body coordinate system (BCS) to a despun coordinate system.
%
%   Process:
%     1. Interpolate right ascension and declination.
%     2. Build transformation matrix.
%
% Calling Sequence
%   BSCS2DESPUN = mms_fdoa_xdespun2gei(ATTITUDE, TIME)
%     Create the transformation matrix using a structure of attitude
%     data returned by mms_fdoa_read_defatt.m. The right ascension
%     and declination of the angular momentum axis 'L' will be
%     interpolated to TIME.
%
%   BSCS2DESPUN = mms_fdoa_xdespun2gei(ATTITUDE, TIME, TYPE)
%     Specify the frame from which to rotate. Choices are:
%       'L' - Angular momentum frame.
%       'P' - Spacecraft Moment of Inertia frame.
%       'w' - Spacecraft spin frame.
%       'z' - Spacecraft body frome.
%
%   [..., RA, DEC] = mms_fdoa_xbcs2despun(__)
%     Also return the interpolated right ascension and declination.
%
% Parameters
%   ATTITUDE        in, required, type=struct
%   TIME            in, required, type=1xN int64 (cdf_time_tt2000)
%   TYPE            in, required, type=char
%
% Returns
%   BCS2DESPUN      out, required, type=3x3xN double
%   RA              out, optional, type=1xN double
%   DEC             out, optional, type=1xN double
%   ATTITUDE        out, optional, type=struct
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-13      Written by Matthew Argall
%   2015-08-21      Extrapolate for TIME outside range of ATTITUDE.TIME. - MRA
%
function [gei2despun, ra, dec] = mms_fdoa_xgei2despun(attitude, time, type)

	if nargin < 3
		type = 'L';
	end

%------------------------------------%
% Interpolate                        %
%------------------------------------%
	% Extract the right ascension and declination
	ra  = attitude.(type)(1,:);
	dec = attitude.(type)(2,:);

	% Unwrap the phase
	ra  = MrPhaseUnwrap(ra,  360.0);
	dec = MrPhaseUnwrap(dec, 360.0);
	
	% Must interpolate with doubles.
	t0       = attitude.tt2000(1);
	att_sse  = double(attitude.tt2000 - t0) * 1e-9;
	time_sse = double(time            - t0) * 1e-9;
	
	% Indicate extrapolating
	if time_sse(1) < att_sse(1)
		nExtrap = find( time_sse >= att_sse(1), 1, 'first' ) - 1;
		mrfprintf( 'logwarn', 'Extrapolating %d points before.', nExtrap );
	end
	if time_sse(end) > att_sse(end)
		nExtrap = length(time_sse) - find( time_sse <= att_sse(end), 1, 'last' );
		mrfprintf( 'logwarn', 'Extrapolating %d points after.', nExtrap );
	end
	
	% Interpolate
	%   - Extrapolate points outside of attitude time range.
	ra  = interp1(att_sse, ra,  time_sse, 'linear', 'extrap');
	dec = interp1(att_sse, dec, time_sse, 'linear', 'extrap');

%------------------------------------%
% Create Transformation Matrix       %
%------------------------------------%
	% Dissect the time
	timevec = MrCDF_Epoch_Breakdown(time)';

	% Seconds since midnight
	ssm = timevec(4,:) * 3600.0 + ...
	      timevec(5,:) * 60.0   + ...
	      timevec(6,:)          + ...
	      timevec(7,:) * 1e-3   + ...
	      timevec(8,:) * 1e-6   + ...
	      timevec(9,:) * 1e-9;

	% Get transformation matrix
	gei2despun = gei2scs(timevec(1,:), timevec(2,:), timevec(3,:), ssm, ra, dec);
end