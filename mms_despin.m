%
% Name
%   mms_dss_despin
%
% Purpose
%   Transform a vector field from a spinning frame to a despun frame.
%
% Calling Sequence
%   [DESPUN] = mms_despin(DATA, TIME, SC, TSTART, TEND, HK_DIR)
%     Despin a vector field DATA with CDF epoch time tags TIME. Data is
%     from MMS spacecraft SC. Look for files that start on or after TSTART
%     and before TEND. Look in HK_DIR for the housekeeping 101 files that
%     contain the sun pulse times. Return the despun data, DESPUN.
%
%   [DESPUN, T_SUN] = mms_despin(__)
%     Also return the sun pulse times. These can be used with
%     mms_dss_despin.m without having to read the entire housekeeping file
%     again. Note that T_SUN is not altered from how it appears in the
%     file. See the 'Unique' and 'Smooth' options in mms_dss_despin.
%
% Parameters
%   DATA:           in, required, type=3xN float
%   TIME            in, required, type=1xN int64/tt2000
%   SC              in, required, type=char
%   TSTART          in, optional, type=char
%   TEND            in, optional, type=char
%   HK_DIR          in, optional, type=char
%
% Returns
%   DESPUN          out, required, type=3xN float
%
% See Also
%   mms_hk_read_sunpulse.m
%   mms_file_search.m
%   mms_dss_despin.m
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-06      Written by Matthew Argall
%
function [data_despun, t_sun] = mms_despin(data, time, sc, tstart, tend, hk_dir)

	% Find the definitive attitude file name
	sunpulse_fname = mms_file_search(sc, 'fields', 'hk', 'l1b_101', ...
	                                 'Directory', hk_dir, ...
	                                 'TStart',    tstart, ...
	                                 'TEnd',      tend);

	% Make sure only one file was found.
	assert( length(sunpulse_fname) == 1, 'More than one Sun Pulse file found. Cannot proceed.' )


	% Read the sun pulse times.
	[~, t_sun] = mms_hk_read_sunpulse(sunpulse_fname{1});

	% Despin the data
	data_despun = mms_dss_despin(t_sun, time, data, ...
	                             'Unique', true, ...
	                             'Smooth', true);
end