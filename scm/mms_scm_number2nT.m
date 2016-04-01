%
% Name
%   mms_scm_number2nT
%
% Purpose
%   Convert a 32-bit integer to a 32-bit double. This converts magnetometer
%   search coil values from numbers to nano-tesla.
%
% Calling Sequence
%   B_NT = mms_scm_number2nT(B_SC);
%     Convert searchcoil magnetometer data from a number B_SC to
%     nano-tesla, B_nT.
%
% Parameters
%   B_SC            in, required, type=3xN int32
%
% Returns
%   B_NT            out, required, type=3xN single
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-03-24      Written by Matthew Argall
%   2015-05-11      Mark
%
function [B_nT] = mms_scm_number2nT(B_sc)
	
	%
	% Translate from a number to nano-tesla (32-bit int to 32-bit float)
	%   - The first  2^15 - 1 numbers are negative
	%   - The number 0
	%   - The second 2^15 - 1 numbers are positive
	%   - Normalize by 2^16 - 1 (65535.0)
	%
	% From Olivier le Contel:
	%
	%   Mark's  L1A CDF files are already signed between 32767 and - 32767.
	%   Furthermore, SCM output voltage is 12.25 V not 10V, so the conversion
	%   factor that we use to go from TM units to Volts, is only 
	%
	%       5.0/0.4030/2.0^16
	%
	B_nT = single( (5.0 / 0.4030) .* double(B_sc) ./ 65535.0);
end