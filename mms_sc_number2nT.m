%
% Name
%   mms_sc_number2nT
%
% Purpose
%   Convert a 32-bit integer to a 32-bit double. This converts magnetometer
%   search coil values from numbers to nano-tesla.
%
% Calling Sequence
%   B_NT = mms_sc_number2nT(B_SC);
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
%
function [B_nT] = mms_sc_number2nT(B_sc)
	
	% Translate from a number to nano-tesla (32-bit int to 32-bit single)
	%   - The first  2^15 - 1 numbers are negative
	%   - The second 2^15     numbers are positive
	%   - Normalize by 2^16 - 1
	B_nT = single( 10.0 .* (double(B_sc) - 32767) ./ 65535.0);
end