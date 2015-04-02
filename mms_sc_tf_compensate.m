%
% NAME
%   mms_scm_tf_compensate
%
% PURPOSE
%   Routine to compute the compensation array for a transfer function.
%
% Calling Sequence
%   COMP = mms_sc_tf_compensate(TRANSFR_FN, F, N, DF)
%     This computes the complex array COMP by interpolation of the
%     transfer function, TF, at given frequencies, F, and puts the complex
%     conjugate values in the upper half of the compensation array
%
%   NOTE: N is the number of frequencies to compensate and should be EVEN !
%
%   NOTE: The NaN value for frequencies outside the range of tr_freq should be
%
%                       comp(isnan(comp)) = 0
%
%                             OR
%                       
%                       comp(isnan(comp)) = Inf 
%
%   Depending on whether the spectra is multiplied- or divided-by the transfer
%   function.
%
%   Here, it is assumed that the spectra is multiplied by the transfer function
%  in 'apply_tranf.m' and so comp(isnan(comp)) = 0.
%
% Parameters
%   TRANSFR_FN      in, required, type=3xN double
%   F               in, required, type=3xN double
%   N               in, required, type=1xN double
%   DF              in, required, type=1xN double
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-03-22      Written by Matthew Argall
%
function [comp] = mms_sc_tf_compensate(transfr_fn, f, N, df)

	% Pivot at the Nyquist frequency. The resulting array will
	% have real frequency components from 1 to N/2 and complex
	% components N/2 + 1 to end.
	pivot = N/2;
	
	% Frequencies at which to interpolate
	freq_out = df * (1:pivot);
	
	% Create the compensation array
	%   - Frequency range extends from f0 to fN
	comp           = zeros(3, length(freq_out));
	comp(1, :)     = interp1(f(1, :), transfr_fn(1, :), freq_out, 'linear');
	comp(2, :)     = interp1(f(2, :), transfr_fn(2, :), freq_out, 'linear');
	comp(3, :)     = interp1(f(3, :), transfr_fn(3, :), freq_out, 'linear');
	comp(:, pivot) = abs( comp(:, pivot));
	
	% Complete teh compensation array
	%  - Append the DC component
	%  - Reflect the complex component beyond N/2
	%  - do not duplicate the Nyquist
	%  - Real and imaginary frequencies have the same spectral components
	comp(isnan(comp)) = 0;
	comp              = [ ones(3,1) comp conj( fliplr( comp(:, 1:end-1) ) )];
end