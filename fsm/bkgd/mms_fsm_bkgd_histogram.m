%
% Name
%   mms_fsm_bkgd_histogram
%
% Purpose
%   Histogram amplitude, phase, and power spectrograms in time to find the
%   occurrence frequency of each quantity.
%
% Calling Sequence
%   DATA = mms_fsm_bkgd_compute(DATA)
%     Using the amplitude, phase, power spectra given as a function of time
%     and frequency in the data structure DATA, histogram each quantity in time
%     to find the occurrence as a function of frequency. Append results to the
%     input structure.
%
%   DATA = mms_fsm_bkgd_compute(DATA, BINSIZE)
%     Specify a bin size for the histograms.
%
% Parameters
%   DATA            in, required, type = struct
%                   Fields are:
%                       't'     - Center times of the time series
%                       'f'     - Frequencies
%                       'amp'   - Wave amplitude
%                       'phase' - Phase
%                       'psd'   - Power spectral density
%   BINSIZE         in, required, type = integer, default=0.1
%
% Returns
%   DATA            in, required, type = struct
%                   New Fields are:
%                       'amp_hist'   - Histogrammed amplitude
%                       'phase_hist' - Histogrammed phase
%                       'psd_hist'   - Histogrammed psd
%                       'amp_bins'   - Histogram bins for amplitude
%                       'phase_bins' - Histogram bins for phase
%                       'psd_bins'   - Histogram bins for psd
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2016-05-31      Written by Matthew Argall
%
function data = mms_fsm_bkgd_histogram(data, binsize)

	if nargin < 2
		binsize = 0.1;
	end

%------------------------------------%
% Histogram                          %
%------------------------------------%
	% Name and number of fields
	names   = fieldnames(data);
	nFields = length(names);
	nTimes  = length(data.t);
	nFreqs  = length(data.f);
	
	% Loop over each data quantity
	for ii = 1 : nFields

	%------------------------------------%
	% Select Data                        %
	%------------------------------------%
		% Amplitude
		if strcmp( names{ii}, 'amp' )
			temp    = log10( data.(names{ii}) );
			range   = [-4.0, 4.0];
		
		% Phase
		elseif strcmp( names{ii}, 'phase' )
			temp    = data.(names{ii});
			range   = [-180.0, 180.0];
		
		% PSD
		elseif strcmp( names{ii}, 'psd' )
			temp    = log10( data.(names{ii}) );
			range   = [-10, 5];
		
		% GAIN
		elseif strcmp( names{ii}, 'gain' )
			temp    = log10( data.(names{ii}) );
			range   = [-3, 3];
		
		% Phase Offset
		elseif strcmp( names{ii}, 'phase_shift' )
			temp    = data.(names{ii});
			range   = [-90,90];
		
		% PSD Ratio
		elseif strcmp( names{ii}, 'psd_rat' )
			temp    = log10( data.(names{ii}) );
			range   = [-5, 5];
		
		% Other
		else
			continue
		end
		
		% Bin sizes
		nBins   = 1 + floor( (range(2) - range(1)) / binsize );
		bins    = linspace(range(1), range(2), nBins);

	%------------------------------------%
	% Histogram                          %
	%------------------------------------%
		% Allocate memory
		h = zeros(3, nBins, nFreqs);

		% Loop over all frequencies
		for jj = 1 : nFreqs
			h(1,:,jj) = histc( temp(1,:,jj), bins );
			h(2,:,jj) = histc( temp(2,:,jj), bins );
			h(3,:,jj) = histc( temp(3,:,jj), bins );
		end

	%------------------------------------%
	% Save                               %
	%------------------------------------%
		data.( [names{ii} '_hist'] )  = uint32( reshape( h, 3, 1, nBins, nFreqs ) );
		data.( [names{ii} '_bins'] )  = single( bins );
	end
end