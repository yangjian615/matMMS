%
% Name
%   mms_fsm_bkgd_fit_area
%
% Purpose
%   Find the average occurrence within a distribution of signal values. Consider
%   this to be the noise floor of the instrument.
%
% Calling Sequence
%   H_FLOOR = mms_fsm_bkgd_compute(H, BINS)
%     Compute the average occurrence H_FLOOR of each distribution in the 3D
%     histogram H with histogram bins BINS.
%
% Parameters
%   H               in, required, type = integer
%   BINS            in, required, type = double
%
% Returns
%   H_FLOOR         out, required, type=single
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2016-05-31      Written by Matthew Argall
%   2016-09-14      Renamed from mms_fsm_bkgd_noisefloor to mms_fsm_bkgd_fit_area - MRA
%
function h_floor = mms_fsm_bkgd_fit_area(h, bins)
	% Size of the histogram
	h   = double( squeeze(h) );
	szh = size(h);

	% Compute area:
	%   - Of each bin (Occurrence * bin value)
	%   - Total area at each frequency
	%   - Relative area at each frequency
	dBin     = median( diff( bins ) );
	area     = h .* dBin;
	area_rel = cumsum(area, 2) ./ repmat( sum(area, 2), 1, szh(2), 1 );
	
	% Find the location where the relative accumulated area reaches 50%.
	[~, ipct] = min( abs( area_rel - 0.5 ), [], 2 );

	% Get rid of some data
	clear area area_rel
	
	% Reconstruct indices
	%   - IPCT are indices into the columns
	%   - We want indices into the entire array.
	%   - nRows*nColumns*(iDepth-1) + nRows*(iColumn-1) + iRow
	%   - abs( area_rel(idx) - 0.5 ) should approx equal pct [the ~ value above]
	irow   = repmat([1:3]', 1, 1, szh(3));
	idepth = repmat( reshape( 1:szh(3), 1, 1, szh(3) ), 3, 1, 1 );
	idx_lo = szh(1)*szh(2)*(idepth-1) + szh(1) * max( (ipct-2), 0 )    + irow;
	idx    = szh(1)*szh(2)*(idepth-1) + szh(1) *      (ipct-1)         + irow;
	idx_hi = szh(1)*szh(2)*(idepth-1) + szh(1) * min( ipct, szh(2)-1 ) + irow;

	% Clear data
	clear irow idepth

	% Take the weighted average of +/- 1 bins from the center
	%   - Weight by occurrence
	%   - Average amplitude
	h_floor = ( h( idx_lo ) .* bins( max(ipct-1, 1) ) + ...
	            h( idx    ) .* bins( ipct   )         + ...
	            h( idx_hi ) .* bins( min(ipct+1, szh(2)) ) ) ./ ...
	          ( h( idx_lo ) + h( idx ) + h( idx_hi ) );

	% METHOD 1 -- OBSOLETE
	if false
		% Choose the largest value
		[~,inds] = max(h, [], 3);
		[~,idx]  = ind2sub( szh, inds );
		h_floor  = bins(inds);

		% Squeeze to remove a shallow dimension
		h_floor = squeeze( h_floor );
	end
end