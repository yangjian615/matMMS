%
% Name
%   mms_fsm_bkgd_read
%
% Purpose
%   Read orbit (for survey) or burst background calibration files.
%
% Calling Sequence
%   DATA = mms_fsm_bkgd_read( FILES )
%     Read FSM background calibration data from files with names provided in
%     FILES into the data structure DATA.
%
% Parameters
%   FILES           in, required, type = char/cell
%
% Returns
%   DATA            out, required, type=struct
%                   Has the following tags:
%                     'f'          - Frequency of each distribution
%                     'comp'       - Vector component
%                     'flag'       - Operations flag
%                     'amp_bins'   - Amplitude bins of each histogram distribution
%                     'phase_bins' - Phase bins of each histogram distribution
%                     'psd_bins'   - PSD bins of each histogram distribution
%                     'amp_hist'   - Amplitude histogram with dimensions [f, bin, flag, comp]
%                     'phase_hist' - Phase histogram with dimensions [f, bin, flag, comp]
%                     'psd_hist'   - PSD histogram with dimensions [f, bin, flag, comp]
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2016-07-15      Written by Matthew Argall
%
function data = mms_fsm_bkgd_read(files)
	
	if ischar(files)
		files = { files };
	end

%------------------------------------%
% Variable Names                     %
%------------------------------------%

	% Parse file name to get variable name prefix and suffix
	[sc, instr, mode, level, ~, ~, optdesc] = mms_dissect_filename(files{1});
	mag_instr = regexp(optdesc, '-', 'split');
	prefix    = [sc '_' instr '_' ];
	suffix    = ['_' mode '_' level];

	% Create Variable Names
	amp_hist_vname   = [ prefix 'amp'   '_hist'  suffix ];
	phase_hist_vname = [ prefix 'phase' '_hist'  suffix ];
	psd_hist_vname   = [ prefix 'psd'   '_hist'  suffix ];

%------------------------------------%
% Read Data                          %
%------------------------------------%
	% Create space for data
	nFiles     = length(files);
	amp_hist   = cell(1, 0);
	phase_hist = cell(1, 0);
	psd_hist   = cell(1, 0);
	hflag      = zeros(1, 0, 'uint8');

	% Loop over files
	for ii = 1 : length(files)
		% Read the data
		[ psd_temp,   f, psd_bins, flag_temp, comp ] = MrCDF_Read( files{ii}, psd_hist_vname );
%		[ amp_temp,   ~, amp_bins ]                  = MrCDF_Read( files{ii}, amp_hist_vname );
%		[ phase_temp, ~, phase_bins ]                = MrCDF_Read( files{ii}, phase_hist_vname );
		
		% Find unique flags
		%   - Identify flags that have already been read
		%   - Determine where the data associated with the flag is located
		[tf_member, iflag] = ismember(flag_temp, hflag);
		
		% Loop over each flag
		for jj = 1 : length( flag_temp )
			% Sum data for flags that already exist
			if tf_member(jj)
				idx             = iflag(jj);
%				amp_hist{idx}   = amp_hist{idx}   + squeeze(  amp_temp(:,:,jj,:));
%				phase_hist{idx} = phase_hist{idx} + squeeze(phase_temp(:,:,jj,:));
				psd_hist{idx}   = psd_hist{idx}   + squeeze(  psd_temp(:,:,jj,:));
			
			% Store data for new flags
			else
				hflag      = [ hflag      flag_temp(jj)                 ];
%				amp_hist   = [ amp_hist   squeeze(  amp_temp(:,:,jj,:)) ];
%				phase_hist = [ phase_hist squeeze(phase_temp(:,:,jj,:)) ];
				psd_hist   = [ psd_hist   squeeze(  psd_temp(:,:,jj,:)) ];
			end
		end
	end
	
	clear amp_temp phase_temp psd_temp
	
	amp_bins   = [];
	phase_bins = [];
	
%------------------------------------%
% Collect Data                       %
%------------------------------------%
	% Order data as [comp, flag, bins, freq]
	%   - Later, this is how the histogram will be fit to the noise floor.
	%   - Historical from how the orbit/brst files are processed.
	data = struct( 'f',          f,          ...
	               'comp',       comp,       ...
	               'flag',       hflag,      ...
	               'amp_bins',   single([]), ... % amp_bins,   ...
	               'phase_bins', single([]), ... % phase_bins, ...
	               'psd_bins',   psd_bins,   ...
	               'amp_hist',   uint32([]), ... % permute( cat( 4, amp_hist{:}   ), [3,4,2,1] ), ...
	               'phase_hist', uint32([]), ... % permute( cat( 4, phase_hist{:} ), [3,4,2,1] ), ...
	               'psd_hist',   permute( cat( 4, psd_hist{:}   ), [3,4,2,1] )  ...
	             );
end