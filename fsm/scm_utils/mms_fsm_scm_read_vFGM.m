%
% Name
%   mms_fsm_scm_read_vFGM
%
% Purpose
%   Read SCM data and sort it into fast and slow survey portions.
%
% Calling Sequence
%   DATA = mms_fsm_scm_read(FILES, FGM)
%     Read SCM data files with names FILES over the course of an orbit, using FGM L1A
%     and L2Pre data obtained from mms_fsm_bkgd_fgm_read.m. Return data in the
%     structure DATA.
%
%   DATA = mms_fsm_scm_read(__, TRANGE)
%     If survey data is being read, the TRANGE specifies the time range of a single
%     orbit, starting at the beginning of slow survey and ending at the end of the
%     following fast survey period. Times should be formatted as 'yyyy-mm-ddTHH:MM:SS'
%     or as empty strings.
%
%   DATA = mms_fsm_scm_read(__, FC)
%     Give the bass-cut frequency at which to high-pass filter the data.
%
% Parameters
%   FILES           in, required, type = char/cell
%   FGM             in, required, type = structure
%                   Required fields are::
%                       't'     - TT2000 epoch times
%                       'flag'  - Flag indicating operational modes
%                                    0 = Slow    1 = Fast
%                                    0 = Lo-     2 = Hi-Range
%                                    0 = Before  4 = After perigee (for slow survey)
%                                    0 = Deck32  8 = Deck64
%   TRANGE          in, optional, type = 1x2 cell
%   FC              in, optional, type = double, default = []
%
% Returns
%   DATA            out, required, type=struct
%                   Fields are:
%                       't'      - TT2000 epoch times
%                       'b'      - Magnetic field in 123 coordinates
%                       'flag'   - Flag indicating operational mode:
%                                    0 = Slow    1 = Fast
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-06-01      Written by Matthew Argall
%   2015-10-01      Renamed from mms_fsm_bkgd_scm_read to mms_fsm_scm_read. FC now
%                     defaults to the empty array. - MRA
%
function scm = mms_fsm_scm_read_vFGM(file, fgm, trange, fc)

	if nargin() < 4
		fc = [];
	end

%------------------------------------%
% Read Data                          %
%------------------------------------%
	
	% Read the L2Pre files
	[b_123, t] = mms_fsm_scm_read_l1b(file, trange);
	
	% Categorize the data
	scm = mms_fsm_scm_read_categorize_vFGM(t, b_123, fgm);
	
	% Clear the data
	clear t b_123

%------------------------------------%
% Interpolate Over NaNs              %
%------------------------------------%

	% Find NaNs
	tf_nan = isnan(scm.b(1,:));
	nNaN   = sum(tf_nan);
	if nNaN > 0
		% Issue warning
		mrfprintf('logwarn', '%d NaN values found. Inteprolating.', nNaN);
		
		% Convert time to doubles (required by interp1)
		t_ssm = MrCDF_epoch2ssm(scm.t);

		% Interpolate
		scm.b(1,tf_nan) = interp1( t_ssm(~tf_nan), scm.b(1,~tf_nan), t_ssm(tf_nan) );
		scm.b(2,tf_nan) = interp1( t_ssm(~tf_nan), scm.b(2,~tf_nan), t_ssm(tf_nan) );
		scm.b(3,tf_nan) = interp1( t_ssm(~tf_nan), scm.b(3,~tf_nan), t_ssm(tf_nan) );
		
		% Clear data
		clear t_ssm
	end

%------------------------------------%
% Isolate Modes                      %
%------------------------------------%
	
	% Flag
	%   - Bit 1 set:     Fast and Brst
	%   - Bit 1 not set: Slow 
	%   - iFast serves as iBrst when iSlow is empty
	iFast = find( bitget(scm.flag, 1, 'uint8') );
	iSlow = find( bitget(scm.flag, 1, 'uint8') == 0 );

%------------------------------------%
% High-Pass Filter                   %
%------------------------------------%
	if ~isempty(fc)
		% Slow
		if ~isempty(iSlow)
			% Create filter
			%   - butter() + filtfilt() was producing NaNs. This was the solution
			%   - https://www.mathworks.com/matlabcentral/answers/62553-filtfilt-function-returning-nan-at-certain-frequencies
			fN = scm.sr( iSlow(1) ) / 2.0;
			d  = fdesign.highpass('N,F3dB', 9, fc/fN);
			h  = design(d, 'butter');
		
			% Apply the filter
			scm.b(1,iSlow) = filtfilt( h.sosMatrix, h.ScaleValues, double( scm.b(1,iSlow) ));
			scm.b(2,iSlow) = filtfilt( h.sosMatrix, h.ScaleValues, double( scm.b(2,iSlow) ));
			scm.b(3,iSlow) = filtfilt( h.sosMatrix, h.ScaleValues, double( scm.b(3,iSlow) ));
		end

		% Fast & Brst
		if ~isempty(iFast)
			% Create filter
			%   - butter() + filtfilt() was producing NaNs. This was the solution
			%   - https://www.mathworks.com/matlabcentral/answers/62553-filtfilt-function-returning-nan-at-certain-frequencies
			fN = scm.sr( iFast(1) ) / 2.0;
			d  = fdesign.highpass('N,F3dB', 9, fc/fN);
			h  = design(d, 'butter');
		
			% Apply the filter
			scm.b(1,iFast) = filtfilt( h.sosMatrix, h.ScaleValues, double( scm.b(1,iFast) ));
			scm.b(2,iFast) = filtfilt( h.sosMatrix, h.ScaleValues, double( scm.b(2,iFast) ));
			scm.b(3,iFast) = filtfilt( h.sosMatrix, h.ScaleValues, double( scm.b(3,iFast) ));
		end
	end
end