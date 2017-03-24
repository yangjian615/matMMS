%
% Name
%   mms_fsm_fgm_compensate
%
% Purpose
%   Apply David Fischer's filter models to AFG/DFG data. Timing delays that are removed
%   in the L1A process are added back into the data before the models are applied.
%   Aliasing effects are not removed here, so the data must be low-pass filtered later.
%
% Calling Sequence
%   [TIME, DATA] = mms_fsm_bkgd_compute(MODEL, TIME, DATA, FS)
%     Apply the filter model MODEL to fluxgate magnetometer TIME and DATA
%     sampled at frequency FS.
%
%   [TIME, DATA] = mms_fsm_bkgd_compute(MODEL, TIME, DATA, DELAY, FS)
%      Additionally, restore a time delay DELAY that was removed during the
%      initial create of TIME and DATA (i.e. from L0 to L1A). DELAY should
%      be in seconds.
%
%   [TIME, DATA] = mms_fsm_bkgd_compute(MODEL, TIME, DATA, DELAY, FS, FSNEW)
%     Data is upsampled from FS to a new sampling frequency FSNEW.
%
% Parameters
%   TIME            in, required, type = in64 (cdf_time_tt2000)
%   DATA            in, required, type = Nx3 single
%   DELAY           in, optional, type = double
%   FS              in, required, type = int16
%   FSNEW           in, optional, type = int16
%
% Returns
%   TIME            out, required, type = in64
%   DATA            out, required, type = Nx3 single
%
% MATLAB release(s) MATLAB 9.0.0.341360 (R2016a)
% Required Products None
%
% History:
%   2017-01-31      Written by Matthew Argall
%
function [time, data] = mms_fsm_fgm_compensate(model, time, data, delay, fs, fsnew)

	if nargin == 4
		fs    = delay;
		delay = 0;
	elseif nargin == 5
		fsnew = fs;
	end
	
	% Compensation works up to 1024 Hz.
	%   - Beyond that, we need a second upsampling step
	if fsnew > 1024
		fstemp = 1024;
		fsout  = fsnew;
	else
		fstemp = fsnew;
		fsout  = 0;
	end

%------------------------------------%
% Restore Time Delay                 %
%------------------------------------%
	if delay ~= 0
		time = time + round( delay*1e9 );
	end

%------------------------------------%
% Frequency Compensation             %
%------------------------------------%

	% Upsample to the target product sampling frequency (1kHz) without antialising filter
	mrfprintf('logtext', ['    Upsampling from ',num2str(fs),' Hz to ',num2str(fsnew),' Hz without aliasing filter'])
	[time, data] = mms_resample( time, data, fs, fstemp, 0, 'upsample');

	% Alternative: with antialiasing...
	% disp(['    Upsampling ',instrument,' from ',num2str(xfgfs),' Hz to ',num2str(fsnew),' Hz with aliasing filter'])
	% [time, data] = mms_resample( time, data, fs, fsnew, xfgfs/2);
	% Alternatingly we can also apply the antialiasing filter later on without upsampling with this line
	% [time, data] = mms_resample( time, data, fsnew, fsnew, xfgfs/2);

	% Apply instrument models
	mrfprintf( 'logtext', ['    Compensating @', num2str(fsnew)] );
	data = compensate_xfg( data, model );
	
	% Upsample beyond 1024 S/s
	if fsout ~= 0
		[time, data] = mms_resample( time, data, fstemp, fsout);
	end

	% Remove all leftover data, as some was removed by filtering.
	[time, data] = cut_even_length( time, data );
end