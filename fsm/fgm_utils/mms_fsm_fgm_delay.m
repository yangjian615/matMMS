%
% Name
%   mms_fsm_fgm_delay
%
% Purpose
%   Return the time delay imposed by the electronics of the AFG and DFG instrumetns.
%
% Calling Sequence
%   DELAY = mms_fsm_bkgd_compute(SC, INSTR, MODE)
%     Return the time delay DELAY for MMS spacecraft SC, instrument INSTR, and
%     operatinoal mode MODE. Values for mode are
%       0 = ADCA   (afg)       1 = ADCB   (afg)
%       0 = DEC64  (dfg)       1 = DEC32  (dfg)
%
% Parameters
%   SC              in, required, type = char
%   INSTR           in, required, type = char
%   MODE            in, required, type = char
%
% Returns
%   DELAY           out, required, type = double
%
% MATLAB release(s) MATLAB 9.0.0.341360 (R2016a)
% Required Products None
%
% History:
%   2017-01-31      Written by Matthew Argall
%
function delay = mms_fsm_fgm_delay(sc, instr, fs, mode)
	
	%
	% Timing correction applied at the L1A level
	%   - L1A files have removed this delay (Mark Chutter)
	%   - For the merged product to work with David Fischer's
	%     models, we need to add it back in.
	%

	% AFG
	if strcmp(instr, 'afg')
		% IWF supplied delay for each AFG sample rate, each OBS in seconds
		% subtracting gives sample center time
		%   delay_afg_8   = [0.94946, 0.94907, 0.94759, 0.94916];
		switch fs
			case 128
				delay = [0.01196, 0.01157, 0.01009, 0.01166];
			case 32
				delay = [0.19946, 0.19907, 0.19759, 0.19916];
			case 16
				delay = [0.44946, 0.44907, 0.44759, 0.44916];
			case 8
				delay = [0.94946, 0.94907, 0.94759, 0.94916];
			otherwise
				error( ['Invalid sampling rate ' num2str(fs) 'Hz for ' upper(instr) '.'] );
		end

		% Select from spacecraft
		delay = delay(isc);
	
	% DFG
	elseif strcmp(instr, 'dfg')
		% IWF supplied delay for each DFG sample rate in seconds
		% [128/s, 32/s, 16/s, 8/s] for DEC 32, DEC64 is additional 7.81 ms
		% subtracting gives sample center time
		%   delay_dfg_rate = [0.00834, 0.19584, 0.44584, 0.94584];
		switch fs
			case 128
				delay = 0.00834;
			case 32
				delay = 0.19584;
			case 16
				delay = 0.44584;
			case 8
				delay = 0.94584;
		end
		
		% DEC64 is an additional 7.81 ms
		if mode == 0
			delay = delay + 0.00781;
		end
	else
		error( ['Invalid FGM instrument: "' upper(instr) '".'] )
	end
end