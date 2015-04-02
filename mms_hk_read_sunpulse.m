%
% Name
%   mms_hk_read_sunpulse
%
% Purpose
%   Read houskeeping 0x101 sun pulse and related data.
%
% Calling Sequence
%   [HK_TT2000, SUN_PULSE, FLAG, PERIOD] = mms_hk_sunpulse(FILENAME)
%     Read housekeeping 0x101 packet times HK_TT2000, time of last sun
%     pulse in tt2000, SUN_PULSE, sun pulse flag, FLAG, and external sun
%     pulse period in microseconds, PERIOD.
%
% Parameters
%   FILENAME        in, required, type=char
%
% Returns
%   HK_TT2000       out, required, type=Nx1 int64/tt2000
%   SUN_PULSE       out, required, type=Nx1 int64/tt2000
%   FLAG            out, optional, type=Nx1
%                   Possible values include:
%                     0: s/c sun pulse
%                     1: s/c pseudo sun pulse
%                     2: s/c CIDP generated speudo sun pulse
%   PERIOD          out, optional, type=Nx1
%                   Only returned when FLAG=0 and only on the second and
%                     subsequent received sun pulses from the s/c.
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-03-25      Written by Matthew Argall
%
function [hk_tt2000, sun_pulse, flag, period] = mms_hk_read_sunpulse(filename)

	% Ensure the file exists
	assert( exist(filename, 'file') == 2, ['HK file does not exist ' filename '".']);
	
	% Get the spacecraft identifier
	[~, name] = fileparts(filename);
	sc        = name(1:4);
	assert( ~isempty(regexp(sc, 'mms[1-4]', 'once')), ['Unknown observatory "' sc '".']);
	
	
	epoch_name     = 'Epoch';
	sunpulse_name  = [sc '_101_sunpulse'];
	flag_name      = [sc '_101_sunssps'];
	period_name    = [sc '_101_iifsunper'];
	
	% Read the data
	hk_data = spdfcdfread(filename, ...
		                    'Variables',    {epoch_name sunpulse_name flag_name period_name}, ...
											  'KeepEpochAsIs', true, ...
												'CombineRecords', true);
	% Extract the data
	hk_tt2000 = hk_data{1};
	sun_pulse = hk_data{2};
	flag      = hk_data{3};
	period    = hk_data{4};
end