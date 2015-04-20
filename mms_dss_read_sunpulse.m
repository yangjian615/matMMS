%
% Name
%   mms_hk_read_sunpulse
%
% Purpose
%   Read houskeeping 0x101 sun pulse and related data.
%
% Calling Sequence
%   SUNPULSE = mms_hk_sunpulse(FILENAME)
%     Read housekeeping 0x101 packet times, 'Epoch', time of last sun
%     pulse in tt2000, 'SunPulse', sun pulse flag, 'Flag', and external sun
%     pulse period in microseconds, 'Period'. Return a structure SUNPULSE
%     with those as fields.
%
%   [__] = mms_hk_sunpulse(__, 'ParamName', ParamValue)
%     Any parameter name-value pair given below.
%
% Parameters
%   SC              in, required, type=char
%   TSTART          in, required, type=char
%   TEND            in, required, type=char
%   'Directory'     in, optional, type=boolean, default=pwd();
%                   Directory in which to find sun pulse data.
%   'UniquePackets' in, optional, type=boolean, default=false
%                   If set, information from unique packets is returned.
%                     Packets can can contain overlapping data.
%   'UniquePulse'   in, optional, type=boolean, default=false
%                   Return unique sun pulse times. Sun pulse times are
%                     reported multiple times per packet, even if unchanged.
%                     Setting this parameter true automatically sets
%                     'UniqPackets' to true.
%
% Returns
%   SUNPULSE        out, required, type=struct
%                   Fields are:
%                     'Epoch'    - Packet times
%                     'SunPulse' - Sun pulse times
%                     'Period'   - Period (micro-sec) of revolution. Only
%                                    returned when FLAG=0 and only on the
%                                    second and subsequent received sun
%                                    pulses from the s/c.
%                     'Flag'     - Status flag
%                                    0: s/c sun pulse
%                                    1: s/c pseudo sun pulse
%                                    2: s/c CIDP generated speudo sun pulse
%
% Examples
%   Read data from mms3 on 2015-04-18 between 13:00 and 15:00:
%     >> sc       = 'mms3';
%     >> tstart   = '2015-04-18T13:00:00';
%     >> tend     = '2015-04-18T15:00:00';
%     >> sunpulse = mms_hk_read_sunpulse(sc, tstart, tend, 'UniquePackets', true);
%        sunpulse = Epoch: [1x7200 int64]
%                    Flag: [1x7200 uint8]
%                  Period: [1x7200 uint8]
%                SunPulse: [1x7200 int64]
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-03-25      Written by Matthew Argall
%   2015-04-12      Added the 'UniqPackets' and 'UniqPulse' parameters. - MRA
%   2015-04-14      Return a structure. - MRA
%   2015-04-19      Inputs are more general. Search for file on file system.
%                     Removed call to cdflib.inquireVar, as it does not like
%                     CDF_TIME_TT2000 variables. - MRA
%
function [sunpulse] = mms_hk_read_sunpulse(sc, tstart, tend, varargin)

	% Defaults
	uniq_packets = false;
	uniq_pulse   = false;
	sunpulse_dir = pwd();

	nOptArgs = length(varargin);
	for ii = 1 : 2 : nOptArgs
		switch varargin{ii}
			case 'Directory'
				sunpulse_dir = varargin{ii+1};
			case 'UniquePackets'
				uniq_packets = varargin{ii+1};
			case 'UniquePulse'
				uniq_pulse   = varargin{ii+1};
			otherwise
				error( ['Parameter name not recognized: "' varargin{ii} '".'] );
		end
	end
	
	% Must get unique packets before determining unique sun pulse times.
	if uniq_pulse && ~uniq_packets
		uniq_packets = true;
	end

%------------------------------------%
% Find the File                      %
%------------------------------------%
	% Create a file pattern
	fpattern = mms_construct_filename(sc, 'fields', 'hk', 'l1b_101', ...
	                                  'Tokens',    true, ...
	                                  'Directory', sunpulse_dir);

	% MMS#_DEFATT_%Y%D_%Y%D.V00
	[filename, nFiles] = MrFile_Search( fpattern, ...
	                                    'TStart',    tstart, ...
	                                    'TEnd',      tend, ...
	                                    'TimeOrder', '%Y%M%d');
	
	% Make sure the file exists
	assert( nFiles > 0, ...
	        'Sunpulse file not found or does not exist.');

%------------------------------------%
% Read the Data                      %
%------------------------------------%
	
	% Variable names
	sunpulse_name  = [sc '_101_sunpulse'];
	flag_name      = [sc '_101_sunssps'];
	period_name    = [sc '_101_iifsunper'];
	
	% Read the data
	[sun_pulse, hk_tt2000] = MrCDF_nRead(filename, sunpulse_name, ...
	                                     'sTime',       tstart, ...
	                                     'eTime',       tend, ...
	                                     'ColumnMajor', true);
	flag = MrCDF_nRead(filename, flag_name, ...
	                   'sTime',       tstart, ...
	                   'eTime',       tend, ...
	                   'ColumnMajor', true);
	period = MrCDF_nRead(filename, period_name, ...
	                     'sTime',       tstart, ...
	                     'eTime',       tend, ...
	                     'ColumnMajor', true);

%------------------------------------%
% Unique Data                        %
%------------------------------------%
	% Return only the unique packets?
	%   - Can be overlap with hk packets
	if uniq_packets
		% Unique packet times (returns sorted order)
		[hk_tt2000, iUniq, iHK] = unique(hk_tt2000);
		
		% Associated data
		sun_pulse = sun_pulse(iUniq);
		flag      = flag(iUniq);
		period    = period(iUniq);
	end
	
	% Return unique sun pulse times
	%   - Sun pulse is reported multiple times per packet, even if unchanaged.
	if uniq_pulse
		% Select the unique sun pulse times
		[sun_pulse, iUniq, iPulse] = unique(sun_pulse);
		
		% Associated data
		flag   = flag(iUniq);
		period = period(iUniq);
	end

%------------------------------------%
% Output                             %
%------------------------------------%
	% Create a structure of the output.
	sunpulse = struct( 'Epoch',    hk_tt2000, ...
	                   'Flag',     flag,      ...
	                   'Period',   period,    ...
	                   'SunPulse', sun_pulse );
end