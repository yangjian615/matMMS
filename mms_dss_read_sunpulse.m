%
% Name
%   mms_hk_read_sunpulse
%
% Purpose
%   Read houskeeping 0x101 sun pulse and related data.
%
% Calling Sequence
%   SUNPULSE = mms_hk_read_sunpulse(FILENAMES)
%     Read HK101 digital sun sensor data from files FILENAMES.
%
%   [__] = mms_hk_read_sunpulse(__, TSTART, TEND)
%     Select data between the time interval beginning at TSTART and
%     ending at TEND. TSTART and TEND should be formatted as
%     'yyyy-mm-ddThh:mm:ss'.
%
%   [__] = mms_hk_read_sunpulse(..., 'ParamName', ParamValue)
%     Any parameter name-value pair given below.
%
% Parameters
%   DSS_FILES       in, required, type=char
%   TSTART          in, required, type=char, default=''
%   TEND            in, required, type=char, default=''
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
%     >> dss_file = 'mms3_fields_hk_l1b_101_20150418_v0.3.4.cdf';
%     >> tstart   = '2015-04-18T13:00:00';
%     >> tend     = '2015-04-18T15:00:00';
%     >> sunpulse = mms_hk_read_sunpulse(dss_file, tstart, tend, 'UniquePackets', true);
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
%   2015-05-25      Take file names as inputs. - MRA
%
function sunpulse = mms_hk_read_sunpulse(dss_files, tstart, tend, varargin)

	% Defaults
	uniq_packets = false;
	uniq_pulse   = false;

	nOptArgs = length(varargin);
	for ii = 1 : 2 : nOptArgs
		switch varargin{ii}
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
% Check Files                        %
%------------------------------------%
	% Number of files given
	if iscell(dss_files)
		nFiles = length(dss_files);
	else
		assert(ischar(dss_files) && isrow(dss_files), 'DSS_FILES must be a file name or cell array of file names.')
		nFiles = 1;
	end

	% Dissect the file name
	[sc, instr, mode, level, ~, ~, optdesc] = mms_dissect_filename(dss_files);

	% AFG or DFG
	assert( min( strcmp(instr, 'fields') ) == 1, 'Only FIELDS files are allowed.' );

	% Instr, Level, Mode
	assert( min( strcmp(mode,    'hk') )  == 1, 'All files must be house-keeping (HK) files.' );
	assert( min( strcmp(level,   'l1b') ) == 1, 'Only L1B files are allowed.' );
	assert( min( strcmp(optdesc, '101') ) == 1, 'Only HK 101 files are allowed.' );

	% We now know all the files match, so keep on the the first value.
	if nFiles > 1
		sc      = sc{1};
		instr   = instr{1};
		mode    = mode{1};
		level   = level{1};
		optdesc = optdesc{1};
	end

%------------------------------------%
% Read the Data                      %
%------------------------------------%
	
	% Variable names
	sunpulse_name  = mms_construct_varname(sc, optdesc, 'sunpulse');
	flag_name      = mms_construct_varname(sc, optdesc, 'sunssps');
	period_name    = mms_construct_varname(sc, optdesc, 'iifsunper');
	
	% Read the data
	[sun_pulse, hk_tt2000] = MrCDF_nRead(dss_files, sunpulse_name, ...
	                                     'sTime',       tstart, ...
	                                     'eTime',       tend, ...
	                                     'ColumnMajor', true);
	flag = MrCDF_nRead(dss_files, flag_name, ...
	                   'sTime',       tstart, ...
	                   'eTime',       tend, ...
	                   'ColumnMajor', true);
	period = MrCDF_nRead(dss_files, period_name, ...
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

	% Filter period = 0
	izero = find(period == 0);
	if ~isempty(izero)
		warning('DSS:Read_SunPulse', 'Period = 0 found with flags %d. Removing.', flag(izero(1)));
		flag(izero)      = [];
		period(izero)    = [];
		sun_pulse(izero) = [];
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