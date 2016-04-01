%
% Name
%   mms_fdoa_epoch2tt2000
%
% Purpose
%   Parse information for the header of MMS definitive attitude files.
%
% Calling Sequence
%   TT2000 = mms_fdoa_epoch2tt2000(ATT_EPOCH)
%     Convert an FDOA UTC string or cell array of strings to CDF TT2000
%     times.
%
%   TT2000 = mms_fdoa_epoch2tt2000(__, 'ParamName', ParamValue)
%     Any of the parameter name-value pairs below.
%
% Parameters
%   ATT_EPOCH       in, required, type=char/cell
%   'AttTAI'        in, optional, type=boolean, default=false
%                   Indicate that Attitude TAI times were given. These are
%                     seconds since the MMS reference epoch.
%   'EphemTAI'      in, optional, type=boolean, default=false
%                   Indicate that Ephemeris TAI times were given. These are
%                     days since the MMS reference epoch.
%
% Returns
%   TT2000          out, required, type=int64 (cdf_time_tt2000)
%
% Examples
%   Convert UTC to TT2000
%     >> tt2000 = mms_fdoa_epoch2tt2000( '2015-079T00:59:04.312' )
%       tt2000 = 480085211496000000
%     >> utc = spdfencodett2000(tt2000)
%       utc = '2015-03-20T00:59:04.312000000'
%
%   Convert TAI to TT2000
%     >> tt2000 = mms_fdoa_epoch2tt2000( 1805504379.312, 'AttTAI', true )
%       tt2000 = 480085211496000000
%     >> utc = spdfencodett2000(tt2000)
%       utc = '2015-03-20T00:59:04.312000000'
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-09      Written by Matthew Argall
%   2015-04-10      Attempt to speed up parsing epoch strings. - MRA
%   2015-04-11      Converted TAI parameter to 'AttTAI' and 'EphemTAI'. - MRA
%
function tt2000 = mms_fdoa_epoch2tt2000(fdoa_epoch, varargin)
	%
	% COMMENT           The time tags are given twice. First as calendar UTC time,
	% COMMENT           and then as TAI, expressed as elapsed SI seconds since
	% COMMENT           the mission reference epoch, 1958-001T00:00:00 UTC.
	%
	% COMMENT   Time (UTC)    Elapsed Sec
	% 2015-079T00:59:04.312 1805504379.312
	%
	
	% Defaults
	tai_ephem = false;
	tai_att   = false;
	
	nOptArg = length(varargin);
	for ii = 1 : 2 : nOptArg
		switch varargin{ii}
			case 'EphemTAI'
				tai_ephem = varargin{ii+1};
			case 'AttTAI'
				tai_att = varargin{ii+1};
			otherwise
				error( ['Optional parameter not recognized: ' varargin{ii} '".'] );
		end
	end
	
	% Cannot use AttTAI and EphemTAI together
	assert( ~tai_ephem || ~tai_att, 'EphemTAI and AttTAI cannot both be true.' );

%------------------------------------%
% TAI to TT2000                      %
%------------------------------------%
	if tai_ephem || tai_att
		% Convert to double precision
		if ischar(fdoa_epoch)
			fdoa_epoch = str2double(fdoa_epoch);
		elseif iscell(fdoa_epoch)
			fdoa_epoch = cellfun(@str2double, fdoa_epoch);
		end
		
		% Ephemeris epoch is given in number of days.
		%   - Convert to seconds.
		if tai_ephem
			fdoa_epoch = fdoa_epoch * 86400.0;
		end
		
		% Convert 1958-001T00:00:00 UTC to TT2000 (TAI)
		ref_epoch = MrCDF_Epoch_Parse('1958-01-01T00:00:00.000000000', 'CDF_TIME_TT2000');
		
		% Convert TAI seconds to nanoseconds
		tai = int64(fdoa_epoch * 1e9);
		
		% Convert to TT2000
		tt2000 = tai + ref_epoch;

%------------------------------------%
% UTC to TT2000                      %
%------------------------------------%
	% Convert the UTC string
	else
		% Make sure a cell or string was given
		assert(iscell(fdoa_epoch) || ( ischar(fdoa_epoch) && isrow(fdoa_epoch) ), ...
		       'ATT_EPOCH must be a string or cell array of strings')
		
		%
		% Expected:
		%   yyyy-mm-ddTHH:MM:SS.mmmuuunnn
		%
		% Times:
		%   - Attitude:  2015-079T00:59:04.312
		%   - Ephemeris: 2015-073/03:00:25.000
		%
		
		% Find missing number of zeros.
		if iscell(fdoa_epoch)
			len = length( fdoa_epoch{1} );
		else
			len = length( fdoa_epoch );
		end
		nZeros = 29 - len;
		zz     = repmat('0', 1, nZeros);
		
		% Separate year, doy, and time
		epoch_char   = char(fdoa_epoch);
		year         = epoch_char(:,  1:4);
		doy          = epoch_char(:,  6:8);
		time         = epoch_char(:, 10:end);
		
		% Compuate month and day
		[month, day] = MrDOY2MonthDay( str2num(doy), str2num(year) );
		clear epoch_char doy
		
		% Gather pieces into cell arrays of strings
		year  = cellstr(year);
		time  = cellstr(time);
		month = num2str(month, '%02d');
		day   = num2str(day,   '%02d');
		
		% UTC
		utc = strcat(year, '-', month, '-', day, 'T', time, zz);
		
		% Convert to tt2000
		tt2000 = spdfparsett2000(utc);
	end
end