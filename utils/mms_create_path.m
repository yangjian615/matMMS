%
% Name
%   mms_create_path
%
% Purpose
%   Create an SDC-like file path. The path will have the form:
%
%     root/sc/instr/mode/level[/optdesc]/year/month[/day]
%
% Calling Sequence
%   THEPATH = mms_create_path(ROOT, SC, INSTR, MODE, LEVEL, TSTART);
%     Create a file path THEPATH beginning in root directory ROOT with
%     spacecraft identifier SC, instrument identifier INSTR, telemetry
%     mode MODE, data level LEVEL, and data start time TSTART.
%
%   [__] = mms_create_path(..., OPTDESC);
%     Provide the file name/path optional descriptor, if it exists.
%
% Parameters
%   SC:             in, required, type=char
%   INSTR:          in, required, type=char
%   MODE:           in, required, type=char
%   LEVEL:          in, required, type=char
%   TSTART:         in, required, type=char, default='%Y%M%d' or '%Y%M%d%H%m%S'
%   OPTDESC         in, optional, type=char, default=''
%
% MATLAB release(s) MATLAB 8.2.0.701 (R2013b)
% Required Products None
%
% History:
%   2016-04-01      Written by Matthew Argall
%
function thePath = mms_create_path(root, sc, instr, mode, level, tstart, optdesc)

	% Defaults
	if nargin < 6
		tstart = '';
	end
	if nargin < 7
		optdesc  = '';
	end
	
	% Default TStart
	%   - Default to using tokens
	%   - FPI srvy files resemble brst files
	%   - EDP srvy files resemble brst files except L2 DCE
	if isempty(tstart)
		if strcmp( mode, 'brst' )
			tstart = '%Y%M%d%H%m%S';
		else
			switch instr
				case 'cal'
					tstart = ''
				case 'edp'
					if strcmp(optdesc, 'dce') && strcmp(level, 'l2')
						tstart = '%Y%M%d%H%m%S';
					else
						tstart = '%Y%M%d%H%m%S';
					end
				case 'fpi'
					tstart = '%Y%M%d%H%m%S';
				otherwise
					tstart = '%Y%M%d';
			end
		end
	end
	

	% Year, month and day directories
	if tstart(1) == '%'
		year  = tstart(1:2);
		month = tstart(3:4);
		day   = tstart(5:6);
	else
		tvec  = mms_parse_time(tstart);
		year  = tvec{1};
		month = tvec{2};
		day   = tvec{3};
	end

	% Do not use the day directory if we are in burst mode
	if ~strcmp(mode, 'brst')
		day = '';
	end
	
	% Create the file path
	thePath = fullfile(root, sc, instr, mode, level, optdesc, year, month, day);
end