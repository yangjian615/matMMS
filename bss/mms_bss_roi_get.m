%
% Name
%   mms_bss_roi_read
%
% Purpose
%   Read MMS region of interest (ROI) start and end times from a file. ROI times can
%   be obtained from IDL and SPEDAS using:
%       IDL> print, time_string( mms_get_roi( ['2015-09-18/00:00:00', '2015-09-18/24:00:00'] ) )
%
%   Or can be obtained and written to a file using
%       IDL> unh_bss_roi_write()
%
%   Times are returned as CDF TT2000 times.
%
% Calling Sequence
%   TT2000 = mms_bss_roi_read( TSTART );
%     Search for an ROI that contains the time TSTART. If TSTART lies outside of an ROI,
%     the very next ROI will be selected. TSTART should be a date-time string formatted
%     as 'yyyy-mm-ddTHH:MM:SS'
%
%   TT2000 = mms_bss_roi_read( TSTART, TEND );
%     Search for ROIs contained in the interval [TSTART, TEND]. TEND should be formatted
%     the same as TSTART.
%
%   TT2000 = mms_bss_roi_read( TSTART, TEND, 'ParamName', ParamValue );
%     Modify the ROI selections using parameter name-value pairs, as defined below.
%
% Parameters
%   TSTART:             in, required, type=char
%   TEND:               in, required, type=char
%   'Filename':         in, optinoal, type=char, default='/home/argall/mms_roi_times.txt'
%                       Name of the ASCII file containing ROI start and end times.
%   'Next':             in, optinoal, type=boolean, default=false
%                       If true, the search for the ROI will be performed as normal,
%                           then the ROI after the one found will be returned.
%   'Orbit':            in, optinoal, type=boolean, default=false
%                       If true, the search for the ROI will be performed as normal,
%                           then the start time will moved to the end of the previous
%                           ROI, providing one orbit, starting at the beginning of slow
%                           survey and ending at the end of fast survey.
%   'Previous':         in, optinoal, type=boolean, default=false
%                       If true, the search for the ROI will be performed as normal,
%                           then the ROI before the one found will be returned.
%   'SlowSurvey':       in, optinoal, type=boolean, default=false
%                       If true, the search for the ROI will be performed as normal,
%                           then the interval between the ROI found and the previous ROI
%                           will be returned. (The data rate mode outside the ROI is
%                           named "Slow Survey").
%
% Returns
%   TT2000:             in, required, type=int64 Nx2
%
% MATLAB release(s) MATLAB 8.2.0.701 (R2013b)
% Required Products None
%
% History:
%   2016-05-20      Written by Matthew Argall
%   2016-07-15      Added the "Orbit" optional parameter. Providing TEND returns
%                       all ROI times between TSTART and TEND, not just the first. - MRA
%
function t_roi = mms_bss_roi_get( tstart, tend, varargin )

	% Read in ROI times only once
	global bss_roi_times
	
%------------------------------------%
% Defaults                           %
%------------------------------------%

	% Initial Values
	if nargin() < 2
		tend = '';
	end
	
	filename    = '';
	tf_next     = false;
	tf_previous = false;
	tf_slow     = false;
	tf_orbit    = false;
	
	% Check inputs for optinoal arguments
	nOptArgs = length(varargin);
	for ii = 1 : 2 : nOptArgs
		switch varargin{ii}
			case 'Filename'
				filename = varargin{ii+1};
			case 'Next'
				tf_next = varargin{ii+1};
			case 'Previous'
				tf_previous = varargin{ii+1};
			case 'SlowSurvey'
				tf_slow = varargin{ii+1};
			case 'Orbit'
				tf_orbit = varargin{ii+1};
			otherwise
				error(['Optional argument not recognized: "' varargin{ii} '".']);
		end
	end
	
	% Restrictions
	assert( (tf_previous + tf_next)  ~= 2, '"Next" and "Previous" are mutually exclusive.' )
	assert( (tf_slow     + tf_orbit) ~= 2, '"Slow" and "Orbit" are mutually exclusive.' )
	if ~isempty(tend)
		assert( ~tf_previous && ~tf_next, '"Previous" and "Next" cannot be used with TEND.')
	end
	
%------------------------------------%
% Check Inputs                       %
%------------------------------------%
	
	% Is TSTART a string?
	if ischar(tstart)
		assert( MrTokens_IsMatch(tstart, '%Y-%M-%dT%H:%m:%S'), 'TSTART not formatted correctly.');
		t0_tt2000 = MrCDF_Epoch_Parse([tstart '.000000000'], 'CDF_TIME_TT2000');
	elseif isa(tstart, 'int64')
		t0_tt2000 = tstart;
	else
		error('TSTART must be class char or int64');
	end

	% Is TEND a string?
	if ischar(tend)
		if ~isempty(tend)
			assert( MrTokens_IsMatch(tend, '%Y-%M-%dT%H:%m:%S'), 'TEND not formatted correctly.');
			t1_tt2000 = MrCDF_Epoch_Parse([tend '.000000000'], 'CDF_TIME_TT2000');
		end
	elseif isa(tstart, 'int64')
		t1_tt2000 = tend;
	else
		error('TEND must be class char or int64');
	end
	
	% Filename
	if isempty(filename)
		filename = '/home/argall/mms_roi_times.txt';
	end
	
%------------------------------------%
% Read File                          %
%------------------------------------%
	
	% Read the file
	if isempty(bss_roi_times)
		bss_roi_times = mms_bss_roi_read(filename);
	end
	
%------------------------------------%
% Select ROI                         %
%------------------------------------%
	% Only TSTART is given:
	if isempty(tend)
		% If TSTART lies within an ROI, pick that ROI
		iROI = find(bss_roi_times(1,:) <= t0_tt2000 & ...
		            bss_roi_times(2,:) >= t0_tt2000, 1, 'first');
		
		% If the time does not lie within an ROI, select the next one
		if isempty(iROI)
			iROI = find(bss_roi_times(1,:) >= t0_tt2000, 1, 'first');
		end
		
		% Adjust the ROI selection
		if tf_next
			iROI = iROI + 1;
		elseif tf_previous
			iROI = iROI - 1;
		end
	
	% TSTART and TEND are given
	else
		% Find ROIs containing both
		iROI = find( bss_roi_times(1,:) <= t1_tt2000 & ...
		             bss_roi_times(2,:) >= t0_tt2000 );
	end
	
%------------------------------------%
% Select ROI                         %
%------------------------------------%
	% ROI times
	if tf_slow
		t_roi = [ bss_roi_times(2,iROI-1); bss_roi_times(1,iROI) ];
	elseif tf_orbit
		t_roi = [ bss_roi_times(2,iROI-1); bss_roi_times(2,iROI) ];
	else
		t_roi = bss_roi_times(:,iROI);
	end
end