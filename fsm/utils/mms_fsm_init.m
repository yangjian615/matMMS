%
% Name
%   mms_fsm_init
%
% Purpose
%   Create a set of global variables specifying specific data paths.
%
% Calling Sequence
%   mms_fsm_init();
%     Establish data processing paths as global variables.
%
%   mms_fsm_init('ParamName', ParamValue);
%     Provide any of the parameter-name parameter-value pairs below.
%
% Parameters
%   'cal_path_root':    in, required, type=char, default='/nfs/edi/cals'
%                       Root of an SDC-like directory structure in which to find
%                         calibration data.
%   'data_path_root':   in, required, type=char, default='/nfs'
%                       Root of an SDC-like directory structure data data finds its
%                         final resting place.
%   'dropbox_root':     in, required, type=char, default='/nfs/fsm/temp'
%                       Root of an SDC-like directory structure where data is initially
%                         deposited. It is later moved to "data_path_root".
%   'hk_root':          in, required, type=char, default='/nfs/hk'
%                       Root of an SDC-like directory structure in which to find
%                         calibration data.
%   'unh_data_root':    in, required, type=char, default='/nfs/fsm'
%                       Root of an SDC-like directory structure at data is save locally.
%
% MATLAB release(s) MATLAB 8.2.0.701 (R2013b)
% Required Products None
%
% History:
%   2016-04-01      Written by Matthew Argall
%
function [] = mms_fsm_init( varargin )

	% Set global variables
	global cal_path_root data_path_root dropbox_root hk_root log_path_root unh_data_root

%------------------------------------%
% Defaults                           %
%------------------------------------%

	% Defaults
	cal_path_root  = '/nfs/fsm/cals';
	data_path_root = '/nfs';
	dropbox_root   = '/nfs/fsm/temp';
	hk_root        = '/nfs/hk';
	log_path_root  = '/nfs/fsm/logs';
	unh_data_root  = '/nfs/fsm';

%------------------------------------%
% Optional Arguments                 %
%------------------------------------%
	nOptArgs = length(varargin);
	for ii = 1 : 2 : nOptArgs
		switch vargin{ii}
			case 'cal_path_root'
				cal_path_root = varargin{ii+1};
			case 'data_path_root'
				data_path_root = varargin{ii+1};
			case 'dropbox_root'
				dropbox_root = varargin{ii+1};
			case 'hk_root'
				hk_root = varargin{ii+1};
			case 'log_path_root'
				log_path_root = varargin{ii+1};
			case 'unh_data_root'
				unh_data_root = varargin{ii+1};
			otherwise
				error('MMS:INIT', 'Input parameter not recognized: "%s"', varargin{ii});
		end
	end

%------------------------------------%
% Search Environment Variables       %
%------------------------------------%

	% Search for environment variables
	temp = getenv('CAL_PATH_ROOT');
	if ~isempty(temp)
		cal_path_root = temp;
	end
	
	temp = getenv('DATA_PATH_ROOT');
	if ~isempty(temp)
		data_path_root = temp;
	end
	
	temp = getenv('DROPBOX_ROOT');
	if ~isempty(temp)
		dropbox_root = temp;
	end
	
	temp = getenv('HK_ROOT');
	if ~isempty(temp)
		hk_root = temp;
	end
	
	temp = getenv('LOG_PATH_ROOT');
	if ~isempty(temp)
		log_path_root = temp;
	end
	
	temp = getenv('UNH_DATA_ROOT');
	if ~isempty(temp)
		unh_data_root = temp;
	end
end