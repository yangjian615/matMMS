%
% Construct an MMS file name. The file name format is:
%
%   scId_instrumentId_mode_dataLevel_optionalDataProductDescriptor_startTime_vX.Y.Z.cdf
%
% :Examples:
%   Dissect and construct a filename
%     filename = 'mms3_dfg_hr_l2p_duration-1h1m_20010704_v1.0.0.cdf';
%     [sc, instr, mode, level, desc, tstart, version] = mms_dissect_filename(filenames);
%     fname = mms_construct_filename(sc, instr, mode, level, desc, tstart, version)
%     fname =
%              mms3_dfg_hr_l2p_duration-1h1m_20010704_v1.0.0.cdf
%
% Calling Sequence:
%   FNAME = mms_construct_filename(SC, INSTR, MODE, LEVEL)
%     Construct an MMS file name FNAME using the spacecraft ID, SC,
%     instrumnet identifier INSTR, telemetry mode MODE and data level
%     LEVEL. Other mandatory components of the file name will be the
%     wildcard character '*' and optional components will be the empty
%     string ''.
%
%   FNAME = mms_construct_filename(..., 'ParamName', ParamValue)
%     Any parameter name and value described below.
%
% :Params:
%   SC:                 in, required, type=char
%                       mms1, mms2, mms3, mms4
%   INSTR:              in, required, type=char
%                       Instrument or investigation identifier
%                           hpca
%                           aspoc
%                           epd
%                           epd-eis
%                           epd-feeps
%                           fpi
%                           des
%                           dis
%                           des-dis
%                           fields
%                           edi
%                           adp
%                           sdp
%                           adp-sdp
%                           afg
%                           dfg
%                           dsp
%                           afg-dfg
%                           scm
%   MODE:               in, required, type=string
%                       Instrument telemetry mode:
%                           fast
%                           slow
%                           brst
%                           srvy
%   LEVEL:              in, required, type=string
%                       Level of the data product:
%                           l1a
%                           l1b
%                           l2
%                           ql
%                           l2pre
%                           l2plus
%   'Directory':        in, required, type=char
%                       A directory to be appended to the file name.
%   'OptDesc':          in, optional, type=char, default=''
%                       Optional data product descriptor. Should be short
%                           (3-8 characters). Hyphens used to separate
%                           multiple components.
%   'TStart':           in, required, type=char
%                       Start time of the data product, formatted as:
%                           'yyyymmddhhmmss'. Least significant fields can
%                           be dropped when files start on regular hourly
%                           or minute boundaries.
%   'Version':          in, required, type=char
%                       Version number in the form: "vX.Y.Z"
%                           X - Interface number. Increments represent
%                               significant changes that will break code or
%                               require code changes in analysis software.
%                           Y - Quality number. Represents change in
%                               quality of the, such as calibration or
%                               fidelity. Should not impact software.
%                           Z - Bug fix/Revision number. Minor changes to
%                               the contents of the file due to
%                               reprocessing of missing data. Dependent
%                               data products should be reprocessed.
%
function fname = mms_construct_filename(sc, instr, mode, level, varargin)
	
%------------------------------------%
% Inputs                             %
%------------------------------------%
	
	% Defaults
	tokens    = false;
	tstart    = '';
	optDesc   = '';
	directory = '';
	version   = '';

	% Check for optional arguments
	nOptArgs = length(varargin);
	for ii = 1 : 2 : nOptArgs
		switch varargin{ii}
			case 'Directory'
				directory = varargin{ii+1};
			case 'OptDesc'
				optDesc = varargin{ii+1};
			case 'Tokens'
				tokens = varargin{ii+1};
			case 'Version'
				version = varargin{ii+1};
			case 'TStart'
				tstart = varargin{ii+1};
			otherwise
				error( ['Unknown parameter "' varargin{ii} '".'] );
		end
	end

	% Use MrToken tokens or a wildcard character?
	if isempty(tstart)
		if tokens
			if strcmp(mode, 'brst')
				tstart = '%Y%M%d%H%m%S';
			else
				switch instr
					% FPI
					case 'fpi'
						tstart = '%Y%M%d%H%m%S';
			
					% FSM
					case 'fsm'
						if ~isempty( regexp(optdesc, '^cal-(afg|dfg|scm)$') )
							tstart = '%Y%M%d%H%m%S';
						else
							tstart = '%Y%M%d';
						end
			
					% EDP
					case 'edp'
						if strcmp(optdesc, 'dce')
							if strcmp(level, 'l2')
								tstart = '%Y%M%d%H%m%S';
							else
								tstart = '%Y%M%d';
							end
						else
							tstart = '%Y%M%d';
						end
				
					% Everthing else
					otherwise
						tstart = '%Y%M%d';
				end
			end
		else
			tstart = '*';
		end
	end
	
	% Use a wildcard to represent the version.
	if isempty(version)
		version = '*';
	end

	% Separate the optional descriptor from the start time.
	if ~isempty(optDesc)
		optDesc = [optDesc '_'];
	end

%------------------------------------%
% Create File Name                   %
%------------------------------------%

	% Construct the file name.
	fname = strcat(sc, '_', instr, '_', mode, '_', level, '_', optDesc, tstart, '_v', version, '.cdf');

	% Return a single string if only one file name was made
	if iscell(fname) && length(fname) == 1
		fname = fname{1};
	end

	% Prepend a directory?
	if ~isempty(directory)
		fname = fullfile(directory, fname);
	end
end