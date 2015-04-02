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
% 				fname =
% 				mms3_dfg_hr_l2p_duration-1h1m_20010704_v1.0.0.cdf
%
% Calling Sequence:
%   FNAME = mms_construct_filename(SC, INSTR, MODE, LEVEL)
%     Construct an MMS file name FNAME using the spacecraft ID, SC,
%     instrumnet identifier INSTR, telemetry mode MODE and data level
%     LEVEL. Other mandatory components of the file name will be the
%     wildcard character '*' and optional components will be the empty
%     string ''.
%
%   FNAME = mms_construct_filename(..., TSTART)
%     Specify the time at which the data begins TSTART. The time is
%     formatted as 'YYYYMMDD' with optional hour, minute, and second
%     fields, if necessary.
%
%   FNAME = mms_construct_filename(..., VERSION)
%     Specify the version of the file: 'X.Y.Z'.
%
%   FNAME = mms_construct_filename(..., DESC)
%     Add an optional descriptor. Parts should be separated by the hyphen
%     '-'.
%
%   FNAME = mms_construct_filename(..., DIRECTORY)
%     Indicate a directory DIRECTORY that should be prepended to the file
%     name.
%
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
%   DESC:               in, optional, type=char, default=''
%                       Optional data product descriptor. Should be short
%                           (3-8 characters). Hyphens used to separate
%                           multiple components.
%   TSTART:             in, required, type=char
%                       Start time of the data product, formatted as:
%                           'yyyymmddhhmmss'. Least significant fields can
%                           be dropped when files start on regular hourly
%                           or minute boundaries.
%   VERSION:            in, required, type=char
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
function fname = mms_construct_filename(sc, instr, mode, level, tstart, version, desc, directory)


%------------------------------------%
% Defaults                           %
%------------------------------------%
	nin = nargin();
	if nin < 8
		directory = '';
	end
	if nin < 7 || isempty(desc)
		desc = '';
	end
	if nin < 6
		version = '*';
	end
	if nin < 5
		tstart = '*';
	end

%------------------------------------%
% Create File Name                   %
%------------------------------------%

	% Separate the optional descriptor from the start time.
	if ~isempty(desc)
		desc = [desc '_'];
	end

	% Construct the file name.
	fname = [sc '_' instr '_' mode '_' level '_' desc tstart '_' version '.cdf'];
	
	% Prepend a directory?
	if ~isempty(directory)
		fname = fullfile(directory, fname);
	end
end