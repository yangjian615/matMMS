%
% Name
%   mms_fsm_create_l1b
%
% Purpose
%   Merge MMS fluxgate and search coil magnetometer data in the frequency
%   domain.
%
% Calling Sequence
%   [B_MERGED, T_MERGED] = mms_fsm_merge(SC. TSTART, TEND);
%     Given the MMS spacecraft number SC (e.g. 'mms1'), and a time
%     interval, [TSTART, TEND), gather all of the required search coil and
%     fluxgate magnetometer data and merge them into a single dataset
%     B_MERGED with time stamps T_MERGED.
%
%   [..., B_FG, T_FG] = mms_fsm_merge(__);
%     Also return the calibrated FGM magnetic field B_FG and its time
%     stamps T_FG.
%
%   [..., B_SC, T_SC] = mms_fsm_merge(__);
%     Also return the *UN*calibrated SCM magnetic field B_SC and its time
%     stamps T_SC.
%
% Parameters
%   SC:             in, required, type=char
%   TSTART:         in, required, type=char
%   TEND:           in, required, type=char
%   'Duration':     in, required, type=double
%                   The duration of each merging interval. Sets the
%                     frequency resolution of the final dataset.
%   'f_max':        in, required, type=double, default=Nyquist frequency
%                   The maximum of the frequency range to merge.
%   'f_min':        in, required, type=double, default=df
%                   The minimum ( > 0 ) of the frequency range to merge.
%   'fg_dir':       in, required, type=char, default=pwd();
%                   Directory in which to find FGM data.
%   'fg_cal_dir':   in, required, type=char, default=pwd();
%                   Directory in which to find FGM calibration data.
%   'sc_dir':       in, required, type=char, default=pwd();
%                   Directory in which to find SCM data.
%   'sc_cal_dir':   in, required, type=char, default=pwd();
%                   Directory in which to find SCM calibration data.
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-12-08      Written by Matthew Argall
%
function scm_files = mms_scm_l2_process(sc, mode, date_start, date_end)

	% Check inputs
	if nargin < 1 || isempty(sc)
		sc = {'mms1', 'mms2', 'mms3', 'mms4'};
	end
	if nargin < 2 || isempty(sc)
		mode = {'srvy', 'brst'};
	end
	if nargin < 3 || isempty(date_start)
		date_start = datestr(now()-3, 'yyyymmdd');
	end

	if ischar(sc)
		sc = { sc };
	end
	if ischar(mode)
		mode = { mode };
	end
	
	
	% If we are processing survey data, we must merge fast and
	% slow separately, then combine them after.
	[tf_srvy, ind] = ismember('srvy', mode);
	if tf_srvy
		mode{ind} = [];
		mode      = [ mode 'fast' 'slow'];
	end

	% How much to process
	nsc   = length(sc);
	nmode = length(mode);

%------------------------------------%
% Handle Inputs                      %
%------------------------------------%
	% Parse input times
	start_vec = mms_parse_time(date_start);
	end_vec   = mms_parse_time(date_end);
	
	% Convert to ISO
	tstart = datestr(start_vec, 'yyyy-mm-ddTHH:MM:SS');
	tend   = datestr(end_vec,   'yyyy-mm-ddTHH:MM:SS');

%------------------------------------%
% Each Mode and SC                   %
%------------------------------------%
	for ii = 1 : nsc
	for jj = 1 : nmode
		if strcmp(mode{jj}, 'brst')
			timeorder = '%Y%M%d%H%m%S';
			subdirs = {'%Y', '%M', '%d'};
		else
			timeorder = '%Y%M%d';
			subdirs = {'%Y', '%M'};
		end

	%------------------------------------%
	% Find Attitude & DSS File           %
	%------------------------------------%
		sdc_root = '/nfs';
		att_dir  = fullfile(sdc_root, 'ancillary', sc{ii}, 'defatt');
		hk_dir   = fullfile(sdc_root, 'hk');
	
		% Attitude
		att_ftest = fullfile( att_dir, [upper(sc{ii}) '_DEFATT_%Y%D_%Y%D.V*'] );
		[att_file, count] = MrFile_Search(att_ftest,              ...
		                                  'Closest',      true,   ...
		                                  'TimeOrder',    '%Y%D', ...
		                                  'TStart',       tstart, ...
		                                  'TEnd',         tend, ...
		                                  'VersionRegex', 'V[0-9]{2}');
		assert(count > 0, ['No definitive attitude file found: "' att_ftest '".']);
	
		% DSS Data File
		[dss_file, count, str] = mms_file_search(sc{ii}, 'fields', 'hk', 'l1b', ...
		                                         'SDCroot', hk_dir,         ...
		                                         'OptDesc', '101',          ...
		                                         'TStart',  tstart,         ...
		                                         'TEnd',    tend);
		assert(count > 0, ['DSS file not found: "' str '".']);

	%------------------------------------%
	% Find SCM File                      %
	%------------------------------------%
		scm_cal_dir = '/home/argall/data/mms/scm_cal/';
		
		if strcmp(mode, 'fast')
			duration = 64.0;
			optdesc  = 'scf';
		elseif strcmp(mode, 'slow')
			duration = 64.0;
			optdesc = 'scs';
		elseif strcmp(mode, 'brst')
			duration = 2.0
			optdesc = 'scb';
		else
			error( ['Unexpected mode: "' mode '".'] )
		end
	
		% SCM Data File
		[scm_file, count, str] = mms_file_search(sc{ii}, 'scm', mode{jj}, 'l1a', ...
		                                         'SubDirs',   subdirs,     ...
		                                         'OptDesc',   optdesc,     ...
		                                         'TimeOrder', timeorder,   ...
		                                         'TStart',    tstart,      ...
		                                         'TEnd',      tend);
		assert(count > 0, ['SCM file not found: "' str '".']);
	
		% Determine the flight model
		switch sc{ii}
			case 'mms1'
				fm = 'scm1';
			case 'mms2'
				fm = 'scm2';
			case 'mms3'
				fm = 'scm4';
			case 'mms4'
				fm = 'scm3';
			otherwise
				error(['Invalid spacecraft ID: "' sc{ii} '".'])
		end
	
		% SCM Cal File
		scm_ftest = fullfile(scm_cal_dir, [sc{ii} '_' fm '_caltab_%Y%M%d%H%M%S_v*.txt']);
		[scm_cal_file, count] = MrFile_Search(scm_ftest, ...
		                                      'Closest',      true, ...
		                                      'TimeOrder',    '%Y%M%d%H%M%S', ...
		                                      'TStart',       tstart, ...
		                                      'TEnd',         tend, ...
		                                      'VersionRegex', 'v[0-9]');
		assert(count > 0, ['No SCM calibration file found: "' scm_ftest '".']);
	

	%------------------------------------%
	% Process Data                       %
	%------------------------------------%
		[att, hdr] = mms_fdoa_read_defatt(att_file);
		sunpulse   = mms_dss_read_sunpulse(dss_file, '', '', 'UniquePulse', true);
		[t, b_gsm, b_gse, ~, b_dmpa, b_bcs, ~, b_omb, ~] ...
			= mms_scm_l2_create(scm_file, scm_cal_file, tstart, tend, ...
			                    'Attitude', att, ...
			                    'Duration', duration, ...
			                    'SunPulse', sunpulse, ...
			                    'zMPA',     hdr.zMPA(:,1));
		
		clear att hdr sunpulse
		
	%------------------------------------%
	% Write Data to File                 %
	%------------------------------------%
		% Structure
		scm = struct( 'tt2000', t',      ...
		              'b_omb',  single(b_omb)',  ...
		              'b_bcs',  single(b_bcs)',  ...
		              'b_dmpa', single(b_dmpa)', ...
		              'b_gse',  single(b_gse)',  ...
		              'b_gsm',  single(b_gsm)'   ...
		            );
		clear t b_omb b_bcs b_dmpa b_gse b_gsm
		
		% Parent files
try
		parents      = { scm_file{:} scm_cal_file att_file{:} dss_file{:} };
		[~, parents] = cellfun(@fileparts, parents, 'UniformOutput', false);
catch ME
	keyboard
end
		
		% Write to file
		optdesc   = sprintf('-%ds', duration);
		scm_files = mms_scm_l2_write( parents, scm, optdesc );
		
		if nargout() == 0
			disp(['SCM file written to "' scm_file '".'])
		end
	end
	end

	disp('Finished!')
end