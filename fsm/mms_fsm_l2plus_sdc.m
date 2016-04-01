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
function fsm_file = mms_fsm_l2plus_sdc(fgm_file, scm_file, scm_cal_file, att_file, dss_file, duration)

	% Dissect the searchcoil file
	[sc, ~, mode, ~, tstart] = mms_dissect_filename(scm_file);
	if strcmp(mode, 'brst')
		duration = 10;
	else
		duration = 150;
	end

%------------------------------------%
% Read Attitude and zMPA             %
%------------------------------------%
	[defatt, att_hdr] = mms_fdoa_read_defatt(att_file);
	zmpa              = att_hdr.zMPA(:,1);

%------------------------------------%
% Process Data                       %
%------------------------------------%
	[t_fsm, b_omb] = mms_fsm_l2plus_create(fgm_file, scm_file, scm_cal_file, zmpa, duration);

%------------------------------------%
% SMPA & BCS                         %
%------------------------------------%

	% Rotation Matrices
	omb2smpa = mms_fg_xomb2smpa();
	bcs2smpa = mms_fg_xbcs2smpa(zmpa);
	smpa2bcs = permute(bcs2smpa, [2,1]);
	
	% OMB --> SMPA --> BCS
	b_smpa = mrvector_rotate(omb2smpa, b_omb);
	b_bcs  = mrvector_rotate(smpa2bcs, b_smpa);

	% Delete data
	clear omb2smpa bcs2smpa b_omb

%------------------------------------%
% Despin                             %
%------------------------------------%
	% Read DSS data
	sunpulse = mms_dss_read_sunpulse(dss_file, '', '', 'UniquePulse', true);
	
	% Rotation matrix
	smpa2dmpa = mms_dss_xdespin(sunpulse, t_fsm);
	
	% SMPA --> DMPA
	b_dmpa = mrvector_rotate(smpa2dmpa, b_smpa);

	% Delete data
	clear sunpulse smpa2dmpa

%------------------------------------%
% DMPA --> GSE & GSM                 %
%------------------------------------%

	% Rotation Matrix to GEI
	[b_gsm, b_gse] = mms_rot_despun2gsm(t_fsm, b_dmpa, defatt);
	
	% Delete data
	clear defatt
	
%------------------------------------%
% Write to File                      %
%------------------------------------%
	% Parent files
try
	parents      = { fgm_file scm_file scm_cal_file att_file{:} dss_file };
	[~, parents] = cellfun(@fileparts, parents, 'UniformOutput', false);
catch ME
	keyboard
end

	% Create a data structure
	data = struct( 'tt2000', t_fsm',  ...
	               'b_bcs',  single(b_bcs)',  ...
	               'b_dmpa', single(b_dmpa)', ...
	               'b_gse',  single(b_gse)',  ...
	               'b_gsm',  single(b_gsm)'   ...
	             );
	clear t_fsm b_bcs b_smpa b_gse b_gsm

	% Write the file
	fsm_file = mms_fsm_l2plus_write( parents, data );
end