%
% Name
%   mms_edi_view
%
% Purpose
%   Take l1a EDI  data, despin, and rotate into GSE. Then, plot the 
%   firing vectors and gun positions in BPP.
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-14      Written by Matthew Argall
%

%------------------------------------%
% Inputs                             %
%------------------------------------%
	sc_data_dir  = '/Users/argall/Documents/Work/Data/MMS/SCM/';
	sc_cal_dir   = '/Users/argall/Documents/Work/Data/MMS/SCM_Cal/';
	sunpulse_dir = '/Users/argall/Documents/Work/Data/MMS/HK/';
	att_dir      = '/Users/argall/Documents/Work/Data/MMS/Attitude/';
	sc           = 'mms2';
	instr        = 'scm';
	mode         = 'comm';
	optdesc      = 'sc256';
	duration     = 32.0;
	tstart       = '2015-03-17T14:40:00';
	tend         = '2015-03-17T15:00:00';

%------------------------------------%
% Find Files                         %
%------------------------------------%
	% SCM L1A Data File
	[l1a_fname, count, str] = mms_file_search(sc, instr, mode, 'l1a', ...
	                                          'TStart',    tstart, ...
	                                          'TEnd',      tend, ...
	                                          'OptDesc',   optdesc, ...
	                                          'Directory', sc_data_dir);
	assert(count > 0, ['SCM file not found: "' str '".']);
	
	% SCM Cal File
	cal_fname = fullfile(sc_cal_dir, [sc '_' instr sc(4) '_caltab_%Y%M%d%H%m%S_v*.txt']);
	[cal_fname, nFiles] = MrFile_Search( cal_fname, 'VersionRegex', '([0-9])');
	assert(nFiles > 0, ['SCM cal file not found: "' str '".']);

%------------------------------------%
% Calibrated Mag in BCS              %
%------------------------------------%
	[t, ~, b_omb] = mms_sc_bcs(l1a_fname, cal_fname, tstart, tend, duration);
