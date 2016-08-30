%
% Name
%   mms_fsm_test_cals
%
% Purpose
%   Plot calibration data for the FSM data product..
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2016-06-05      Written by Matthew Argall
%
%***************************************************************************

% Data file to plot
file = '/nfs/fsm/temp/mms1_fsm_srvy_l2plus_cal_20151108150948_v0.0.0.cdf';

%------------------------------------%
% Get Data                           %
%------------------------------------%

% Parse file name to get variable name prefix and suffix
[sc, instr, mode, level] = mms_dissect_filename(file);
prefix = [sc '_' instr '_'];
suffix = ['_' mode '_' level];

% Create Variable Names
fgm_amp_vname   = [prefix 'amp'   '_omb_fgm' suffix];
fgm_phase_vname = [prefix 'phase' '_omb_fgm' suffix];
fgm_psd_vname   = [prefix 'psd'   '_omb_fgm' suffix];
scm_amp_vname   = [prefix 'amp'   '_omb_scm' suffix];
scm_phase_vname = [prefix 'phase' '_omb_scm' suffix];
scm_psd_vname   = [prefix 'psd'   '_omb_scm' suffix];
flag_vname      = [prefix 'flag'             suffix];

% Read data
[amp_fgm, t, f] = MrCDF_Read(file, fgm_amp_vname);
phase_fgm       = MrCDF_Read(file, fgm_phase_vname);
psd_fgm         = MrCDF_Read(file, fgm_psd_vname);
amp_scm         = MrCDF_Read(file, scm_amp_vname);
phase_scm       = MrCDF_Read(file, scm_phase_vname);
psd_scm         = MrCDF_Read(file, scm_psd_vname);


keyboard

%------------------------------------%
% Plot Data                          %
%------------------------------------%
nTime = length(t);
nFreq = length(f);
x     = repmat( MrCDF_epoch2ssm(t)', 1, nFreq);
y     = repmat(f',  nTime, 1);

surf( double(x), double(y), double(amp_fgm(:,:,ii)) );




for ii = 1 : 3
	fig = figure();
	
	subplot(3,2,1);
	surf( double(x), double(y), double(amp_fgm(:,:,ii)) );

	subplot(3,2,2)
	contour( t, f, phase_fgm(:,:,ii) );

	subplot(3,2,3)
	contour( t, f, pwr_fgm(:,:,ii) );
	
	subplot(3,2,4);
	contour( t, f, amp_scm(:,:,ii) );

	subplot(3,2,5)
	contour( t, f, phase_scm(:,:,ii) );

	subplot(3,2,6)
	contour( t, f, pwr_scm(:,:,ii) );
end