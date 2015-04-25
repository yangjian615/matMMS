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






% Inputs
edi_data_dir = '/Users/argall/Documents/Work/Data/MMS/EDI/';
dfg_data_dir = '/Users/argall/Documents/Work/Data/MMS/DFG/';
fg_cal_dir   = '/Users/argall/Documents/Work/Data/MMS/FG_Cal/';
att_dir      = '/Users/argall/Documents/Work/Data/MMS/Attitude/';
sunpulse_dir = '/Users/argall/Documents/Work/Data/MMS/HK/';
sc           = 'mms3';
instr        = 'edi';
mode         = 'slow';
level        = 'l1a_efield';
tstart       = '2015-04-18T14:00:00';
tend         = '2015-04-18T14:00:15';

%------------------------------------%
% Get DFG and EDI Data in DMPA       %
%------------------------------------%


% Read EDI data in DMPA
[gd12_dmpa, gd21_dmpa] = mms_edi_gse(sc, 'edi', 'slow', 'l1a_efield', tstart, tend, ...
                                     'SunPulseDir', sunpulse_dir, ...
                                     'DataDir',     edi_data_dir);