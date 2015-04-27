%
% Name
%   mms_edi_tutorial
%
% Purpose
%   Demonstrate how to use the primary EDI functions.
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-14      Written by Matthew Argall
%

get_data = false;
quality  = [];          % Can be vector containing any of [0, 1, 2, 3]

if get_data

	%
	% General flow
	%   mms_edi_gse
	%      Get data in BCS coordinates
	%        mms_edi_bcs -- see below
	%      Rotate from BCS to SMPA
	%        mms_fdoa_read_defatt -- get attitude data and z-MPA axis
	%        mms_fdoa_xbcs2dmpa -- create rotation matrix from BCS to DMPA
	%      Despin
	%        a) mms_fdoa_xdespin -- use attitude data to despin
	%        b) mms_dss_read_sunpulse -- read sunpulse data
	%           mms_dss_xdespin -- use sunpulse data to despin
	%      Rotate to GSE
	%        mms_fdoa_xgei2despun -- Use attitude data to get transformation to GEI
	%        gei2gse -- Hapgood rotation from GEI to GSE
	%
	%   mms_edi_bcs
	%     Read data
	%       mms_edi_efield -- see below
	%     Get gun positions
	%       mms_instr_origins_ocs -- get gun & detector positions
	%     Convert firing angles to firing vectors
	%     Filter by quality
	%     Rotate from EDI to BCS
	%       mms_instr_xxyz2ocs    -- create rotation from instrument frame to BCS
	%
	%   mms_edi_efield
	%     Get data from CDF file.
	%


	% Inputs
	edi_data_dir = '/Users/argall/Documents/Work/Data/MMS/EDI/';
	dfg_data_dir = '/Users/argall/Documents/Work/Data/MMS/DFG/';
	fg_cal_dir   = '/Users/argall/Documents/Work/Data/MMS/FG_Cal/';
	att_dir      = '/Users/argall/Documents/Work/Data/MMS/Attitude/';
	sunpulse_dir = '/Users/argall/Documents/Work/Data/MMS/HK/';
	sc           = 'mms4';
	instr        = 'edi';
	mode         = 'slow';
	level        = 'l1a_efield';
	tstart       = '2015-04-22T17:03:15';
	tend         = '2015-04-22T17:03:28';
	
%------------------------------------%
% Read Data From File                %
%------------------------------------%
	%
	% mms_edi_read_efiled.m
	%   - Read data directly from file.
	%
	
	% Read EDI data directoy from the data file
	%   gd12        Fields are:
	%                 'epoch_gd12'    -  TT2000 Epoch time.
	%                 'polar_gd12'    -  polar firing angle.
	%                 'azimuth_gd12'  -  azimuth firing angles.
	%                 'q_gd12'        -  Quality flag.
	%                 'tof_gd12'      -  Time of flight.
	[gd12, gd21] = mms_edi_read_efield(sc, instr, mode, level, tstart, tend, edi_dir);
	
%------------------------------------%
% Get Data in BCS                    %
%------------------------------------%
	%
	% mms_fg_bcs.m
	%   1. Get gun positions in BCS
	%   2. Call mms_edi_read_efiled.m
	%   3. Select quality
	%   4. Convert firing angles to firing vectors and rotate to BCS
	%   5. Compute location of gun on virtual spacecraft
	%
	
	% Read EDI data in BCS
	%   GD12_BCS    Fields are:
	%                 't_gd12'        -  TT2000 Epoch time for gun 1 and detector 2.
	%                 'gun_gd12_bcs'  -  Gun1 position in BCS.
	%                 'det_gd12_bcs'  -  Detector2 position in BCS.
	%                 'fv_gd12_bcs'   -  Firing vectors from gun1 in BCS.
	%                 'q_gd12'        -  Quality flag.
	%                 'tof_gd12'      -  Time of flight.
	[gd12_bcs, gd21_bcs] = mms_edi_bcs(sc, 'edi', 'slow', 'l1a_efield', tstart, tend, ...
	                                   'DataDir',     edi_data_dir, ...
	                                   'Quality',     quality);
	
%------------------------------------%
% Get Data in DMPA and GSE           %
%------------------------------------%
	%
	% mms_fg_gse.m
	%   1. Call mms_fg_bcs.m
	%   2. Rotate to SMPA
	%   3. Despin firing vectors
	%   4. Despin positions
	%   5. Rotate to GSE
	%
	% Can Return GSE, DMPA, and BCS
	%
	
	
	% Read EDI data in DSL
	%   GD12_DMPA   Fields are:
	%                 'epoch_gd12'     -  TT2000 Epoch time for gun 1 and detector 2.
	%                 'gun_gd12_dmpa'  -  Gun1 position in DMPA.
	%                 'det_gd12_dmpa'  -  Detector2 position in DMPA.
	%                 'gun1_gd12_dmpa' -  Gun1 position on virtual spacecraft in DMPA.
	%                 'fv_gd12_dmpa'   -  Firing vectors from gun1 in DMPA.
	%                 'q_gd12'         -  Quality flag.
	%                 'tof_gd12'       -  Time of flight.
	[gd12_dmpa, gd21_dmpa] = mms_edi_gse(sc, 'edi', 'slow', 'l1a_efield', tstart, tend, ...
	                                    'SunPulseDir', sunpulse_dir, ...
	                                    'DataDir',     edi_data_dir, ...
	                                    'Quality',     quality);
end