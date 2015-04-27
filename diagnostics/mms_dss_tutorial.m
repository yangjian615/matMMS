%
% Name
%   mms_dss_tutorial
%
% Purpose
%   Demonstrate how to use DSS functions.
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-27      Written by Matthew Argall
%


%
% General flow
%   Read sunpulse data
%     mms_dss_read_sunpulse
%
%   Create transformation matrix
%     mms_dss_xdespin -- build transformation matrix
%       mms_dss_sunpulse2phase -- convert period & time to spin phase
%

% Inputs
sunpulse_dir = '/Users/argall/Documents/Work/Data/MMS/HK/';
sc           = 'mms4';
tstart       = '2015-04-22T17:03:15';
tend         = '2015-04-22T17:03:30';

%------------------------------------%
% Read Data From File                %
%------------------------------------%

% Read sun pulse data
%   SUNPULSE  Fields are:
%               'Epoch'    - Packet times
%               'SunPulse' - Sun pulse times
%               'Period'   - Period (micro-sec) of revolution Only
%                              returned when FLAG=0 and only on the
%                              second and subsequent received sun
%                              pulses from the s/c.
%               'Flag'     - Status flag
%                              0: s/c sun pulse
%                              1: s/c pseudo sun pulse
%                              2: s/c CIDP generated speudo sun pulse
sunpulse = mms_dss_read_sunpulse(sc, tstart, tend, ...
                                 'UniquePulse', true, ...
                                 'Directory', sunpulse_dir);

%------------------------------------%
% Create Tranformation Matrix        %
%------------------------------------%

% Supply EPOCH_TIME -- the time at which you have data
EPOCH_TIME = zeros(0, 0, 'int64');
smpa2dmpa  = mms_dss_xdespin( sunpulse, EPOCH_TIME );

