%
% Name
%   mms_adsm_compare_l1b
%
% Purpose
%   Compare the AFG, DFG, and SCM L1B data.
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-07-31      Written by Matthew Argall
%

%------------------------------------%
% Inputs                             %
%------------------------------------%
sc      = 'mms2';
instr   = 'dfg';
mode    = 'f128';
optdesc = '';
tstart  = '2015-06-30T14:40:00';
tend    = '2015-06-30T15:00:00';

%------------------------------------%
% Find Files                         %
%------------------------------------%
% DFG L1A Data File
[fname, count, str] = mms_file_search(sc, instr, mode, 'l1a',  ...
                                          'TStart',    tstart, ...
                                          'TEnd',      tend);
assert(count > 0, ['Mag L1A file not found: "' str '".']);

% DSS L1B Data File
[dss_fname, count, str] = mms_file_search(sc, 'fields', 'hk', 'l1b', ...
                                          'OptDesc',   '101',        ...
                                          'TStart',    tstart,       ...
                                          'TEnd',      tend,         ...
                                          'SDCroot',   '/nfs/hk/');
assert(count > 0, ['DSS L1A file not found: "' str '".']);

%------------------------------------%
% Read Data                          %
%------------------------------------%
% Magnetometer data
switch instr
	case 'afg'
		mag_l1a = mms_fg_read_l1a( fname, tstart, tend );
	case 'dfg'
		mag_l1a = mms_fg_read_l1a( fname, tstart, tend );
	case 'dfg'
		mag_l1a = mms_sc_read_l1a( fname, tstart, tend );
	otherwise
		error( 'INSTR must be afg, dfg, or scm.' )
end

% Sunpulse data
sunpulse = mms_dss_read_sunpulse( dss_fname, tstart, tend, 'UniquePulse', true );

%------------------------------------%
% Find a Spin                        %
%------------------------------------%

% Find the first interval
idss  = find( sunpulse.SunPulse > mag_l1a.tt2000(1), 1, 'first' );
ispin = find( mag_l1a.tt2000 >= sunpulse.SunPulse(idss) & ...
              mag_l1a.tt2000 <= sunpulse.SunPulse(idss+1) );

% Extract one spin
t = mag_l1a.tt2000(ispin);
b = mag_l1a.b_123(:, ispin);

% Convert to seconds
t = double( t - t(1) ) * 1e-9;

%------------------------------------%
% Model Fit                          %
%------------------------------------%

% FitType object
myFit = fittype( 'a + b*sin(w*x) + c*cos(w*x)',       ...
                 'Coefficients', {'a', 'b', 'c', 'w'} );

% Model fit
myFit = fit(t', double(b(1,:)'), myFit, ...
            'StartPoint', [ 0, max(b(1,:)), max(b(1,:)), 2*pi/(t(end)-t(1)) ] );

keyboard;
plot(t, b(1,:));
hold on
plot(myFit);
hold off

%------------------------------------%
% Plot                               %
%------------------------------------%
