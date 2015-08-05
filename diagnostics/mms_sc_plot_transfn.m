%
% Name
%   mms_sc_plot_transfn
%
% Purpose
%   Plot the SCM transfer function gain and phase
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
sc     = 'mms1';
instr  = 'scm';

% Directories
sc_cal_dir = fullfile('/home', 'argall', 'data', 'mms', 'scm_cal');

%------------------------------------%
% Find Files                         %
%------------------------------------%

if strcmp(sc, 'mms1')
	fm = 'scm1';
elseif strcmp(sc, 'mms2')
	fm = 'scm2';
elseif strcmp(sc, 'mms3')
	fm = 'scm4';
elseif strcmp(sc, 'mms4')
	fm = 'scm3';
else
	error( ['Unknown spacecraft ID "' sc '".'] );
end

% SCM Cal File
ftest = fullfile(sc_cal_dir, [sc '_' fm '_caltab_%Y%M%d%H%m%S_v*.txt']);
[cal_fname, nFiles] = MrFile_Search( ftest, 'VersionRegex', '([0-9])');
assert(nFiles > 0, ['SCM cal file not found: "' ftest '".']);

%------------------------------------%
% Read Data                          %
%------------------------------------%

% Read the cal file
[transfr_fn, f] = mms_sc_read_caltab(cal_fname);

% Gain
gain = 20.0 * log10( abs( transfr_fn ) );

% Phase
phase = atan( imag(transfr_fn) ./ real(transfr_fn) ) * (180 / pi);

% Quadrants
Q2 = imag(transfr_fn) > 0 & real(transfr_fn) < 0;
Q3 = imag(transfr_fn) < 0 & real(transfr_fn) < 0;

% Shift phase
phase( Q2 ) = phase( Q2 ) + 180.0;
phase( Q3 ) = phase( Q3 ) - 180.0;

%------------------------------------%
% Plot Gain and Phase                %
%------------------------------------%
f_smpa = figure();

% Gain
subplot(2,1,1)
semilogx( f, gain );
title([ 'Transfer Function Gain'] );
xlabel( 'Frequency (Hz)' );
ylabel( 'Gain (dB)' );
legend('X', 'Y', 'Z');

% Phase
subplot(2,1,2)
semilogx( f, phase );
title([ 'Transfer Function Phase'] );
xlabel( 'Frequency (Hz)' );
ylabel( 'Phase (deg)' );