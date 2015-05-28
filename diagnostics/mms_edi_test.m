%
% Name
%   mms_edi_test
%
% Purpose
%   Plot all EDI beams in the given time interval in BPP using the average
%   magnetic field for that interval. All beams are projected into the same
%   BPP.
%
% See Also
%   mms_edi_test
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-14      Written by Matthew Argall
%

get_data = false;
quality  = 3;

% Inputs
if get_data
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
	tend         = '2015-04-22T17:03:30';
	
	%------------------------------------%
	% Get DFG and EDI Data in DSL        %
	%------------------------------------%
	% Read DFG data in DSL
	[t, b_bcs] = mms_fg_bcs(sc, 'dfg', 'f128', tstart, tend, ...
	                        'CalDir',  fg_cal_dir, ...
	                        'DataDir', dfg_data_dir);

	% EDI data
	[gd12, gd21] = mms_edi_read_efield(sc, instr, mode, level, tstart, tend, edi_data_dir);
end

% Indices for times
i_gd12 = 1:length(gd12.azimuth_gd12);
i_gd21 = 1:length(gd21.azimuth_gd21);

%------------------------------------%
% Firing Angles                      %
%------------------------------------%

% Create the figure
fig_fa = figure();

% Firing Angles from GD12
subplot(2,1,1)
plot(i_gd12, gd12.azimuth_gd12, i_gd12, gd12.polar_gd12);
title('Firing Angles GD12')
ylabel('Degrees')
ylim([-180,180])
legend('Azimuth', 'Polar')

% Firing Angles from GD21
subplot(2,1,2)
plot(i_gd21, gd21.azimuth_gd21, i_gd21, gd21.polar_gd21);
title('Firing Angles GD21')
ylabel('Degrees')
ylim([-180,180])
legend('Azimuth', 'Polar')

%------------------------------------%
% Firing Vectors                     %
%------------------------------------%

fig_fv = figure();

% Cartesian coordinates
fv_gd12      = zeros(3, length(gd12.epoch_gd12));
fv_gd12(1,:) = sind( gd12.polar_gd12 ) .* cosd( gd12.azimuth_gd12 );
fv_gd12(2,:) = sind( gd12.polar_gd12 ) .* sind( gd12.azimuth_gd12 );
fv_gd12(3,:) = cosd( gd12.polar_gd12 );

fv_gd21      = zeros(3, length(gd21.epoch_gd21));
fv_gd21(1,:) = sind( gd21.polar_gd21 ) .* cosd( gd21.azimuth_gd21 );
fv_gd21(2,:) = sind( gd21.polar_gd21 ) .* sind( gd21.azimuth_gd21 );
fv_gd21(3,:) = cosd( gd21.polar_gd21 );

% Plot vectors from gun1
subplot(2,1,1)
plot(i_gd12, fv_gd12);
title('Firing Vectors BCS');
legend('X', 'Y', 'Z');

subplot(2,1,2);
plot(i_gd21, fv_gd21)
title('Firing Vectors BCS')
legend('X', 'Y', 'Z')


%------------------------------------%
% |B| and ToF                        %
%------------------------------------%

%
% Plot B from ToF
%   2 pi / ToF = q B / m
%
%   B = (2 pi m) / (q ToF)
%
b_mag      = mrvector_magnitude(t);
b_tof_gd12 = 1e9 * ( 2 * pi * 9.11e-31 ) ./ ( 1.602e-19 * double(gd12.tof_gd12) * 1e-6 );
b_tof_gd21 = 1e9 * ( 2 * pi * 9.11e-31 ) ./ ( 1.602e-19 * double(gd21.tof_gd21) * 1e-6 );

fig_tof = figure();

% Plot |B|
p1 = plot(b_bcs, b_mag);
hold on
title('|B| from DFG and ToF')

% Plot B-ToF
p2 = scatter(gd12.epoch_gd12, b_tof_gd12, [], 'blue');
p3 = scatter(gd21.epoch_gd21, b_tof_gd21, [], 'red');
legend([p1 p2 p3], 'B_{DFG}', 'B_{GD12}', 'B_{GD21}')
hold off

