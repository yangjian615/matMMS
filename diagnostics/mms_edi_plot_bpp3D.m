%
% Name
%   mms_edi_plot_bpp3D
%
% Purpose
%   Plot all EDI beams in the given time interval in BPP using the average
%   magnetic field for that interval. All beams are projected into the same
%   BPP. Results are shown in 3D.
%
% See Also
%   mms_edi_view
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-27      Written by Matthew Argall
%

get_data = true;
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
	[t, ~, b_dmpa] = mms_fg_gse(sc, 'dfg', 'f128', tstart, tend, ...
	                           'CalDir',      fg_cal_dir, ...
	                           'SunPulseDir', sunpulse_dir, ...
	                           'DataDir',     dfg_data_dir);
	
	% Read EDI data in DSL
	[gd12_dmpa, gd21_dmpa] = mms_edi_gse(sc, 'edi', 'slow', 'l1a_efield', tstart, tend, ...
	                                    'SunPulseDir', sunpulse_dir, ...
	                                    'DataDir',     edi_data_dir, ...
	                                    'Quality',     quality);
end

%------------------------------------%
% View in Bavg Plane                 %
%------------------------------------%
% Average magnetic field
b_avg = mean(b_dmpa, 2);
b_std = std(b_dmpa, 0, 2);

% Create a single transformation matrix for all beams
xyz2bpp_bavg = mms_edi_xxyz2bpp(b_avg);

% Rotate firing vectors and gun positions into bpp
fv_gd12_bpp  = mrvector_rotate(xyz2bpp_bavg, gd12_dmpa.fv_gd12_dmpa);
pos_gun1_bpp = mrvector_rotate(xyz2bpp_bavg, gd12_dmpa.gun1_dmpa);

fv_gd21_bpp  = mrvector_rotate(xyz2bpp_bavg, gd21_dmpa.fv_gd21_dmpa);
pos_gun2_bpp = mrvector_rotate(xyz2bpp_bavg, gd21_dmpa.gun2_dmpa);

%------------------------------------%
% Spacecraft Coordinates             %
%------------------------------------%

% Create the spacecraft
sc_radius   = mms_instr_origins_instr('EDI1_GUN', 'EDI2_DETECTOR');
sc_radius   = sqrt( sum( sc_radius.^2 ) );
sc_npts     = 100;
sc_sph      = zeros(3, sc_npts);
sc_sph(3,:) = sc_radius;
sc_sph(1,:) = 2 * pi * (1:1:sc_npts) / sc_npts;
[sc_x, sc_y, sc_z]  = sph2cart(sc_sph(1,:), sc_sph(2,:), sc_sph(3,:));
sc_xyz = [sc_x; sc_y; sc_z];

% Rotate the spacecraft into BPP
sc_bpp = mrvector_rotate( xyz2bpp_bavg, sc_xyz );

% Clear temporary variables
clear sc_radius sc_nps sc_sph sc_x sc_y sc_z

%------------------------------------%
% Firing Vectors                     %
%------------------------------------%

range = [-10, 10];

% Remove NaNs

iGood        = find( ~isnan(fv_gd12_bpp(1,:)) );
nGood_gd12   = length(iGood);
fv_gd12_bpp  = fv_gd12_bpp(:,iGood);
pos_gun1_bpp = pos_gun1_bpp(:,iGood);

iGood        = find( ~isnan(fv_gd21_bpp(1,:)) );
nGood_gd21   = length(iGood);
fv_gd21_bpp  = fv_gd21_bpp(:,iGood);
pos_gun2_bpp = pos_gun2_bpp(:,iGood);

% Skip this round
if nGood_gd12 == 0 && nGood_gd21 == 0
	error('No beams found. Cannot plot.');
end

%
% Knowing two vectors that point from the origin to two points in space
%   r0 = (x0, y0, z0)
%   r1 = (x1, y1, z1)
%
% A vector, v, parallel to the vector a, which connects r0 to r1, can be formed
%   v = r1 - r0
%     = (x1 - x0, y1 - y0, z1 - z0)
%     = (a, b, c)
%
% A vector from the origin to any other point on the line  can be formed by adding some
% multiplicative factor, t, of v to r0
%   r = r0 + vt
%
% It follows that
%   x = x0 + at
%   y = y0 + bt
%   z = z0 + ct
%
% Or
%   t = (x - x0) / a
%     = (y - y0) / b
%     = (z - z0) / c
%
% Using the x-component of the firing vectors and gun positions, we can determine t
% for the range of the plot
%   t(1) = ( pos_gun(1) - range(1) ) / fv(1)
%   t(1) = ( pos_gun(2) - range(2) ) / fv(2)
%
% Then, we can find the (x,y,z) coordinates that define the end points of our line.
%   x = gun_pos(1) + fv(1) * t
%   y = gun_pos(2) + fv(2) * t
%   z = gun_pos(3) + fv(3) * t
%

% GD12
%   - find (x0,x1), (y0,y1), and (z0,z1) for each beam.
t_gd12      = zeros(2, nGood_gd12);
x_gd12      = zeros(2, nGood_gd12);
y_gd12      = zeros(2, nGood_gd12);
z_gd12      = zeros(2, nGood_gd12);
t_gd12(1,:) = ( pos_gun1_bpp(1,:) - range(1) ) ./ fv_gd12_bpp(1,:);
t_gd12(2,:) = ( pos_gun1_bpp(1,:) - range(2) ) ./ fv_gd12_bpp(1,:);
x_gd12(1,:) = pos_gun1_bpp(1,:) + fv_gd12_bpp(1,:) .* t_gd12(1,:);
x_gd12(2,:) = pos_gun1_bpp(1,:) + fv_gd12_bpp(1,:) .* t_gd12(2,:);
y_gd12(1,:) = pos_gun1_bpp(2,:) + fv_gd12_bpp(2,:) .* t_gd12(1,:);
y_gd12(2,:) = pos_gun1_bpp(2,:) + fv_gd12_bpp(2,:) .* t_gd12(2,:);
z_gd12(1,:) = pos_gun1_bpp(3,:) + fv_gd12_bpp(3,:) .* t_gd12(1,:);
z_gd12(2,:) = pos_gun1_bpp(3,:) + fv_gd12_bpp(3,:) .* t_gd12(2,:);

% GD21
t_gd21      = zeros(2, nGood_gd21);
x_gd21      = zeros(2, nGood_gd21);
y_gd21      = zeros(2, nGood_gd21);
z_gd21      = zeros(2, nGood_gd21);
t_gd21(1,:) = ( pos_gun2_bpp(1,:) - range(1) ) ./ fv_gd21_bpp(1,:);
t_gd21(2,:) = ( pos_gun2_bpp(1,:) - range(2) ) ./ fv_gd21_bpp(1,:);
x_gd21(1,:) = pos_gun2_bpp(1,:) + fv_gd21_bpp(1,:) .* t_gd21(1,:);
x_gd21(2,:) = pos_gun2_bpp(1,:) + fv_gd21_bpp(1,:) .* t_gd21(2,:);
y_gd21(1,:) = pos_gun2_bpp(2,:) + fv_gd21_bpp(2,:) .* t_gd21(1,:);
y_gd21(2,:) = pos_gun2_bpp(2,:) + fv_gd21_bpp(2,:) .* t_gd21(2,:);
z_gd21(1,:) = pos_gun2_bpp(3,:) + fv_gd21_bpp(3,:) .* t_gd21(1,:);
z_gd21(2,:) = pos_gun2_bpp(3,:) + fv_gd21_bpp(3,:) .* t_gd21(2,:);

% Clear
clear t_gd12 t_gd21

%------------------------------------%
% Plot Spacecraft Outline            %
%------------------------------------%

% Title information
ttl_tstart = MrCDF_Epoch_Encode(t(1));
ttl_tend   = MrCDF_Epoch_Encode(t(end));
ttl_date   = ttl_tstart{1}(1:10);
ttl_tstart = ttl_tstart{1}(12:19);
ttl_tend   = ttl_tend{1}(12:19);
ttl = sprintf( ['EDI eField Mode\n' ...
                'Time Interval: %s  %s - %s\n' ...
                'nBeams: GD12 = %d, GD21 = %d\n' ...
                'Bavg = [%0.2f, %0.2f, %0.2f] +/- [%0.2f, %0.2f, %0.2f]'], ...
              ttl_date, ttl_tstart, ttl_tend, nGood_gd12, nGood_gd21, ...
              b_avg, b_std );

% Append quality to title
if ~isempty(quality)
	ttl = sprintf('%s\nQuality = %d', ttl, quality);
end

% S/C Outline
plot3( sc_bpp(1,:), sc_bpp(2,:), sc_bpp(3,:) );
hold on
grid on;
title(ttl);
xlabel('x (m)');
ylabel('y (m)');

%---------------------------------------%
% Plot Gun Positions & Firing Vectors   %
%---------------------------------------%
% Create a scatter plot of gun positions
s_gd12 = scatter3(pos_gun1_bpp(1,:), pos_gun1_bpp(2,:), pos_gun1_bpp(3,:), [], 'blue');
s_gd21 = scatter3(pos_gun2_bpp(1,:), pos_gun2_bpp(2,:), pos_gun2_bpp(3,:), [], 'red');

% Create lines
l_gd12 = [];
l_gd21 = [];
for ii = 1 : nGood_gd12
	l_gd12 = line(x_gd12(:,ii), y_gd12(:,ii), z_gd12(:,ii), 'Color', 'blue');
end
for ii = 1 : nGood_gd21
	l_gd21 = line(x_gd21(:,ii), y_gd21(:,ii), z_gd21(:,ii), 'Color', 'red');
end
legend([l_gd12, l_gd21], 'GD12', 'GD21');

xlim(range)
ylim(range)
zlim(range)
hold off
