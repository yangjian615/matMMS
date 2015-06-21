% common_science_constants

c      =  2.99792458e+08;  % speed of light (SI)
e_mass =  9.10938188e-31;  % electron mass, kg (SI)
q      =  1.60217650e-19;  % coulomb (SI)
mV2V   =  1.0e-3;          % mV > V (SI)
nT2T   =  1.0e-9;          % nT > T (SI)
C_V_T  =  mV2V / nT2T;     % Combining constants to save flops; potentially used 100Ks of times

twoPi   = 2.0 * pi;
halfPi  = pi / 2.0;
deg2rad = pi / 180.0;
rad2deg = 180.0 / pi;

Re = 6371; % kilometers

q_over_m = -q / e_mass;
v_1keV_electron  = 18727897.366; % m/s, 18755373. m/s relativistic: difference of 0.147%
v_500eV_electron = 13252328.354; % m/s, 13262052. m/s relativistic: difference of 0.073%
