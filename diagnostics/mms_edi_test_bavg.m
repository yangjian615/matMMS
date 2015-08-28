%
% Name
%   mms_edi_test_bavg
%
% Purpose
%   Create a MATLAB save file of inputs needed for Bestarg.
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-08-27      Written by Matthew Argall
%
function [] = mms_edi_test_bavg()

sc     = 'mms2';
tstart = '2015-05-09T00:00:00';
tend   = '2015-05-09T24:00:00';

% Defaults
beam_quality    = 3;
dt              = 5.0;

%------------------------------------%
% Find Files                         %
%------------------------------------%

% FG L1B Data File
instr   = 'dfg';
mode    = 'srvy';
level   = 'l1b';
optdesc = '';
[fg_l1b_file, count, str] = mms_file_search(sc, instr, mode, level, ...
                                            'TStart',    tstart, ...
                                            'TEnd',      tend, ...
                                            'OptDesc',   optdesc);
assert(count > 0, ['DFG L1B file not found: "' str '".']);

% EDI Slow L1A E-Field Data File
%    - Find last, so file descriptors are saved.
instr   = 'edi';
mode    = 'slow';
level   = 'l1a';
optdesc = 'efield';
[edi_slow_file, slw_cnt, str] = mms_file_search(sc, instr, mode, level, ...
                                                'TStart',    tstart, ...
                                                'TEnd',      tend, ...
                                                'OptDesc',   optdesc);

% EDI Fast L1A E-Field Data File
%    - Find last, so file descriptors are saved.
instr   = 'edi';
mode    = 'fast';
level   = 'l1a';
optdesc = 'efield';
[edi_fast_file, fst_cnt, str] = mms_file_search(sc, instr, mode, level, ...
                                                'TStart',    tstart, ...
                                                'TEnd',      tend, ...
                                                'OptDesc',   optdesc);

%------------------------------------%
% Read Data                          %
%------------------------------------%

% EDI Fast data
edi_slow = [];
if slw_cnt > 0
	edi_slow = mms_edi_l1b_efield_create( edi_slow_file, tstart, tend, ...
	                                        'CS_BCS',   true,          ...
	                                        'Quality',  beam_quality  );
end

% EDI Slow data
edi_fast = [];
if fst_cnt > 0
	%
	% Occasionally, a file will be found with zero records. This is a
	% problem with L1A processing which produces files when no data
	% is available.
	%
	try
		edi_fast = mms_edi_l1b_efield_create( edi_slow_file, tstart, tend, ...
		                                        'CS_BCS',   true,          ...
		                                        'Quality',  beam_quality  );
	% Not the error and continue processing slow survey data.
	catch ME
		mrfprintf('logerr', ME);
	end
end

% Check after we check for empty fast survey file.
assert( slw_cnt + fst_cnt > 0, 'Unable to find fast or slow survey EDI files.' );

% Combine slow and fast survey data
edi = mms_edi_srvy_combine( edi_slow, edi_fast );
clear edi_slow edi_fast

% FGM data
fg_l1b = mms_fg_read_l1b(fg_l1b_file, tstart, tend);

%------------------------------------%
% Compuate Average B                 %
%------------------------------------%

% Time range that we have data
if isempty(edi.tt2000_gd12) && isempty(edi.tt2000_gd12)
	error( 'No EDI E-field data available.' );
elseif isempty(edi.tt2000_gd12)
	t0 = edi.tt2000_gd21(1);
	t1 = edi.tt2000_gd21(end);
elseif isempty(edi.tt2000_gd12)
	t0 = edi.tt2000_gd12(1);
	t1 = edi.tt2000_gd12(end);
else
	t0 = min( [edi.tt2000_gd12(1)   edi.tt2000_gd21(1)  ] );
	t1 = max( [edi.tt2000_gd12(end) edi.tt2000_gd21(end)] );
end

% Breakdown into time vectors
tvec = MrCDF_Epoch_Breakdown( [t0, t1] );

% Round down to the nearest DT seconds and recompute
tvec(:,6)     = tvec(:,6) - mod(tvec(:,6), dt);
tvec(:,7:end) = 0.0;
tedge         = MrCDF_Epoch_Compute(tvec);
t0            = tedge(1);
t1            = tedge(2) + int64(dt * 1d9);

% Find FGM data within this time interval
%   - Prune out |B|
ifg      = find( fg_l1b.tt2000 >= t0 & fg_l1b.tt2000 <= t1 );
t_fg     = fg_l1b.tt2000(ifg);
b_fg_bcs = fg_l1b.b_bcs(1:3, ifg);

% Compute the averaged magnetic field
b_avg_bcs = mms_edi_bavg(t_fg, b_fg_bcs, edi.tt2000_gd12, edi.tt2000_gd21, dt);

%------------------------------------%
% Plot Results                       %
%------------------------------------%
fig = figure()

subplot(4,1,1)
plot( t_fg, mrvector_magnitude(b_fg_bcs), ...
      edi.tt2000_gd12, mrvector_magnitude(b_avg_bcs.b_gd12), ...
      edi.tt2000_gd21, mrvector_magnitude(b_avg_bcs.b_gd21), ...
      b_avg_bcs.t_avg, mrvector_magnitude(b_avg_bcs.b_avg) );

subplot(4,1,2)
plot( t_fg, b_fg_bcs(1,:), ...
      edi.tt2000_gd12, b_avg_bcs.b_gd12(1,:), ...
      edi.tt2000_gd21, b_avg_bcs.b_gd21(1,:), ...
      b_avg_bcs.t_avg, b_avg_bcs.b_avg(1,:) );

subplot(4,1,3)
plot( t_fg, b_fg_bcs(2,:), ...
      edi.tt2000_gd12, b_avg_bcs.b_gd12(2,:), ...
      edi.tt2000_gd21, b_avg_bcs.b_gd21(2,:), ...
      b_avg_bcs.t_avg, b_avg_bcs.b_avg(2,:) );

subplot(4,1,4)
plot( t_fg, b_fg_bcs(3,:), ...
      edi.tt2000_gd12, b_avg_bcs.b_gd12(3,:), ...
      edi.tt2000_gd21, b_avg_bcs.b_gd21(3,:), ...
      b_avg_bcs.t_avg, b_avg_bcs.b_avg(3,:) );
