%
% Name
%   mms_edi_ev_compare
%
% Purpose
%   Compare EDI E-field and VxB drift velocity with that
%   of EDP/DFG.
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-04-14      Written by Matthew Argall
%

get_data      = true;
tf_EdotB_zero = true;

if get_data
%------------------------------------%
% Inputs                             %
%------------------------------------%
	sc         = 'mms2';
	tstart     = '2015-07-22T09:00:00';
	tend       = '2015-07-22T17:00:00';
	edi_dir    = '/nfs/edi/';
	att_dir    = fullfile('/nfs', 'ancillary', sc, 'defatt');
	eph_dir    = fullfile('/nfs', 'ancillary', sc, 'defeph');
	
%------------------------------------%
% Find Files                         %
%------------------------------------%
	% EDI QL File
	[edi_fname, count, str] = mms_file_search(sc, 'edi', 'srvy', 'ql', ...
	                                          'Directory', edi_dir,    ...
	                                          'OptDesc',   'efield',   ...
	                                          'TStart',    tstart,     ...
	                                          'TEnd',      tend);
	assert(count > 0, ['EDI file not found: "' str '".']);
	
	% DFG QL File
	[fg_fname, count, str] = mms_file_search(sc, 'dfg', 'srvy', 'ql', ...
	                                         'TStart',  tstart, ...
	                                         'TEnd',    tend);
	assert(count > 0, ['DFG file not found: "' str '".']);
	
	% EDP Fast QL File
	[edp_ffast, fcnt, str] = mms_file_search(sc, 'edp', 'fast', 'ql', ...
	                                          'OptDesc', 'dce2d',          ...
	                                          'TimeOrder', '%Y%M%d%H%m%S', ...
	                                          'TStart',  tstart,           ...
	                                          'TEnd',    tend);
	
	% EDP slow QL File
	[edp_fslow, scnt, str] = mms_file_search(sc, 'edp', 'slow', 'ql', ...
	                                          'OptDesc', 'dce2d',          ...
	                                          'TimeOrder', '%Y%M%d%H%m%S', ...
	                                          'TStart',  tstart,           ...
	                                          'TEnd',    tend);
	assert(fcnt + scnt > 0, ['No EDP fast or slow survey files found.']);
	
	% Attitude files
	att_ftest = fullfile( att_dir, [upper(sc) '_DEFATT_%Y%D_%Y%D.V*'] );
	[defatt_fname, count] = MrFile_Search(att_ftest, ...
	                                      'Closest',      true, ...
	                                      'TimeOrder',    '%Y%D', ...
	                                      'TStart',       tstart, ...
	                                      'TEnd',         tend, ...
	                                      'VersionRegex', 'V[0-9]{2}');
	assert(count > 0, ['No definitive attitude file found: "' att_ftest '".']);
	
	% Ephemeris files
	eph_ftest = fullfile( eph_dir, [upper(sc) '_DEFEPH_%Y%D_%Y%D.V*'] );
	[defeph_fname, count] = MrFile_Search(eph_ftest, ...
	                                      'Closest',      true, ...
	                                      'TimeOrder',    '%Y%D', ...
	                                      'TStart',       tstart, ...
	                                      'TEnd',         tend, ...
	                                      'VersionRegex', 'V[0-9]{2}');
	assert(count > 0, ['No definitive ephemeris file found: "' att_ftest '".']);

%------------------------------------%
% Read Data                          %
%------------------------------------%
	% Get data
	edi_ql   = mms_edi_read_ql_efield(edi_fname,  tstart, tend);
	fg_ql    = mms_fg_read_ql(fg_fname,           tstart, tend);
	attitude = mms_fdoa_read_defatt(defatt_fname, tstart, tend);
	ephem    = mms_fdoa_read_defeph(defeph_fname, tstart, tend);
	
	% Read EDP
	if fcnt > 0
		edp_fast = mms_edp_read_ql(edp_ffast, tstart, tend);
	end
	if scnt > 0
		edp_slow = mms_edp_read_ql(edp_fslow, tstart, tend);
	end
	
	% Combine fast and slow
	if fcnt > 0 && scnt > 0
		% Collect and sort times
		edp_ql.tt2000          = [edp_fast.tt2000 edp_slow.tt2000];
		[edp_ql.tt2000, isort] = sort( edp_ql.tt2000 );
		
		% Collect E-field data
		edp_ql.E_dsl = [edp_fast.E_dsl edp_slow.E_dsl];
		edp_ql.E_dsl = edp_ql.E_dsl(:, isort);
	
		% Clear data
		clear edp_fast edp_slow isort
	
	% FAST
	elseif fcnt > 0
		edp_ql = edp_fast;
		clear edp_fast
		
	% SLOW
	else
		edp_ql = edp_slow;
		clear edp_slow
	end
	
	% Reference time
	t0 = min( [edi_ql.tt2000(1), fg_ql.tt2000(1), edp_ql.tt2000(1)] );
	t_edi = MrCDF_epoch2sse( edi_ql.tt2000,   t0 );
	t_edp = MrCDF_epoch2sse( edp_ql.tt2000,   t0 );
	t_fg  = MrCDF_epoch2sse( fg_ql.tt2000,    t0 );
	t_att = MrCDF_epoch2sse( attitude.tt2000, t0 );
	t_eph = MrCDF_epoch2sse( ephem.tt2000,    t0 );
	
	% Remove fill values
	ifill = find( edp_ql.E_dsl(1,:) == single(-1e31) );
	if ~isempty(ifill)
		edp_ql.E_dsl(:, ifill) = NaN;
	end
	ifill = find( edi_ql.E_dmpa(1,:) == single(-1e31) );
	if ~isempty(ifill)
		edi_ql.E_dmpa(:, ifill) = NaN;
		edi_ql.v_ExB_dmpa(:, ifill) = NaN;
	end
	

%------------------------------------%
% Compute ExB Drift Velocity         %
%------------------------------------%
	% Interpolate EDP onto DFG
	E = zeros(3, length(t_fg));
	for ii = 1 : 3
		E(ii,:) = interp1( t_edp, edp_ql.E_dsl(ii,:), t_fg );
	end

	% Compute Ez using E dot B = 0?
	if tf_EdotB_zero
		E(3,:) = -( E(1,:) .* fg_ql.b_dmpa(1,:) + E(2,:) .* fg_ql.b_dmpa(2,:) ) ./ fg_ql.b_dmpa(3,:);
	end

	% Compute the cross product
	%   - v = ExB/|B|^2
	%   - Scale: mV/m * nT / (nT)^2  -->  V/m / T * 1e-6  -->  km/s * 1e3
	v_ExB = 1e-3 * mrvector_cross( E, fg_ql.b_dmpa(1:3,:) );
	v_ExB = v_ExB ./ repmat( fg_ql.b_dmpa(4,:).^2, 3, 1 );


%------------------------------------%
% Compute VxB s/c E-Field            %
%------------------------------------%
	% Interpolate ephemeris position and velocity to DFG times
	[r_gei, v_gei] = mms_fdoa_interp_ephem(ephem.tt2000, ephem.Position, fg_ql.tt2000, ephem.Velocity);
	
	% Rotate from GEI to DMPA
	gei2dmpa = mms_fdoa_xgei2despun(attitude, fg_ql.tt2000, 'L');
	v_dmpa   = mrvector_rotate(gei2dmpa, v_gei);
	
	% Compute the cross product
	%   - E = -1.0 * (V x B)
	%   - Convert: km/s * nT  -->  m/s * T * 1e-6  -->  V / m * 1e-6  -->  mV/m * 1e-3
	E_VxB = -1e-3 * mrvector_cross( v_dmpa, fg_ql.b_dmpa(1:3,:) );


%------------------------------------%
% Co-Rotation E-Field                %
%------------------------------------%
	% Rotation vector of Earth in GEI
	w_gei = [0, 0, 2.0*pi / (24.0*60.0*60.0)];
	
	% Rotation vector of Earth in DMPA
	r_dmpa = mrvector_rotate( gei2dmpa, r_gei );
	w_dmpa = mrvector_rotate( gei2dmpa, w_gei );
	
	% Co-rotation velocity at the position of the spacecraft
	v_CoRot = mrvector_cross( w_dmpa, r_dmpa );
	
	% Co-rotation electric field
	%   - E = -(v x B)
	%   - Convert: km/s * nT  -->  m/s * T * 1e-6  -->  V / m * 1e-6  -->  mV/m * 1e-3
	E_CoRot = -1e-3 .* mrvector_cross( v_CoRot, fg_ql.b_dmpa(1:3,:) );
end

%------------------------------------%
% Plot B, E, E_VxB                   %
%------------------------------------%

%
% SMPA
%
fig_efield = figure();

% Magnetic Field
subplot(4,1,1)
plot( t_fg, fg_ql.b_dmpa );
title([ 'Magnetic and Electric Fields'] );
ylabel( {'B', '(nT)'} );
legend( {'B_{x}', 'B_{y}', 'B_{z}', '|B|'} );

% Ex
subplot(4,1,2)
plot( t_edp, edp_ql.E_dsl(1,:),  ...
      t_edi, edi_ql.E_dmpa(1,:), ...
      t_fg,  E_VxB(1,:),         ...
      t_fg,  E_CoRot(1,:) );
ylabel( {'Ex', '(mV/m)'} );
legend( {'E_{EDP}', 'E_{EDI}', 'E_{S/C}', 'E_{CoRot}' } );

% Ey
subplot(4,1,3)
plot( t_edp, edp_ql.E_dsl(2,:), ...
      t_edi, edi_ql.E_dmpa(2,:), ...
      t_fg,  E_VxB(2,:),         ...
      t_fg,  E_CoRot(2,:) );
ylabel( {'Ey', '(mV/m)'} );

% Ey
subplot(4,1,4)
plot( t_edp, edp_ql.E_dsl(3,:), ...
      t_edi, edi_ql.E_dmpa(3,:), ...
      t_fg,  E_VxB(3,:),         ...
      t_fg,  E_CoRot(3,:) );
xlabel( ['Time (seconds since ' tstart] );
ylabel( {'Ez', '(mV/m)'} );

%------------------------------------%
% Plot B, V_ExB                      %
%------------------------------------%

%
% BCS
%
fig_vdrift = figure();

% Magnetic Field
subplot(4,1,1)
plot( t_fg, fg_ql.b_dmpa );
title([ 'Magnetic Field and Drift Velocity'] );
ylabel( {'B', '(nT)'} );
legend( {'B_{x}', 'B_{y}', 'B_{z}', '|B|'} );

% Vx
subplot(4,1,2)
plot( t_edi, edi_ql.v_ExB_dmpa(1,:), ...
      t_fg,  v_ExB(1,:),             ...
      t_fg,  v_CoRot(1,:) );
ylabel( {'Vx', '(km/s)'} );
legend( {'V_{EDI}', 'V_{ExB}', 'V_{CoRot}'} );

% Vy
subplot(4,1,3)
plot( t_edi, edi_ql.v_ExB_dmpa(2,:), ...
      t_fg,  v_ExB(2,:),             ...
      t_fg,  v_CoRot(2,:) );
ylabel( {'Vy', '(km/s)'} );

% Vy
subplot(4,1,4)
plot( t_edi, edi_ql.v_ExB_dmpa(3,:), ...
      t_fg,  v_ExB(3,:),             ...
      t_fg,  v_CoRot(3,:) );
xlabel( ['Time (seconds since ' tstart] );
ylabel( {'Vz', '(km/s)'} );
