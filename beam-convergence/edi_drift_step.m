% Help comments
%
%

% function [driftStep, uncertainty] = edi_drift_step ( ...
% 	beams (5s): firingVectors_gse, quality, runnerOrder, position;
% 	mag: Bfield_avg_gse

function [E_dmpa, driftVelocity_dmpa, driftStep_dmpa] = edi_drift_step ( ...
		b_tt2000, ...
		b_avg_dmpa, ...
		edi_gun1_virtual_dmpa, ...
		edi_gd12_fv_dmpa, ...
		edi_gun2_virtual_dmpa, ...
		edi_gd21_fv_dmpa, ...
		edi_gd12_tof, ...
		edi_gd21_tof, ...
		E_dmpa)

	common_science_constants
	plot_edi_drift_step = false;

	I3   =  [ 1 0 0 ;  0 1 0 ; 0 0 1 ];
	B2n  = norm (b_avg_dmpa, 2);
	Bu   = b_avg_dmpa / B2n;
	OCSz = [ 0.0; 0.0; 1.0 ];
	r    = cross (OCSz, Bu); % rotate OCSz to B
	rx   = [ ...
		  0.0  r(3) -r(2);
		-r(3)   0.0  r(1);
		 r(2) -r(1)   0.0 ];
	cosTheta = Bu (3); % dot (Bu, OCSz): ||OCSz|| = ||Bu|| = 1; % Always Bu (3)
	DMPA2BPP  = rx + (cosTheta * I3) + ((1-cosTheta) * (r*r') / sum (r.^2));
	B_bpp = DMPA2BPP * b_avg_dmpa;
	% OR
	B_bpp = [ 0.0; 0.0; B2n ];

	% -~-~-~-~-~-~-~-~-~
	% Rotate GDU positions and firing vectors for each GDU into BPP
	edi_gun1_virtual_bpp = DMPA2BPP * edi_gun1_virtual_dmpa;
	edi_gun2_virtual_bpp = DMPA2BPP * edi_gun2_virtual_dmpa;
	edi_gd12_fv_bpp      = DMPA2BPP * edi_gd12_fv_dmpa;
	edi_gd21_fv_bpp      = DMPA2BPP * edi_gd21_fv_dmpa;

	% -~-~-~-~-~-~-~-~-~
	% Find the most probable beam convergence for the drift step "target", S*
	% Theory:
	% The beams were rotated into BPP. Those beams are parallel to BPP (perpendicular to B),
	% though shifted in BBPz according to the tilt of the spacecraft in BPP.
	% If we assume a value of zero for the BPPz component of the beam, then we have but to solve
	% for the intersections of the beams in 2D. Here we gather the info that we'll need later:
	% slope and y-intercept. After processing all the beams, we'll calculate the intersections.
	% y = mx + b => m = dy/dx =; b = y - m*x, where x,y = GDU pos

	% pre-alloc?
	gd12_m_bpp = zeros (1, length (edi_gun1_virtual_dmpa), 'double');
	gd12_b_bpp = zeros (1, length (edi_gun1_virtual_dmpa), 'double');
	gd21_m_bpp = zeros (1, length (edi_gun2_virtual_dmpa), 'double');
	gd21_b_bpp = zeros (1, length (edi_gun2_virtual_dmpa), 'double');

	% -~-~-~-~-~-~-~-~-~
	% A gd12 beam originates in gun 1 and is detected in det 2.
	% Find all the beam slopes and y-intercepts
	gd12_m_bpp = edi_gd12_fv_bpp (2,:) ./ edi_gd12_fv_bpp (1,:);
	gd12_b_bpp = edi_gun1_virtual_bpp (2,:) - gd12_m_bpp.* edi_gun1_virtual_bpp (1,:);

	gd21_m_bpp = edi_gd21_fv_bpp (2,:) ./ edi_gd21_fv_bpp (1,:);
	gd21_b_bpp = edi_gun2_virtual_bpp (2,:) - gd21_m_bpp.* edi_gun2_virtual_bpp (1,:);

	% All of the beams contribute to the drift step; concatenate the data

	targetBeam2DSlope  = [ gd12_m_bpp, gd21_m_bpp];
	targetBeam2Db      = [ gd12_b_bpp, gd21_b_bpp];

	% Find the target in BPP, using BPP FV convergence
	% preAlloc beamIntercepts based on nBeams: (nBeams - 1) * nBeams / 2

	% -~-~-~-~-~-~-~-~-~
	nBeams          = uint32 (length (targetBeam2DSlope));
	beamIntercepts  = zeros (2, (nBeams-1) * nBeams / 2, 'double');
	nBeamIntercepts = 0;

	for i = 1: nBeams-1
		for j = i+1: nBeams
			XY = [ ...
				-targetBeam2DSlope(i) 1 ;
				-targetBeam2DSlope(j) 1 ];
			b = [ targetBeam2Db(i) ; targetBeam2Db(j) ];
			nBeamIntercepts = nBeamIntercepts + 1;
			beamIntercepts (:, nBeamIntercepts) = XY \ b;
		end
	end

	edi_stats_alpha = (1.0 - 0.84) / 2.0;

	[ ~, ~, GrubbsBeamIntercepts ] = ... 
		edi_drift_step_Grubbs (beamIntercepts (1,:), edi_stats_alpha);
	ibx = find (isnan (GrubbsBeamIntercepts));

	[ ~, ~, GrubbsBeamIntercepts ] = ... 
		edi_drift_step_Grubbs (beamIntercepts (2,:), edi_stats_alpha);
	iby = find (isnan (GrubbsBeamIntercepts));

	iIntersectXYOutliers = union (ibx, iby);
	GrubbsBeamIntercepts = beamIntercepts;
	GrubbsBeamIntercepts (:, iIntersectXYOutliers) = [];

	nGrubbsBeamIntercepts          = length (GrubbsBeamIntercepts);
	GrubbsBeamInterceptMean        = mean (GrubbsBeamIntercepts, 2);
	GrubbsBeamInterceptStdDev      = std (GrubbsBeamIntercepts, 1, 2);
	GrubbsBeamInterceptMean_stdDev = GrubbsBeamInterceptStdDev / sqrt (nGrubbsBeamIntercepts); % x,y mean std dev
%	disp ( sprintf ('nGrubbsBeamIntercepts %4d : GrubbsInterceptMean, StdDev: (%+7.3f, %+7.3f) (%+7.3f, %+7.3f)', ...
%		nGrubbsBeamIntercepts, GrubbsBeamInterceptMean, GrubbsBeamInterceptStdDev) )

	% -~-~-~-~-~-~-~-~-~
	% -~-~-~-~-~-~-~-~-~
	% -~-~-~-~-~-~-~-~-~
	% now we need the drift step...
	driftStep_bpp = [ GrubbsBeamInterceptMean(1); GrubbsBeamInterceptMean(2); 0.0 ];
	gyroFrequency = (q * B2n * nT2T) / e_mass; % (SI) (|q| is positive here.)
	gyroPeriod    = (twoPi / gyroFrequency);    % (SI) The result is usually on the order of a few ms
	
	driftStep_dmpa = (DMPA2BPP' * driftStep_bpp);
	
	% -~-~-~-~-~-~-~-~-~
	% -~-~-~-~-~-~-~-~-~
	% -~-~-~-~-~-~-~-~-~
	% Drift Velocity
	driftVelocity_bpp  = driftStep_bpp * gyroPeriod;
	driftVelocity_dmpa = (DMPA2BPP' * driftVelocity_bpp);
	
	
	% -~-~-~-~-~-~-~-~-~
	% -~-~-~-~-~-~-~-~-~
	% -~-~-~-~-~-~-~-~-~
	% Electric Field
	
	% driftStep * B^2 / gyroPeriod = E_bpp x B_bpp =>
	% 	E_bpp = cross ( (driftStep * B2n*B2n / gyroPeriod), B_bpp) / dot (B_bpp, B_bpp);
	%OR
	% 	E_bpp = cross ( (driftStep * B2n*B2n / gyroPeriod), B_bpp) / B2n*B2n;
	%OR
	E_bpp = cross ( (driftStep_bpp          / gyroPeriod), B_bpp*1.0e-9); % B_bpp is in nT, and all these calcs are in SI
	E_dmpa = (DMPA2BPP' * E_bpp);
	E_dmpa = E_dmpa * 1.0e3; % convert V/m -> mV/m... later, try to combine scaling into constant
	% -~-~-~-~-~-~-~-~-~
	% -~-~-~-~-~-~-~-~-~
	% -~-~-~-~-~-~-~-~-~

	if plot_edi_drift_step
		edi_drift_step_plot ( ...
			b_tt2000, ...
			edi_gun1_virtual_bpp, edi_gun2_virtual_bpp, ...
			edi_gd12_fv_bpp, edi_gd21_fv_bpp, ...
			DMPA2BPP, ...
			GrubbsBeamInterceptMean, GrubbsBeamInterceptMean_stdDev);
	end

end