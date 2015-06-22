% Help comments
%
%

function [ driftStep driftStepSigma driftVelocity drift_E ] = ...
	edi_drift_step ( ...
		B_tt2000, ...
		B_dmpa, ...
		gd_virtual_dmpa, ...
		gd_fv_dmpa, ...
		obsID, ...
		gd_ID );

% 	disp 'Entering edi_drift_step'
	myLibScienceConstants
	plot_edi_drift_step = false;

	% I3   =  [ 1 0 0 ;  0 1 0 ; 0 0 1 ];
	B2n  = norm (B_dmpa, 2);
	Bu   = B_dmpa / B2n;
	% 	OCSz = [ 0.0; 0.0; 1.0 ];
	% 	r    = cross (OCSz, Bu); % rotate OCSz to B
	% 	rx   = [ ...
	% 		  0.0  r(3) -r(2);
	% 		-r(3)   0.0  r(1);
	% 		 r(2) -r(1)   0.0 ];
	% 	cosTheta = Bu (3); % dot (Bu, OCSz): ||OCSz|| = ||Bu|| = 1; % Always Bu (3)
	% 	DMPA2BPP  = rx + (cosTheta * I3) + ((1-cosTheta) * (r*r') / sum (r.^2));

	DMPA2BPPy = [ 0.0; 1.0; 0.0 ];
	DMPA2BPPx = cross (DMPA2BPPy, Bu);
	DMPA2BPPx = DMPA2BPPx / norm (DMPA2BPPx, 2);
	DMPA2BPPy = cross (Bu, DMPA2BPPx);
	DMPA2BPP  = [ DMPA2BPPx'; DMPA2BPPy'; Bu' ];

	B_bpp = DMPA2BPP * B_dmpa;
	% OR
	B_bpp = [ 0.0; 0.0; B2n ];

	% -~-~-~-~-~-~-~-~-~
	% Rotate GDU positions and firing vectors for each GDU into BPP
	gd_virtual_bpp = DMPA2BPP * gd_virtual_dmpa;
	gd_fv_bpp      = DMPA2BPP * gd_fv_dmpa;
% keyboard
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
	gd_m_bpp = zeros (1, length (gd_virtual_dmpa), 'double');
	gd_b_bpp = zeros (1, length (gd_virtual_dmpa), 'double');

	% -~-~-~-~-~-~-~-~-~
	% A gd12 beam originates in gun 1 and is detected in det 2.
	% Find all the beam slopes and y-intercepts
	gd_m_bpp = gd_fv_bpp (2,:) ./ gd_fv_bpp (1,:);
	gd_b_bpp = gd_virtual_bpp (2,:) - gd_m_bpp.* gd_virtual_bpp (1,:);

	% Find the target in BPP, using BPP FV convergence
	% preAlloc beamIntercepts based on nBeams: (nBeams - 1) * nBeams / 2
	% -~-~-~-~-~-~-~-~-~
	nBeams           = uint32 (length (gd_m_bpp));
	nBeamIntercepts  = (nBeams-1) * nBeams / 2;
	beamIntercepts   = zeros (2, nBeamIntercepts, 'double');
	interceptWeights = zeros (1, nBeamIntercepts, 'double');

	driftStep      = [NaN; NaN; NaN];
	driftStepSigma = [NaN; NaN; NaN];
	driftVelocity  = [NaN; NaN; NaN];
	drift_E        = [NaN; NaN; NaN];

	nBeamIntercepts  = 0;
	if nBeams > 3
		for i = 1: nBeams-1
			for j = i+1: nBeams
				XY = [ ...
					-gd_m_bpp(i) 1 ;
					-gd_m_bpp(j) 1 ];
				b = [ gd_b_bpp(i) ; gd_b_bpp(j) ];
				nBeamIntercepts = nBeamIntercepts + 1;
				beamIntercepts (:, nBeamIntercepts) = XY \ b;
			end
		end

		macroBeamCheckAngle = atan(tand(5)); % Beams closer than 5° get zero intersection weight
		nBeamIntercepts = 0.0;
		for i = 1: nBeams-1
			for j = i+1: nBeams
				nBeamIntercepts = nBeamIntercepts + 1;

				interceptAngle (nBeamIntercepts) = atan ( ...
					(gd_m_bpp(j) - gd_m_bpp(i)) / ...
					(1.0 + gd_m_bpp(i) * gd_m_bpp(j)) );
% [interceptAngle(nBeamIntercepts)  macroBeamCheckAngle]
				if (abs (interceptAngle (nBeamIntercepts)) > macroBeamCheckAngle)
% disp 'abs (interceptAngle (nBeamIntercepts)) > macroBeamCheckAngle'
					interceptWeights (1, nBeamIntercepts) = sin (interceptAngle (nBeamIntercepts))^2;
% 					interceptWeights (1, nBeamIntercepts)
				else
					interceptWeights (1, nBeamIntercepts) = NaN;
				end
			end
		end
% size(interceptWeights)
% nansum (interceptWeights)
		if nansum (interceptWeights) > 0.0
% disp 'nansum (interceptWeights) > 0.0'
			% We will check the upper bound, P0=84%, alpha=P1=0.16 ~> +- 0.08 (lower|upper),
			% so div by 2.0, and pass just the upper bound.
			P0 = 0.84;
			edi_stats_alpha = (1.0 - P0) / 2.0;

			[ GrubbsBeamInterceptMean(1,1), GrubbsBeamInterceptStdDev(1,1), GrubbsBeamIntercepts ] = ...
				edi_drift_step_Grubbs (beamIntercepts (1,:), interceptWeights, edi_stats_alpha);
			ibx = find (isnan (GrubbsBeamIntercepts));

			[ GrubbsBeamInterceptMean(2,1), GrubbsBeamInterceptStdDev(2,1), GrubbsBeamIntercepts ] = ...
				edi_drift_step_Grubbs (beamIntercepts (2,:), interceptWeights, edi_stats_alpha);
			iby = find (isnan (GrubbsBeamIntercepts));

			iIntersectXYOutliers = union (ibx, iby);
			GrubbsBeamIntercepts = beamIntercepts;
			GrubbsBeamIntercepts (:, iIntersectXYOutliers) = [];

			nGrubbsBeamIntercepts          = length (GrubbsBeamIntercepts);
			GrubbsBeamInterceptMean_stdDev = GrubbsBeamInterceptStdDev / sqrt (nGrubbsBeamIntercepts); % x,y mean std dev
%			disp ( sprintf ('nGrubbsBeamIntercepts %4d : GrubbsInterceptMean, StdDev: (%+7.3f, %+7.3f) (%+7.3f, %+7.3f)', ...
%				nGrubbsBeamIntercepts, GrubbsBeamInterceptMean, GrubbsBeamInterceptStdDev) )

			% -~-~-~-~-~-~-~-~-~
			% now we need the drift step...
			driftStep_bpp = [ GrubbsBeamInterceptMean(1); GrubbsBeamInterceptMean(2); 0.0 ];
			gyroFrequency = (q * B2n * nT2T) / e_mass; % (SI) (|q| is positive here.)
			gyroPeriod    = (twoPi / gyroFrequency);    % (SI) The result is usually on the order of a few ms
			% driftStep * B^2 / gyroPeriod = E_bpp x B_bpp =>
			% 	E_bpp = cross ( (driftStep * B2n*B2n / gyroPeriod), B_bpp) / dot (B_bpp, B_bpp);
			%OR
			% 	E_bpp = cross ( (driftStep * B2n*B2n / gyroPeriod), B_bpp) / B2n*B2n;
			%OR
			E_bpp = cross ( (driftStep_bpp          / gyroPeriod), B_bpp*1.0e-9); % B_bpp is in nT, and all these calcs are in SI
			driftStep      = (DMPA2BPP' * driftStep_bpp);
			driftStepSigma = (DMPA2BPP' * [GrubbsBeamInterceptStdDev; 0.0]);
			driftVelocity  = driftStep * gyroPeriod;
			drift_E        = (DMPA2BPP' * E_bpp) * 1.0e3; % convert V/m -> mV/m

			% -~-~-~-~-~-~-~-~-~
			if plot_edi_drift_step
				edi_drift_step_plot ( ...
					obsID, ...
					B_tt2000, ...
					gd_virtual_bpp, ...
					gd_fv_bpp, ...
					DMPA2BPP, ...
					GrubbsBeamIntercepts, GrubbsBeamInterceptMean, GrubbsBeamInterceptMean_stdDev, ...
					P0);
			end
		end
	end
end
