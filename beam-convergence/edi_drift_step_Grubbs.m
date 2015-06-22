function [ dataMean, dataSigma, dataOut ] = GrubbsTestForOutliers (dataIn, dataWeight, thisAlpha);
% dataIn must be reasonably approximated by a normal distribution.
% References:
% "Procedures for Detecting Outlying Observations in Samples,"
%   F.E. Grubbs; Technometrics, 11-1:1--21; Feb., 1969, and
% http://www.itl.nist.gov/div898/handbook/eda/section3/eda35h1.htm
% Note here... http://www.itl.nist.gov/div898/handbook/eda/section3/eda35h.htm
% "It is not appropriate to apply a test for a single outlier sequentially in order
% to detect multiple outliers."

	nArgsIn    = nargin;
	nVarArgout = nargout;
	assert (nArgsIn == 3, 'GrubbsTestForOutliers: args: data, weights, alpha.');

% 	disp 'Entering GrubbsTestForOutliers'
	dataOut = dataIn;

	while 1
		sumWeights   = nansum (dataWeight);
		dataMean     = nansum ((dataOut .* dataWeight)) / sumWeights;
		diffFromMean = dataOut - dataMean;
		dataSigma    = sqrt (nansum (diffFromMean.^2 .* dataWeight) / sumWeights);
% keyboard
		[ maxVal iMax ] = max (abs (diffFromMean));
		Gtest = maxVal / dataSigma;
		critical_Z = tDistroCriticalZ (thisAlpha, length (dataOut));
% [ Gtest  critical_Z dataSigma]
		if Gtest > critical_Z
			dataOut    (iMax) = NaN;
			dataWeight (iMax) = 0.0; % Could be zero or NaN ~> affect other calcs?
		else
			break;
		end
	end
end

function critical_Z = tDistroCriticalZ (thisAlpha, DoF)
	% Computes the critical z value for rejecting outliers (Grubbs Test)
	% http://www.itl.nist.gov/div898/handbook/eda/section3/eda35h1.htm, 2015-05-28
	% DoF = degrees of freedom
	% tcrit = critical value of the t distribution with (N-2) DoF and a significance level of alpha/(2N).
	% For two-sided tests (see above), the hypothesis of no outliers is rejected if
	% G = max |Y[i - Ymean| / dataSigma   >   zcrit
	% For one-sided tests, use a significance level of level of alpha/N.
	% Critical value for an upper one-tailed test: 2.032: alpha = 0.05, DoF=8 (8 samples) (ref above)
	% clc
	% thisAlpha=0.05
	% DoF=8
	% tDistCritSq = tinv (thisAlpha / DoF, DoF - 2.0)^2.0
	% critical_Z = (DoF-1.0) * sqrt ( tDistCritSq / (DoF * (DoF - 2.0 + tDistCritSq)) )

	tDistCritSq = tinv (thisAlpha / DoF, DoF - 2.0)^2.0;
	critical_Z = (DoF-1.0) * sqrt ( tDistCritSq / (DoF * (DoF - 2.0 + tDistCritSq)) );
end

% Test data
% UMassSpec = [ 140.0 199.31 199.53 200.19 200.82 201.92 201.95 202.18 245.57 inf NaN ];
% UMassSpec = [ 140.0 199.31 199.53 200.19 200.82 201.92 201.95 202.18 245.57 ];
% H0: there are no outliers in the data
% Ha: the maximum value is an outlier
% Test statistic:      G = 2.4687
% Significance level:  alpha = 0.05
% Critical value for an upper one-tailed test:  2.032
% Reject H0 if G > 2.032
% GoodUMassSpec = GrubbsTestForOutliers (UMassSpec, 0.05)

% aa=0:0.1:15;
% weight=1-cosd(aa);
% plot (aa, weight, 'r');
% expWeight10 = 10.^(weight)-1.0;
% expWeight = exp(weight)-1.0;
% hold on
% plot (aa, expWeight10, 'g');
% plot (aa, expWeight, 'b');
