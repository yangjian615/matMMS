function [ dataMean, dataSigma, dataOut ] = GrubbsTestForOutliers (dataIn, thisAlpha);
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
	assert (nArgsIn == 2, 'GrubbsTestForOutliers: args: data, alpha.');

	dataOut = dataIn;

	while 1
		dataMean  = nanmean (dataOut);
		dataSigma = nanstd (dataOut);
		diffFromMean = abs (dataOut - dataMean);

		[ maxVal iMax ] = max (diffFromMean);
		Gtest = maxVal / dataSigma;
		critical_Z = tDistroCriticalZ (thisAlpha, length (dataOut));

		if Gtest > critical_Z
			dataOut (iMax) = NaN;
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
	tDistCritSq = tinv (thisAlpha / DoF, DoF - 2.0)^2.0;
	critical_Z = (DoF-1.0) * (sqrt (tDistCritSq / (DoF * (DoF - 2.0 + tDistCritSq)) ));
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
