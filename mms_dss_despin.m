%
% Name
%   mms_dss_despin
%
% Purpose
%   Transform a vector field from a spinning frame to a despun frame.
%
% Calling Sequence
%   [DESPUN] = mms_dss_despin(T_DSS, TIME, DATA)
%     Take a vector field DATA as a function of time TIME in a spinning
%     reference frame and transform it into DESPUN, the vector field in a
%     despun coordinate system. Spin phase is determined by using the
%     digital sun sensor's time tags T_DSS as bin edges of a histogram.
%     TIME falls into bins marked by the closest T_DSS rounded down from
%     TIME. The spin phase, then, is 2 * pi / ( TIME - T_DSS(iBin) ).
%
%   [DESPUN] = mms_dss_despin(__, 'ParamName', ParamValue)
%
% Parameters
%   T_DSS:          in, required, type=1xN int64/tt2000
%   TIME            in, required, type=1xN int64/tt2000
%   DATA            in, required, type=3xN double
%   'Offset'        in, optional, type=double, default=0.0
%                   A constant angular offset in radians to be added to the
%                     despinning rotation.
%   'Omega'         in, optional, type=double, default=mean(diff(T_DSS))
%                   Spin frequency of the data.
%   'Smooth'        in, optional, type=double, default=false
%                   Force T_DSS to step evenly from T_DSS(1) to T_DSS(end)
%                     in intervals of exactly "Omega".
%   'Unique'        in, optional, type=double, default=false
%                   Take only the unique values of T_DSS. Often, the sun
%                     pulse is reported multiple times per spin. This will
%                     remove duplicate entries.
%
% Returns
%   DESPUN          out, required, type=3xN double
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-03-31      Written by Matthew Argall
%   2015-04-06      Added the 'Unique', and 'Smooth' options. Use the mean, 
%                     not the median, spin frequency. - MRA
%
function despun = mms_dss_despin(t_dss, time, data, varargin)

	% Defaults
	omega  = [];
	offset = 0.0;
	smooth = false;
	uniq   = false;

	% Check for optional arguments
	nOptArgs = length(varargin);
	for ii = 1 : 2: nOptArgs
		switch varargin{ii}
			case 'Omega'
				omega = varargin{ii+1};
			case 'Offset'
				offset = varargin{ii+1};
			case 'Unique'
				uniq = varargin{ii+1};
			case 'Smooth'
				smooth = varargin{ii+1};
			otherwise
				error( ['Unknown parameter "' varargin{ii} '".'] );
		end
	end
	
	% Number of points to despin
	npts = length(time);
	
%------------------------------------%
% Massage Sun Pulse Time             %
%------------------------------------%
	% Take only the unique elements?
	if uniq
		t_dss = unique(t_dss);
	end
	
	% Convert time to seconds
	t_min     = min([time(1) t_dss(1)]);
	t_sec     = MrCDF_epoch2sse(time,  t_min);
	t_sec_dss = MrCDF_epoch2sse(t_dss, t_min);

	% Spin frequency
	if isempty(omega)
		omega = 2 * pi / mean(diff(t_sec_dss));
	end
	
	% Smooth the results
	%   - Any elements of T_SEC that are before T_SEC_DSS(1) or after
	%     T_SEC_DSS(end) will return an index of 0 and cause an error.
	%   - To ensure T_SEC_DSS encompasses all of T_SEC, smooth out to
	%     T_SEC(end) + OMEGA.
	if smooth
		t_sec_dss = t_sec_dss(1) : omega : ( t_sec(end) + (2 * pi / omega) );
		
	% Add extra points even if we are not smoothing.
	elseif t_sec_dss(end) < t_sec(end)
		t_sec_dss = [ t_sec_dss (t_sec_dss(end)+omega):omega:(t_sec(end)+omega) ];
	end
	
%------------------------------------%
% Determine Spin Phase               %
%------------------------------------%
	
	% Histogram data times using sunpulse times as bin edges.
	[~, inds] = histc(t_sec, t_sec_dss);
	
	% Radians into the spin
	phase = omega .* ( t_sec - t_sec_dss(inds) );
	
	% Add the offset
	if offset ~= 0
		phase = phase + offset;
	end
	
%------------------------------------%
% Despin                             %
%------------------------------------%
	% Create the rotation matrix to despin the data.
	%    |  cos(omega)  sin(omega)  0  |
	%    | -sin(omega)  cos(omega)  0  |
	%    |      0           0       1  |
	spin2despun        = zeros(3, 3, npts);
	spin2despun(1,1,:) = cos(phase);
	spin2despun(2,1,:) = sin(phase);
	spin2despun(1,2,:) = -sin(phase);
	spin2despun(2,2,:) = cos(phase);
	spin2despun(3,3,:) = 1;
	
	% Despin each vector
	despun = mrvector_rotate(spin2despun, data);
end