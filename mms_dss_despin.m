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
%   'Omega'         in, optional, type=double, default=median(diff(T_DSS))
%                   Spin frequency of the data.
%
% Returns
%   DESPUN          out, required, type=3xN double
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-03-31      Written by Matthew Argall
%
function despun = mms_dss_despin(t_dss, time, data, varargin)

	omega  = [];
	offset = 0.0;

	% Check for optional arguments
	nOptArgs = length(varargin);
	for ii = 1 : 2: nOptArgs
		switch varargin{ii}
			case 'Omega'
				omega = varargin{ii+1};
			case 'Offset'
				offset = varargin{ii+1};
			otherwise
				error( ['Unknown parameter "' varargin{ii} '".'] );
		end
	end
	
%------------------------------------%
% Find Closest Sun Pulse             %
%------------------------------------%
	% Convert time to seconds
	t_min     = min([time(1) t_dss(1)]);
	t_sec     = MrCDF_epoch2sse(time,  t_min);
	t_sec_dss = MrCDF_epoch2sse(t_dss, t_min);
	
	% Histogram merged times using sunpulse times as bin edges.
	[~, inds] = histc(t_sec, t_sec_dss);
	
%------------------------------------%
% Determine Spin Phase               %
%------------------------------------%
	% Spin frequency
	if isempty(omega)
		omega = 2 * pi / median(diff(t_sec_dss));
	end
	
	% Radians into the spin
	phase = omega .* ( t_sec - t_sec_dss(inds)' );
	
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