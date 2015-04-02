%
% Name
%   mms_instr_origins_instr
%
% Purpose
%   Return the origins of various instruments relative to the origin of
%   another instrument.
%
% Calling Sequence
%   ORIGIN = mms_instr_origins_ocs(INSTR1, INSTR2);
%			Return the cartesian coordinates [x, y, z] of the origin of the instrument
%			named INSTR1 relative to the orgigin of the instrument named INSTR2.
%			For valid instrument names, see mms_instr_origins_ocs.m
%
%   ORIGIN = mms_instr_origins_ocs(__, 'Spherical', TF);
%			Specify whether the origin should be returned in spherical or
%			cartesian coordinates. If spherical coordinate are returned, ORIGINS
%			is [azimuth, elevation, radius], where azimuth is the angle from
%			INSTR2 x-axis, elevation is the angle up from the OCS xy-plane, and
%			radius is the radial distance from INSTR2 to INSTR1.
%
% Parameters
%   INSTR1          in, required, type=char
%   INSTR2          in, required, type=char
%   'Spherical'     in, optional, type=boolean, default=false
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-03-22      Written by Matthew Argall
%
function origin = mms_instr_origins_instr(instr1, instr2, varargin)
	
	% Default to cartesian coordinates
	spherical = false;
	
	% Optional inputs
	nOptArg = length(varargin);
	if nOptArg == 2
		if strcmp(varargin{1}, 'Spherical')
			spherical = varargin{2};
		else
			error(['Unknown input "' varargin{1} '".']);
		end
	elseif nOptArg ~= 0
		error('Incorrect number of parameters.');
	end
		
	% Get the instrument origins
	instr1_origin = mms_instr_origins_ocs(instr1);
	instr2_origin = mms_instr_origins_ocs(instr2);

	% Position of Instr1's origin with respect to the origin of Instr2
	origin = instr1_origin - instr2_origin;
	
	% Convert to spherical coordinates
	if spherical
		[az, el, r] = cart2sph(origin(1), origin(2), origin(3));
		origin      = [az, el, r];
	end
end