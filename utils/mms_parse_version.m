%
% Name
%   mms_parse_version
%
% Purpose
%   Parse MMS file versions.
%
% Calling Sequence
%   VX = mms_parse_version(VERSION);
%     Parse MMS file versions VERSION and return the x-version number, VX.
%      VERSION is a string or cell array of strings. VX is the same type
%      and length as VERSION.
%
%   [..., VY, VZ] = mms_parse_version(__);
%     Also return the y- and z-version numbers.
%
% Parameters
%   VERSION:        in, required, type=char/cell
%
% Outputs:
%   VX:             out, optional, type=char/cell
%   VY:             out, optional, type=char/cell
%   VZ:             out, optional, type=char/cell
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2016-04-02      Written by Matthew Argall
%
function [vx, vy, vz] = mms_parse_version(version)
	
	% If input is a string, convert it to a cell
	if ~ischar(version) && ~isrow(version)
		error('MMS:badDatatype', 'Version must be a scalar or cell array of char rows.')
	end
	
	% Parse string
	vstr = regexp( version, '^([0-9]+)\.([0-9]+)\.([0-9]+)$', 'tokens' );
	
	% Un-nest cells.
	%   - There is one extra layer if VERSION is a cell array.
	vstr = vertcat( vstr{:} );
	if iscell(version)
		vstr = vertcat( vstr{:} );
		vx   = vstr(:,1);
		vy   = vstr(:,2);
		vz   = vstr(:,3);
	else
		vx   = vstr{:,1};
		vy   = vstr{:,2};
		vz   = vstr{:,3};
	end
end