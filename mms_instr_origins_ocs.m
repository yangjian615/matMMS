%
% Name
%   mms_instr_origins_ocs
%
% Purpose
%   Return the origins of various instruments relative to the Observatory
%   Coordinate System (OCS) origin.
%
% Calling Sequence
%   NAMES = mms_instr_origins_ocs();
%			Return a cell array of strings NAME containing the names of the
%			instruments for which the origins are known.
%
%   ORIGIN = mms_instr_origins_ocs(INSTRUMENT);
%			Return the cartesian coordinates of an instrument's origin ORIGIN in
%			OCS.
%
%   ORIGIN = mms_instr_origins_ocs(__, 'Spherical', TF);
%			Specify whether the origin should be returned in spherical or
%			cartesian coordinates. If spherical coordinate are returned, ORIGIN
%			is [azimuth, elevation, radius], where azimuth is the angle from
%			the OCS x-axis, elevation is the angle up from the OCS xy-plane, and
%			radius is the radial distance from the OCS origin.
%
% Parameters
%   INSTRUMENT      in, optional, type=char
%   'Spherical'     in, optional, type=boolean, default=false
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-03-22      Written by Matthew Argall
%
function origin = mms_instr_origins_ocs(instrument, varargin)

	
	% If an instrument was given
	if nargin > 0
		% Ensure only one instrument name was given
		assert(ischar(instrument) && isrow(instrument), ...
					 'Input instrument must be a scalar string.');
		
		% Case insensitive -- all uppercase
		instr = upper(instrument);
	end
	
	% Defaults
	spherical = false;

	% Which optional inputs were given
	nOptArgs = length(varargin);
	if nOptArgs > 0
		% Name-value pairs -- every other is a name
		for ii = 1 : 2 : nOptArgs
			switch varargin{ii}
				case 'Spherical'
					spherical = varargin{ii+1};
				otherwise
					error(['Parameter name not recognized: "' varargin{ii} '".']);
			end
		end
	end

  % Instrument list
	instr_map                  = containers.Map('KeyType', 'char', 'ValueType', 'any');
	%                                    X           Y            Z
	instr_map('ADP1')          = [  0.0,       -0.161925,    15.745      ];
	instr_map('ADP2')          = [  0.0,       -0.161925,   -15.185      ];
	instr_map('AFG_BOOM')      = [ -0.99147,   -0.99147,     -0.0771     ];
	instr_map('AFG_MECH')      = [  5.18785,   -0.0080,       0.0021262  ];    % From AFG_BOOM origin
	instr_map('DFG_BOOM')      = [ -0.99147,    0.99147,     -0.0771     ];
	instr_map('DFG_MECH')      = [  5.18785,   -0.0080,       0.0021262  ];    % From DFG_BOOM origin
	instr_map('EDI1')          = [ -1.332748,   0.889069,     1.051      ];
	instr_map('EDI1_GUN')      = [ -1.45598,    1.11837,      0.0        ];
	instr_map('EDI1_DETECTOR') = [ -1.35885,    1.03395,      0.0        ];
	instr_map('EDI2')          = [  1.332748,  -0.889069,     1.051      ];
	instr_map('EDI2_GUN')      = [  1.45598,   -1.11837,      0.0        ];
	instr_map('EDI2_DETECTOR') = [  1.35885,   -1.03395,      0.0        ];
	instr_map('IDCS')          = [  0.0,        0.0,          1.0510     ];
	instr_map('DSS')           = [  1.01764,    1.25506,      0.127220   ];
	instr_map('OCS')           = [  0.0,        0.0,          0.0        ];
	instr_map('SC')            = [  0.0,        0.0,          0.1670     ];
	instr_map('SCM_BOOM')      = [ -0.99147,   -0.99147,     -0.0771     ];
	instr_map('SCM_MECH')      = [  4.14785,   -0.00479899,  -0.0332010  ];    % From SCM_BOOM origin
	instr_map('SDP1')          = [  1.342598,   0.865542,     1.050      ];
	instr_map('SDP2')          = [ -1.342598,  -0.865542,     1.050      ];
	instr_map('SDP3')          = [ -0.865542,   1.342598,     1.050      ];
	instr_map('SDP4')          = [  0.865542,  -1.342598,     1.050      ];
	
	% Convert to OCS origin
	instr_map('AFG_123') = instr_map('AFG_MECH') + instr_map('AFG_BOOM');
	instr_map('AFG_XYZ') = instr_map('AFG_MECH') + instr_map('AFG_BOOM');
	instr_map('DFG_123') = instr_map('DFG_MECH') + instr_map('DFG_BOOM');
	instr_map('DFG_XYZ') = instr_map('DFG_MECH') + instr_map('DFG_BOOM');
	instr_map('SCM_XYZ') = instr_map('SCM_MECH') + instr_map('SCM_BOOM');
	instr_map('SCM_123') = instr_map('SCM_MECH') + instr_map('SCM_BOOM');
	
	% Remove non-OCS origins
	instr_map.remove({'AFG_MECH' 'DFG_MECH' 'SCM_MECH'});

	% Return only the instrument names?
	if nargin == 0
		origin = instr_map.keys();
		return
	end

	% Case insensitive version
	tf_has = instr_map.isKey(instr);
	if ~tf_has
		error( ['Invalid instruments given: "' instrument(tf_has == 0)] );
	end
	
	% Retrieve the desired origin(s).
	origin = instr_map(instr);

	% Return sperical coordinates?
	if spherical
		[az, el, r] = cart2sph(origin(1), origin(2), origin(3));
		origin      = [az el r];
	end
end