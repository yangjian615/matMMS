%
% Name
%   mms_fsm_nf_add
%
% Purpose
%   Add an entry to the FSM noise floor calibration data.
%
% Calling Sequence
%   DATA = mms_fsm_nf_add(OLD, NEW)
%     Add an entry to the noise floor calibration. OLD contains all data from
%     all previous noise floor calibration entries. NEW is a new entry to be
%     spliced into or added to OLD. The result is returned in the structure DATA.
%
% Parameters
%   OLD:            in, required, type=struct
%                   Structure with the following fields:
%                     't'     -  Time (cdf_time_tt2000)
%                     'dt'    -  Time interval over which the noise floor is valid.
%                     'nf'    -  Noise floor
%                     'std'   -  Standard deviation of the noise floor
%                     'flag'  -  Operational flag
%   NEW:            in, required, type=struct
%                   Structure with the following fields:
%                     't'     -  Time (cdf_time_tt2000)
%                     'dt'    -  Time interval over which the noise floor is valid.
%                     'nf'    -  Noise floor
%                     'std'   -  Standard deviation of the noise floor
%                     'flag'  -  Operational flag
%
% Returns
%   DATA            out, required, type=struct
%                   Structure with the following fields:
%                     't'     -  Time (cdf_time_tt2000)
%                     'dt'    -  Time interval over which the noise floor is valid.
%                     'nf'    -  Noise floor
%                     'std'   -  Standard deviation of the noise floor
%                     'flag'  -  Operational flag
%                     'comp'  -  Component index (0=x, 1=y, 2=z)
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2016-10-23      Written by Matthew Argall
%
function data = mms_fsm_nf_add(old, new)

	tf_replace = false;
	tf_overlap = false;
	assert( length( new.t ) == 1, 'NEW must contain a single entry.' )
	
%------------------------------------------------------
% Determine where the entry fits                      |
%------------------------------------------------------
	% Number of points
	N_old = length( old.t );
	
	% Convert TT2000 to seconds
	t_ref     = min( old.t(1), new.t(1) );
	t_old_ssm = MrCDF_epoch2ssm( old.t, t_ref );
	t_new_ssm = MrCDF_epoch2ssm( new.t, t_ref );
	
	% Locate T_NEW within T_OLD
	it = MrValue_Locate( t_old_ssm, t_new_ssm );

%------------------------------------------------------
% Look for Overlap                                    |
%------------------------------------------------------
	if it == 0
		it0 = 1;
	elseif it > N_old;
		it0 = N_old;
	else
		it0 = it;
	end
	if N_old > it0
		it1 = it0+1;
	else
		it1 = it0;
	end

	dt0 = old.t(it0) + old.dt(it0) - new.t;
	dt1 = new.t + new.dt - old.t(it1);

	% Equal to an existing point
	if old.t(it) == new.t
		tf_replace = true;
		mrfprintf( 'logwarn', 'Overwriting old noise floor value.' );
	
	% Overlaps previous point
	elseif dt0 > 0
		% How much overlap?
		dt0   = dt0 * 1e-9;
		nDays = floor( dt1 / 86400.0 );
		nHrs  = floor( mod(dt1, 86400.0) / 3600.0 )
		nMin  = floor( mod(dt1, 3600.0) / 60.0 )
		nSec  = mod(dt1, 60.0)
		
		% Print warning
		tf_overlap = true;
		mrfprintf( 'logwarn', 'Overlaps previous point by %02i days %02i hrs %02i min %07.4f sec.', ...
		           nDays, nHrs, nMin, nSec );
	
	% Overlaps next point
	elseif dt1 > 0
		% How much overlap?
		dt1   = dt1 * 1e-9;
		nDays = floor( dt1 / 86400.0 );
		nHrs  = floor( mod(dt1, 86400.0) / 3600.0 )
		nMin  = floor( mod(dt1, 3600.0) / 60.0 )
		nSec  = mod(dt1, 60.0)
		
		% Print warning
		tf_overlap = true;
		mrfprintf( 'logwarn', 'Overlaps future point by %02i days %02i hrs %02i min %07.4f sec.', ...
		           nDays, nHrs, nMin, nSec );
	end

%------------------------------------------------------
% Splice                                              |
%------------------------------------------------------

keyboard
	if tf_replace
		% Overwrite value
		old.t(it)         = new.t;
		old.dt(it)        = new.dt;
		old.f(:,it)       = new.f;
		old.nf(:,:,:,it)  = new.nf;
		old.std(:,:,:,it) = new.std;
		old.flag(:,it)    = new.flag;
		
		% Create output
		data = old;
	
	elseif it == 0
		data = struct( 't',     cat( 2, new.t,    old.t    ), ...
		               'f',     cat( 2, new.f,    old.f    ), ...
		               'dt',    cat( 2, new.dt,   old.dt   ), ...
		               'nf',    cat( 4, new.nf,   old.nf   ), ...
		               'std',   cat( 4, new.std,  old.std  ), ...
		               'flag',  cat( 2, new.flag, old.flag ) );
	
	elseif it == N_old
		data = struct( 't',     cat( 2, old.t,    new.t     ), ...
		               'dt',    cat( 2, old.dt,   new.dt    ), ...
		               'f',     cat( 2, old.f,    new.f     ), ...
		               'nf',    cat( 2, old.nf,   new.nf    ), ...
		               'std',   cat( 2, old.std,  new.std   ), ...
		               'flag',  cat( 2, old.flag, new.flag  ) );
	
	else
		data = struct( 't',     cat( 2, cat( 2, old.t(1:it),         new.t    ), old.t(it+1:end)         ), ...
		               'dt',    cat( 2, cat( 2, old.dt(1:it),        new.dt   ), old.dt(it+1:end)        ), ...
		               'f',     cat( 2, cat( 2, old.f(:,1:it),       new.f    ), old.f(:,it+1:end)       ), ...
		               'nf',    cat( 4, cat( 4, old.nf(:,:,:,1:it),  new.nf   ), old.nf(:,:,:,it+1:end)  ), ...
		               'std',   cat( 4, cat( 4, old.std(:,:,:,1:it), new.std  ), old.std(:,:,:,it+1:end) ),  ...
		               'flag',  cat( 2, cat( 2, old.flag(:,1:it),    new.flag ), old.flag(:,it+1:end)    ) );
	end
end