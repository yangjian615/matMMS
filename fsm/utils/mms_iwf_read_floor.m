%
% Name
%   mms_iwf_read_floor
%
% Purpose
%   Read SCM noise floor measurements.
%
% Calling Sequence
%   DATA = mms_iwf_read_floor()
%     Read ground-based measurments of the DFG noise floor.
%
% Inputs
%   None.
%
% Outputs
%   DATA        out, required, type=struct
%               Structure with the following fields:
%                 'f_srvy_insitu  - Frequency at which noise floor is measured
%                 'bz_srvy_insitu - In-situ mesurments of Bz noise floor from MMS1 in srvy mode
%                 'f_brst_insitu  - Frequency at which noise floor is measured
%                 'bz_brst_insitu - In-situ mesurments of Bz noise floor from MMS1 in brst mode
%                 'mms1_f_nemi'   - Fruqencies at which noise floor is measured during EMI testing
%                 'mms1_b_nemi'   - Noise floor results from EMI testing for MMS1
%                 'mms2_f_nemi'   - Fruqencies at which noise floor is measured during EMI testing
%                 'mms2_b_nemi'   - Noise floor results from EMI testing for MMS2
%                 'mms3_f_nemi'   - Fruqencies at which noise floor is measured during EMI testing
%                 'mms3_b_nemi'   - Noise floor results from EMI testing for MMS3
%                 'mms4_f_nemi'   - Fruqencies at which noise floor is measured during EMI testing
%                 'mms4_b_nemi'   - Noise floor results from EMI testing for MMS4
%
% MATLAB release(s) 9.0.0.341360 (R2016a)
% Required Products None
%
% History:
%   2017G-10-04      Written by Matthew Argall
%
%***************************************************************************
function data = mms_lpp_read_floor()

	% Directory with noise-floor data
	path = fileparts( mfilename('fullpath') );
	file = fullfile( path, '..', '..', '..', 'fischer', 'Crossover ML estimator', 'XYZ_PSD_MMS_DFG_DEC64.txt' );

%------------------------------------%
% Read the Data                      %
%------------------------------------%
	
	% Read
	fID  = fopen( file );
	data = textscan( fID, '%f %f %f %f', ...
	                 'HeaderLines', 1 );
	fclose(fID);
	
	%
	% Collect data
	%
	
	% Create a structure
	data = struct( 'f',  data{1}', ...
	               'Bx', data{2}', ...
	               'By', data{3}', ...
	               'Bz', data{4}'  ...
	             );
end
