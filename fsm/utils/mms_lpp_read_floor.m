%
% Name
%   mms_lpp_read_floor
%
% Purpose
%   Read SCM noise floor measurements.
%
% Calling Sequence
%   DATA = mms_lpp_read_floor()
%     Read in-situ and ground-based measurments of the SCM noise floor.
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
%   2016-10-04      Written by Matthew Argall
%
%***************************************************************************
function data = mms_lpp_read_floor()

	% Directory with noise-floor data
	dir = '/home/argall/data/lpp_scm/';

%------------------------------------%
% In-Situ Noise Floor                %
%------------------------------------%
	
	mms4_bz_floor_srvy = fullfile(dir, 'insitu_floor', 'noise_mms4_z_201606_srvy.txt');
	mms4_bz_floor_brst = fullfile(dir, 'insitu_floor', 'noise_mms4_z_201606_brst.txt');
	
	%
	% MMS4
	%
	
	% Srvy
	fID     = fopen( mms4_bz_floor_srvy );
	bz_srvy = textscan( fID, '%9f %13f', ...
	                    'HeaderLines', 1 );
	fclose(fID);
	
	% Brst
	fID     = fopen( mms4_bz_floor_brst );
	bz_brst = textscan( fID, '%9f %13f', ...
	                    'HeaderLines', 1 );
	fclose(fID);
	
	%
	% Collect data
	%
	
	% Create a structure
	data = struct( 'f_srvy_insitu',  bz_srvy{1}', ...
	               'bz_srvy_insitu', bz_srvy{2}', ...
	               'f_brst_insitu',  bz_brst{1}', ...
	               'bz_brst_insitu', bz_brst{2}'  ...
	             );
	
	clear bz_srvy bz_brst

%------------------------------------%
% EMI Noise Floor                    %
%------------------------------------%
	components = {'x', 'y', 'z'};
	
	% Loop over spacecraft
	for ii = 1 : 4
		b_temp = cell(1,3);
		
		% Flight model
		switch ii
			case 1
				fm = 1;
			case 2
				fm = 2;
			case 3
				fm = 4;
			case 4
				fm = 3;
			otherwise
				error( 'Loop error' );
		end
		
		% Loop over components
		for jj = 1 : 3
			% Create the file name
			fname = sprintf('mms%i_scm%i_nemi_%s.txt', ii, fm, components{jj});
			file  = fullfile(dir, 'nemi_floor', fname);
			
			% Open the file
			fID = fopen(file);

			% Read the data
			temp = textscan( fID, '%11f %11f %11f %11f', ...
			                 'HeaderLines',         2,   ...
			                 'MultipleDelimsAsOne', true );
			
			% Close the file
			fclose(fID);

			% Temporarily store data
			%   - Convert from pT/sqrt(Hz) to nT/sqrt(Hz)
			b_temp{jj} = temp{4} * 1e-3;
			temp{3}    = [];
		end

		% Add to data structure
		f_tag        = sprintf('mms%i_f_nemi', ii);
		b_tag        = sprintf('mms%i_b_nemi', ii);
		data.(f_tag) = temp{1}';
		data.(b_tag) = [ b_temp{1} b_temp{2} b_temp{3} ]';
	end
end
