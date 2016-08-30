%
% Name
%   mms_bss_roi_read
%
% Purpose
%   Read MMS region of interest (ROI) start and end times from a file. ROI times can
%   be obtained from IDL and SPEDAS using:
%       IDL> print, time_string( mms_get_roi( ['2015-09-18/00:00:00', '2015-09-18/24:00:00'] ) )
%
%   Or can be obtained and written to a file using
%       IDL> unh_bss_roi_write()
%
%   Times are returned as CDF TT2000 times.
%
% Calling Sequence
%   TT2000 = mms_bss_roi_read( FILENAME );
%     Read ROI start and end times from file FILENAME and return them as TT2000 times.
%
% Parameters
%   FILENAME:           in, required, type=char
%
% Returns
%   TT2000:             in, required, type=int64 2xN
%
% MATLAB release(s) MATLAB 8.2.0.701 (R2013b)
% Required Products None
%
% History:
%   2016-05-20      Written by Matthew Argall
%   2016-05-20      Return an array of size 2xN, not Nx2. - MRA
%
function tt2000 = mms_bss_roi_read( filename )

	% File must exist
	assert(exist(filename, 'file') == 2, ['FILENAME does not exist: "' filename '".']);

	% Open the file
	fID  = fopen(filename);

	% Try to read the file
	try
		data = textscan(fID, '%20s %20s', 'HeaderLines', 1);
	catch ME
		mrfprintf('logerr', ME)
		fclose(fID);
		return
	end
	
	% Close the file
	fclose(fID);
	
	% Create string arrays
	t0 = vertcat( data{1}{:} );
	t1 = vertcat( data{2}{:} );
	
	% Reformat time strings
	%   - MrCDF_Epoch_Parse requires 'yyyy-mm-ddTHH:MM:SS.mmmuuunnn'
	t0 = strcat( t0(:,1:19), '.000000000' );
	t1 = strcat( t1(:,1:19), '.000000000' );

	% Convert ISO-8601 to TT2000
	tt2000      = zeros( [2, size( t0, 1 )], 'int64' );
	tt2000(1,:) = MrCDF_Epoch_Parse( t0, 'CDF_TIME_TT2000' );
	tt2000(2,:) = MrCDF_Epoch_Parse( t1, 'CDF_TIME_TT2000' );
end