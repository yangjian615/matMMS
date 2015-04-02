%
% Name
%   mms_anc_read_ephem
%
% Purpose
%   Read definitive attitude file.
%
% Calling Sequence
%   [EPHEM] = mms_anc_read_ephem(FILENAME)
%     Read ephemeris data from file FILENAME.
%
%   [DEFATT, META] = mms_anc_read_deffat(FILENAME)
%     Also return the file's metadata. Each line of the file is stored as
%     an element of a cell array.
%
% Parameters
%   FILENAME        in, required, type=char
%
% Returns
%   EPHEM           out, required, type=1x8 cell
%                   Cells are:
%                      1  Time           (UTC) 'yyyy-ddd/HH:MM:SS.sss'
%                      2  Epoch MMS TAI  Days elapsed since 1958-001/00:00:00 UTC
%                      3  X Km           Spacecraft position
%                      3  Y Km           Spacecraft position
%                      3  Z Km           Spacecraft position
%                      3  X Km/s         Spacecraft velocity
%                      3  Y Km/s         Spacecraft velocity
%                      3  Z Km/s         Spacecraft velocity
%                      3  Mass Kg        Spacecraft mass
%   META            out, optional, type=cell
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-03-25      Written by Matthew Argall
%
function [ephem, meta] = mms_anc_read_ephem(filename)

	% File must exist.
	assert( exist(filename, 'file') == 2, ...
		      ['Ephemeris file does not exist: "' filename '".']);

	% Open the file
	fID = fopen(filename);

%------------------------------------%
% Skip Over Header                   %
%------------------------------------%
	
	% Find the first empty line
	line    = 'temp';
	count   = 1;
	meta    = cell(1, 100);
	meta(:) = {''};
	while ~isempty(line) &&  isempty(regexp(line, '^\s*$', 'once'))
		line        = fgetl(fID);
		meta{count} = line;
		count       = count + 1;
	end
	meta = meta(1:count-1);
	
	% Read the column headers
	line = fgetl(fID);
	line = fgetl(fID);

%------------------------------------%
% Read the Data                      %
%------------------------------------%
	format = '%s %d %d %d %d %d %d %d %d';
	ephem = textscan(fID, format, 'MultipleDelimsAsOne', true);
	
	% Close the file
	fclose(fID);
end

% Definitive Orbit Ephemeris
% Spacecraft = MMS2
% StartTime = 2015-073/03:00:25.000 UTC  (MMS TAI: 20891.125694445)
% StopTime = 2015-074/05:28:55.000 UTC  (MMS TAI: 20892.228819445)
% MMS TAI Reference Epoch = 1958-001/00:00:00 UTC
% CentralBody = Earth
% ReferenceFrame = Mean of J2000
% PrincipalPlane = Equatorial
% Project = PPFA_ProductGeneration_DefinitiveEphem.m
% Source = FDGSS OPS
% FileCreationDate = Mar 16 2015 04:15:12.290 UTC
%  
% Epoch (UTC)             Epoch MMS TAI         X                       Y                       Z                       VX                      VY                      VZ                      Mass             
%                                               Km                      Km                      Km                      Km/Sec                  Km/Sec                  Km/Sec                  Kg               
% 2015-073/03:00:25.000   20891.125694445       8947.173812221          -8931.236682260         -6773.106446976         0.7986656112000         6.1909841475000         2.6588927471000         1353.6114694391  


