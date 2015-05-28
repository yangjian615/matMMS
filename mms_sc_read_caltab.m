%
% Name
%   mms_scm_read_caltab
%
% Purpose
%   Read calibration tables for the MMS SCM instrument.
%
% Calling Sequence
%   [RESPONSE, FREQS] = mms_scm_read_caltab(FILENAME)
%     Read MMS SCM calibration tables from the file FILENAME and return the
%     x-, y-, and z-sensor response functions RESPONSE as a function of the
%     frequencie FREQS.
%
% Parameters
%   FILENAME        in, required, type=char
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-03-21      Written by Matthew Argall
%
function [response, freqs] = mms_sc_read_caltab(filename)

	% Ensure the file exists
	assert( exist(filename, 'file') == 2, ['File does not exist: "' filename '".']);
	
	% Open the file
	fileID = fopen(filename);

%------------------------------------%
% Read Header                        %
%------------------------------------%
	% Skip over the header for now
	line = '';
	while isempty(regexp(line, '# Frequency', 'ONCE'))
		line = fgetl(fileID);
	end

%------------------------------------%
% X-Sensor                           %
%------------------------------------%
	% Read the calibration table for the X-sensor
	count   = 1;
	n_alloc = 1000;
	x_freq  = zeros(1, n_alloc);
	x_value = zeros(1, n_alloc);
	while true
		% Read a line of data
		line = fgetl(fileID);
		if isempty(line)
			continue
		end
		if ~isempty(regexp(line, '#', 'ONCE'))
			break
		end
		
		% Parse the line of data
		c_line         = textscan(line, '', 'MultipleDelimsAsOne', true);
		x_freq(count)  = c_line{1};
		x_value(count) = c_line{2} + 1i * c_line{3};

		% Increase count
		count = count + 1;
		
		% Allocate more memory?
		if count > n_alloc
			x_freq  = [x_freq  zeros(1, n_alloc)];
			x_value = [x_value complex(zeros(1, n_alloc))];
			n_alloc = 2 * n_alloc;
		end
	end
	
	% Trim the data
	x_freq  = x_freq(1:count-1);
	x_value = x_value(1:count-1);
%------------------------------------%
% Y-Sensor                           %
%------------------------------------%
	% Skip over lines
	line = '';
	while isempty(regexp(line,  '# Frequency', 'ONCE'))
		line = fgetl(fileID);
	end
	
	
	% Read the calibration table for the X-sensor
	count   = 1;
	n_alloc = 1000;
	y_freq  = zeros(1, n_alloc);
	y_value = zeros(1, n_alloc);
	while true
		% Read a line of data
		line = fgetl(fileID);
		if isempty(line)
			continue
		end
		if ~isempty(regexp(line, '#', 'ONCE'))
			break
		end
		
		% Parse the line of data
		c_line         = textscan(line, '', 'MultipleDelimsAsOne', true);
		y_freq(count)  = c_line{1};
		y_value(count) = c_line{2} + 1i * c_line{3};
		
		%Increase count
		count = count + 1;
		
		% Allocate more memory?
		if count > n_alloc
			y_freq  = [y_freq  zeros(1, n_alloc)];
			y_value = [y_value complex(zeros(1, n_alloc))];
		end
	end
	
	% Trim the data
	y_freq  = y_freq(1:count-1);
	y_value = y_value(1:count-1);

%------------------------------------%
% Z-Sensor                           %
%------------------------------------%
	% Skip over lines
	line = '';
	while isempty(regexp(line,  '# Frequency', 'ONCE'))
		line = fgetl(fileID);
	end
	
	
	% Read the calibration table for the X-sensor
	count   = 1;
	n_alloc = 1000;
	z_freq  = zeros(1, n_alloc);
	z_value = zeros(1, n_alloc);
	while true
		% Read a line of data
		line = fgetl(fileID);
		if isempty(line)
			continue
		end
		if ~isempty(regexp(line, 'END', 'ONCE'))
			break
		end
		
		% Parse the line of data
		c_line         = textscan(line, '', 'MultipleDelimsAsOne', true);
		z_freq(count)  = c_line{1};
		z_value(count) = c_line{2} + 1i * c_line{3};
		
		%Increase count
		count = count + 1;
		
		% Allocate more memory?
		if count > n_alloc
			z_freq  = [z_freq  zeros(1, n_alloc)];
			z_value = [z_value complex(zeros(1, n_alloc))];
			n_all
		end
	end
	
	% Trim the data
	z_freq  = z_freq(1:count-1);
	z_value = z_value(1:count-1);

%------------------------------------%
% Combine Output                     %
%------------------------------------%
	freqs    = [x_freq;  y_freq;  z_freq];
	response = [x_value; y_value; z_value];
	
	% Read the
	fclose(fileID);	
end
