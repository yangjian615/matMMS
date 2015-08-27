%
% Name
%   mms_edi_srvy_combine
%
% Purpose
%   Combine slow and fast survey data.
%
% Calling Sequence
%   SRVY = mms_edi_srvy_combine(SLOW_SRVY, FAST_SRVY)
%     Combine slow and fast survey EDI data structures, SLOW_SRVY
%     and FAST_SRVY, into a single structure, SRVY. SLOW_SRVY and
%     FAST_SRVY must have the same fields, two of which must be
%     'tt2000_gd21' and 'tt2000_gd21', which represent the CDF
%     TT2000 times for beam data from gun 1 and gun 2, respectively.
%
% Parameters
%   SLOW_SRVY:      in, required, type=struct
%   FAST_SRVY:      in, required, type=struct
%
% Returns
%   SRVY            out, required, type=struct
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-08-25      Written by Matthew Argall
%
function srvy = mms_edi_srvy_combine(slow_srvy, fast_srvy)

	% Is there data (both must be defined)
	tf_slow = ~isempty(slow_srvy);
	tf_fast = ~isempty(fast_srvy);

	% No data
	if ~tf_fast && ~tf_slow
		srvy = [];

	% No fast data
	elseif ~tf_fast
		srvy = slow_srvy;
		
	% No slow data
	elseif ~tf_slow
		srvy = fast_srvy;
	
	% Combine fast and slow
	else
		% Create the combined structure
		srvy = struct();
		
		% Sort and combine times.
		[srvy.('tt2000_gd12'), isort_gd12] = sort( [slow_srvy.tt2000_gd12 fast_srvy.tt2000_gd12] );
		[srvy.('tt2000_gd21'), isort_gd21] = sort( [slow_srvy.tt2000_gd21 fast_srvy.tt2000_gd21] );
		
		% Find unique times
		[~, iuniq_gd12] = unique( srvy.tt2000_gd12 );
		[~, iuniq_gd21] = unique( srvy.tt2000_gd21 );
		
		% Sorted, unique order
		iorder_gd12 = isort_gd12( iuniq_gd12 );
		iorder_gd21 = isort_gd21( iuniq_gd21 );
		
		% Remove the time fields from the slow and fast structures
		rmfield(slow_srvy, { 'tt2000_gd12', 'tt2000_gd21'} );
		rmfield(edi_fast, { 'tt2000_gd12', 'tt2000_gd21'} );
	
		% Get the names of each remaining field.
		fields  = fieldnames( slow_srvy );
		nFields = length( fields );
		
		% Find fields related to gd12
		igd12 = find( ~cellfun( @isempty, regexp( fields, '(gd12|gun1)', 'once' ) ) );
		igd21 = find( ~cellfun( @isempty, regexp( fields, '(gd21|gun2)', 'once' ) ) );

		% Step through each field.
		%   - Both guns have the same fields.
		%   - Each structure has the same fields.
		for ii = 1 : nFields / 2
			%  Combine data; select sorted, unique elements; add to EDI structure
			field_gd12       = fields{ igd12(ii) };
			tmp_data         = [ slow_srvy.( field_gd12 ) fast_srvy.( field_gd12 ) ];

			% Gun positions are scalar in the spinning frame.
			if ~strcmp(field_gd12, 'virtual_gun1_123')
				tmp_data = tmp_data(:, iorder_gd12);
			end
			srvy.(field_gd12) = tmp_data;
			
			% Repeat for GD21
			field_gd21       = fields{ igd21(ii) };
			tmp_data         = [ slow_srvy.( field_gd21 ) fast_srvy.( field_gd21 ) ];
			if ~strcmp(field_gd21, 'virtual_gun2_123')
				tmp_data = tmp_data(:, iorder_gd21);
			end
			srvy.(field_gd21) = tmp_data;
		end
	end
end