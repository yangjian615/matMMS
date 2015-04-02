%
% Construct an MMS file name. The file name format is:
%
%   scId_instrumentId_mode_dataLevel_optionalDataProductDescriptor_startTime_vX.Y.Z.cdf
%
% :Examples:
%   Dissect and construct a filename
%     varname = 'mms2_dfg_lo_l2pre_G';
%     [sc, instr, name, desc] = mms_dissect_filename(varname);
%     vname = mms_construct_varname(sc, instr, name, desc)
% 				vname =
% 				mms2_dfg_lo_l2pre_G
%
% :Params:
%   SC:                 in, required, type=char
%                       Spacecraft identifier:
%                           mms1
%                           mms2
%                           mms3
%                           mms4
%   INSTR:              in, required, type=char
%                       Instrument or investigation identifier
%                           hpca
%                           aspoc
%                           epd
%                           epd-eis
%                           epd-feeps
%                           fpi
%                           des
%                           dis
%                           des-dis
%                           fields
%                           edi
%                           adp
%                           sdp
%                           adp-sdp
%                           afg
%                           dfg
%                           dsp
%                           afg-dfg
%                           scm
%   NAME:               in, required, type=string
%                       Parameter name of the variable.
%   DESC:               in, optional, type=char, default=''
%                       Optional data product descriptor. Should be short
%                           (3-8 characters). Underscores used to separate
%                           multiple components.
%
function fname = mms_construct_varname(sc, instr, name, desc)

%------------------------------------%
% Check Inputs                       %
%------------------------------------%
	% Default description
	if nargin < 4
		desc = '';
	end
	
	% Separate the description from the parameter name
	if ~isempty(desc)
		desc = ['_' desc];
	end
	 
%------------------------------------%
% Create Variable Name               %
%------------------------------------%

	% Construct the file name.
	fname = [sc '_' instr '_' name desc];
end