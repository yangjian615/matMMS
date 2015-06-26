% edi_drift_step_read_ql__EDI__B__EDP_data

UseFileOpenGUI = true;

% ~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~
if UseFileOpenGUI
	[mms_ql__EDI__B__dataFile, mms_ql_dataPath] = uigetfile ('mms*.cdf', 'Select an MMS ql EDI_&_B CDF file');
	if isequal (mms_ql__EDI__B__dataFile,  0) % then no valid file selected
		msgbox ('No valid MMS ql EDI_&_B data file selected.');
	else
		mms_ql__EDI__B__data = [mms_ql_dataPath, mms_ql__EDI__B__dataFile];
	end
end

%{
				 0         0         0         0         0
mms2_edi_slow_ql_efield_20150509_v0.1.2.cdf
v0.1.2	
  'Epoch'                             [1x2 double]    [ 501]    'tt2000'    'T/'     'Full'    'None'    [0]    [-9223372036854775808]
  'Epoch_delta_plus'                  [1x2 double]    [   1]    'tt2000'    'F/'     'Full'    'None'    [0]    [-9223372036854775808]
  'epoch_gd12_beam'                   [1x2 double]    [3488]    'tt2000'    'T/'     'Full'    'None'    [0]    [-9223372036854775808]
  'epoch_gd21_beam'                   [1x2 double]    [3142]    'tt2000'    'T/'     'Full'    'None'    [0]    [-9223372036854775808]
  'mms2_edi_E_dmpa'                   [1x2 double]    [ 501]    'single'    'T/T'    'Full'    'None'    [0]    [      -1.0000000e+30]
  'mms2_edi_v_ExB_dmpa'               [1x2 double]    [ 501]    'single'    'T/T'    'Full'    'None'    [0]    [      -1.0000000e+30]
  'mms2_edi_d_dmpa'                   [1x2 double]    [ 501]    'single'    'T/T'    'Full'    'None'    [0]    [      -1.0000000e+30]
  'mms2_edi_B_dmpa'                   [1x2 double]    [ 501]    'single'    'T/T'    'Full'    'None'    [0]    [      -1.0000000e+30]
  'mms2_edi_pos_virtual_gun1_dmpa'    [1x2 double]    [3488]    'single'    'T/T'    'Full'    'None'    [0]    [      -1.0000000e+30]
  'mms2_edi_pos_virtual_gun2_dmpa'    [1x2 double]    [3142]    'single'    'T/T'    'Full'    'None'    [0]    [      -1.0000000e+30]
  'mms2_edi_fv_gd12_dmpa'             [1x2 double]    [3488]    'single'    'T/T'    'Full'    'None'    [0]    [      -1.0000000e+30]
  'mms2_edi_fv_gd21_dmpa'             [1x2 double]    [3142]    'single'    'T/T'    'Full'    'None'    [0]    [      -1.0000000e+30]
  'mms2_edi_recnum'                   [1x2 double]    [ 501]    'int32'     'T/'     'Full'    'None'    [0]    [         -2147483647]
  'mms2_edi_recnum_gd12'              [1x2 double]    [3488]    'int32'     'T/'     'Full'    'None'    [0]    [         -2147483647]
  'mms2_edi_recnum_gd21'              [1x2 double]    [3142]    'int32'     'T/'     'Full'    'None'    [0]    [         -2147483647]
  'mms2_edi_beam_quality_gd12'        [1x2 double]    [3488]    'int8'      'T/'     'Full'    'None'    [0]    [                -127]
  'mms2_edi_beam_quality_gd21'        [1x2 double]    [3142]    'int8'      'T/'     'Full'    'None'    [0]    [                -127]
  'mms2_edi_d_std_dmpa'               [1x2 double]    [ 501]    'single'    'T/T'    'Full'    'None'    [0]    [      -1.0000000e+30]
  'mms2_edi_B_std_dmpa'               [1x2 double]    [ 501]    'single'    'T/T'    'Full'    'None'    [0]    [      -1.0000000e+30]
  'E_Labl_Ptr'                        [1x2 double]    [   3]    'char'      'T/'     'Full'    'None'    [0]    '  '                  
  'v_Labl_Ptr'                        [1x2 double]    [   3]    'char'      'T/'     'Full'    'None'    [0]    '  '                  
  'd_Labl_Ptr'                        [1x2 double]    [   3]    'char'      'T/'     'Full'    'None'    [0]    '  '                  
  'B_Labl_Ptr'                        [1x2 double]    [   3]    'char'      'T/'     'Full'    'None'    [0]    '  '                  
  'vg_labl_vname'                     [1x2 double]    [   3]    'char'      'T/'     'Full'    'None'    [0]    ' '                   

				 0         0         0         0         0
mms2_edp_comm_ql_dce2d_20150509000000_v0.0.0.cdf'
	'mms2_edp_dce_epoch'      [1x2 double]    [2764672]    'tt2000'    'T/'     'Full'    'None'      [    0]    [-9223372036854775808]
	'LABL_1'                  [1x2 double]    [      1]    'char'      'F/T'    'Full'    'GZIP.6'    [    1]    '     '               
	'mms2_edp_dce_xyz_dsl'    [1x2 double]    [2764672]    'single'    'T/T'    'Full'    'GZIP.6'    [ 5462]    [      -1.0000000e+30]
	'mms2_edp_dce_bitmask'    [1x2 double]    [2764672]    'uint8'     'T/'     'Full'    'GZIP.6'    [65536]    [                 254]
	'mms2_edp_dce_quality'    [1x2 double]    [2764672]    'int16'     'T/'     'Full'    'GZIP.6'    [32768]    [              -32767]
%}

obsID = mms_ql__EDI__B__dataFile (4:4);
YYYY  = str2num (mms_ql__EDI__B__dataFile (30:33));
MM    = str2num (mms_ql__EDI__B__dataFile (34:35));
DD    = str2num (mms_ql__EDI__B__dataFile (36:37));

% ~~~~~~~~~~~~~~~~~~~
if ~isequal (mms_ql__EDI__B__dataFile, 0) % then a valid [we hope] file selected
	disp ([ 'Reading MMS EDI_&_B data... ', mms_ql__EDI__B__data ])
	mms_ql_EDI_dataFile_info = spdfcdfinfo (mms_ql__EDI__B__data);

	edi_BdvE_tt2000 = spdfcdfread (mms_ql__EDI__B__data, ...
		'CombineRecords',        true, ...
		'Variable',              'Epoch', ...
		'ConvertEpochToDatenum', false, ...
		'KeepEpochAsIs',         true);
	edi_B_dmpa = spdfcdfread (mms_ql__EDI__B__data, ...
		'CombineRecords',        true, ...
		'Variable',              ['mms', obsID, '_edi_B_dmpa']);
	edi_BdvE_recnum = spdfcdfread (mms_ql__EDI__B__data, ...
		'CombineRecords',        true, ...
		'Variable',              ['mms', obsID, '_edi_recnum']);
	disp 'Date range of edi_BdvE_tt2000'
	[ datestr(spdftt2000todatenum(edi_BdvE_tt2000(1)),   'yyyy-mm-dd HH:MM:ss'), ' ',...
	  datestr(spdftt2000todatenum(edi_BdvE_tt2000(end)), 'yyyy-mm-dd HH:MM:ss') ]


	edi_gd12_beam_tt2000 = spdfcdfread (mms_ql__EDI__B__data, ...
		'CombineRecords',        true, ...
		'Variable',              'epoch_gd12_beam', ...
		'ConvertEpochToDatenum', false, ...
		'KeepEpochAsIs',         true);
	edi_gd12_virtual_dmpa = spdfcdfread (mms_ql__EDI__B__data, ...
		'CombineRecords',        true, ...
		'Variable',              ['mms', obsID, '_edi_pos_virtual_gun1_dmpa']);
	edi_gd12_fv_dmpa = spdfcdfread (mms_ql__EDI__B__data, ...
		'CombineRecords',        true, ...
		'Variable',              ['mms', obsID, '_edi_fv_gd12_dmpa']);
	edi_gd12_B_index = spdfcdfread (mms_ql__EDI__B__data, ...
		'CombineRecords',        true, ...
		'Variable',              ['mms', obsID, '_edi_recnum_gd12']);
	edi_gd12_quality = spdfcdfread (mms_ql__EDI__B__data, ...
		'CombineRecords',        true, ...
		'Variable',              ['mms', obsID, '_edi_beam_quality_gd12']);
	disp 'Date range of edi_gd12_beam_tt2000'
	[ datestr(spdftt2000todatenum(edi_gd12_beam_tt2000(1)),   'yyyy-mm-dd HH:MM:ss'), ' ',...
	  datestr(spdftt2000todatenum(edi_gd12_beam_tt2000(end)), 'yyyy-mm-dd HH:MM:ss') ]

	edi_gd21_beam_tt2000 = spdfcdfread (mms_ql__EDI__B__data, ...
		'CombineRecords',        true, ...
		'Variable',              'epoch_gd21_beam', ...
		'ConvertEpochToDatenum', false, ...
		'KeepEpochAsIs',         true);
	edi_gd21_virtual_dmpa = spdfcdfread (mms_ql__EDI__B__data, ...
		'CombineRecords',        true, ...
		'Variable',              ['mms', obsID, '_edi_pos_virtual_gun2_dmpa']);
	edi_gd21_fv_dmpa = spdfcdfread (mms_ql__EDI__B__data, ...
		'CombineRecords',        true, ...
		'Variable',              ['mms', obsID, '_edi_fv_gd21_dmpa']);
	edi_gd21_B_index = spdfcdfread (mms_ql__EDI__B__data, ...
		'CombineRecords',        true, ...
		'Variable',              ['mms', obsID, '_edi_recnum_gd21']);
	edi_gd21_quality = spdfcdfread (mms_ql__EDI__B__data, ...
		'CombineRecords',        true, ...
		'Variable',              ['mms', obsID, '_edi_beam_quality_gd21']);
	disp 'Date range of edi_gd21_beam_tt2000'
	[ datestr(spdftt2000todatenum(edi_gd21_beam_tt2000(1)),   'yyyy-mm-dd HH:MM:ss'), ' ',...
	  datestr(spdftt2000todatenum(edi_gd21_beam_tt2000(end)), 'yyyy-mm-dd HH:MM:ss') ]

	% ~~~~~~~~~~~~~~~~~~~
	% mms2_edp_comm_ql_dce2d_20150509000000_v0.0.0.cdf'
	% ~~~~~~~~~~~~~~~~~~~
	if UseFileOpenGUI
		[mms_ql__EDP_dataFile, mms_ql_dataPath] = uigetfile ('mms*.cdf', 'Select an MMS ql EDP CDF file');
		if isequal (mms_ql__EDP_dataFile,  0) % then no valid file selected
			msgbox ('No valid MMS ql EDP data file selected.');
		else
			mms_ql__EDP_data = [mms_ql_dataPath, mms_ql__EDP_dataFile];
		end
	end
	
	if ~isequal (mms_ql__EDP_dataFile, 0) % then a valid [we hope] file selected
		disp ([ 'Reading EDP data... ', mms_ql__EDP_data ])
		mms_ql_EDP_dataFile_info = spdfcdfinfo (mms_ql__EDP_data);
	
		edp_tt2000 = spdfcdfread (mms_ql__EDP_data, ...
			'CombineRecords',        true, ...
			'Variable',              'mms2_edp_dce_epoch', ...
			'ConvertEpochToDatenum', false, ...
			'KeepEpochAsIs',         true);
		% Electric field in DSL coordinates (DSL ~= DMPA)
		edp_dce_xyz_dsl = spdfcdfread (mms_ql__EDP_data, ...
			'CombineRecords',        true, ...
			'Variable',              ['mms', obsID, '_edp_dce_xyz_dsl']);

		dceFillVal = -1.0000000e+30; % inconsisten fill values
		idceFillVal = find (edp_dce_xyz_dsl (:, 1) < -1.0000000e+20);
		edp_dce_xyz_dsl (idceFillVal, :) = [];
		edp_tt2000      (idceFillVal)    = [];
		edp_datenum = spdftt2000todatenum (edp_tt2000);
% 		disp 'Date range of edp_tt2000'
% 		[ datestr(spdftt2000todatenum(edp_tt2000(1)),   'yyyy-mm-dd HH:MM:ss'), ' ',...
% 		  datestr(spdftt2000todatenum(edp_tt2000(end)), 'yyyy-mm-dd HH:MM:ss') ]
	end

	% If CDF_FileInfo.FileSettings.Majority = 'Row' then transpose to col vector matrix.
	% But MATLAB searches faster in cols, so if there is filtering to be done,
	% better to do it along cols, and transpose only at when necessary for math.

	% ~~~~~~~~~~~~~~~~~~~
	iBeqNaN = find (isnan (edi_B_dmpa (:,1)));
	if ~isempty (iBeqNaN)
		iBeqNaN
	end
	edi_B_dmpa_FillVal = -1.0000000e+30;
	iBeqFillVal = find (edi_B_dmpa == edi_B_dmpa_FillVal);
	iBeqBad = union (iBeqNaN, iBeqFillVal);
	edi_BdvE_tt2000 (iBeqBad   ) = [];
	edi_B_dmpa      (iBeqBad, :) = [];

	% There may be NaNs in the B data
	% If any 5 s period has no B data, we can't use those BEAMs as CENTER beams,
	% but we could use them as beams on either side, so the logic gets more complicated...
	% If we remove the bad B data, we must remove the EDI data that corresponds.
	% If we keep the bad B records, then we must check B valid for each EDI record set.
	% The problem with deleting the bad B data records is that we mess up the
	% EDI-record:B-record index is corrupted, because there is NO B-record index ---
	% the EDI index cross references are based on B-record position. Therefore,
	% if we delete a B-record, all the following records get bumped down by one,
	% and this kills the EDI:B record cross-reference.

	clear edi_gd12_quality
	clear edi_gd21_quality

	% Now the problem is that we want to scroll thru the firing vectors (beams)
	% in linear time, choosing the center beam, and finding 4-8 beams on either side
	% FROM EITHER GDU.
	% That means concat time, fv, and B_index. Q is no longer necessary.
	% Then create index of time, sorted in ascending order. Scroll thru the data
	% during analysis, using the sorted time index.

	% Now is the time to change from nx3 data to 3xn.
	edi_gd_beam_tt2000  = [ edi_gd12_beam_tt2000'     edi_gd21_beam_tt2000' ];
	edi_gd_virtual_dmpa = [ edi_gd12_virtual_dmpa'    edi_gd21_virtual_dmpa' ];
	edi_gd_fv_dmpa      = [ double(edi_gd12_fv_dmpa') double(edi_gd21_fv_dmpa') ];
	edi_gd_B_xref      = [ edi_gd12_B_index'         edi_gd21_B_index' ];
	% Better keep track of which GDU fired the beam. Corresponds to concatenation order above.
	edi_gd_ID           = [ zeros(1,length(edi_gd12_fv_dmpa),'uint8')+1, zeros(1,length(edi_gd21_fv_dmpa),'uint8')+2 ];
	% ... and tranpose B records
	B_datenum       = spdftt2000todatenum (edi_BdvE_tt2000);
	edi_BdvE_tt2000 = edi_BdvE_tt2000';
	edi_B_dmpa      = edi_B_dmpa';

	[ ~, iSorted_beam_tt2000 ] = sort (edi_gd_beam_tt2000, 2);
	% 	[ edi_gd_beam_tt2000(iSorted_beam_tt2000(1:20))', ...
	% 	  edi_gd_ID(iSorted_beam_tt2000(1:20))'        ]
	ValidDataLoaded = true; % assume that user knows data is good

	% keep EDp data that is in the range of EDI data
	iEDP_lt_EDI = find (edp_datenum < B_datenum (1));
	iEDP_gt_EDI = find (edp_datenum > B_datenum (end));
	iEDP_EDI_noMatch = union (iEDP_lt_EDI, iEDP_gt_EDI);
	edp_datenum (iEDP_EDI_noMatch) = [];
	edp_dce_xyz_dsl (iEDP_EDI_noMatch, :) = [];
	disp 'Date range of edp_tt2000'
	[ datestr(spdftt2000todatenum(edp_tt2000(1)),   'yyyy-mm-dd HH:MM:ss'), ' ',...
	  datestr(spdftt2000todatenum(edp_tt2000(end)), 'yyyy-mm-dd HH:MM:ss') ]
end % ~isequal (mms_ql__EDI__B__dataFile, 0)
