%
% Name
%   mms_fsm_fgm_read
%
% Purpose
%   Read FGM data.
%
% Calling Sequence
%   DATA = mms_fsm_fgm_read( L1A_FILES, TRANGE )
%     Read FGM L1A files with names given by L1A_FILES. L2 magnetic field data in OMB
%     coordinates is obtained by passing L1A data into Ken Brommund's L2 data processing
%     routines in IDL. Data is then sorted into categories and returned in the structure
%     DATA. Categories are: 1) slow/fast srvy, 2) hi/lo-range 3) pre/post perigee
%     4) DEC 32/64.
%
%   DATA = mms_fsm_fgm_read( ..., TRANGE )
%     If survey data is being read, the TRANGE specifies the time range of a single
%     orbit, starting at the beginning of slow survey and ending at the end of the
%     following fast survey period. Times should be formatted as 'yyyy-mm-ddTHH:MM:SS'
%     or as empty strings.
%
%   DATA = mms_fsm_fgm_read(__, FC)
%     Give the bass-cut frequency at which to high-pass filter the data.
%
% Parameters
%   L1A_FILES       in, required, type = char/cell
%   L2PRE_FILES     in, required, type = char/cell
%   TRANGE          in, required, type = 1x2 cell
%   FC              in, optional, type = double, default = []
%
% Returns
%   DATA            out, required, type=struct
%                   Fields are:
%                     'sc'       - Spacecraft ID
%                     'instr'    - Instrument that recorded the data
%                     'mode'     - Data telemetry mode
%                     'optdesc'  - Optional file name descriptor
%                     't'        - TT2000 epoch times
%                     'b'        - Magnetic field in DMPA coordinates
%                     'flag'     - Bit flag indicating operational mode:
%                                    1: unset = Slow      set = Fast/Brst
%                                    2: unset = Lo-Range  set = Hi-Range
%                                    3: unset = Deck64    set = Deck32     (dfg)
%                                    3: unset = ADCA      set = ADCB       (afg)
%                                    4: unset = ---       set = Pre-Perigee
%                                    5: unset = ---       set = Post-Perigee
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2015-06-01      Written by Matthew Argall
%   2015-10-01      Renamed from mms_fsm_bkgd_fgm_read to mms_fsm_fgm_read. FC now
%                     defaults to the empty array. - MRA
%   2017-01-30      Add 'sc', 'instr', 'mode', and 'optdesc' fields to the output. - MRA
%
function fgm = mms_fsm_fgm_read( l1a_files, trange, fc )

	if nargin() < 4
		fc = [];
	end

%------------------------------------%
% Read Data                          %
%------------------------------------%

	% Read the L1A files
	[pfmode, t_packet] = mms_fsm_fgm_read_l1a(l1a_files, trange);

%------------------------------------%
% Create L2 Temporary Data           %
%------------------------------------%
	% Get the log file
	oLog    = mrstdlog();
	logfile = oLog.filename;
	tf_log  = isempty( regexp(logfile, '^(std|log)') );

	% l2_files
	if ischar(l1a_files)
		[sc, instr, mode, ~, tstart] = mms_dissect_filename(l1a_files);
	else
		[sc, instr, mode, ~, tstart] = mms_dissect_filename(l1a_files{1});
	end
	
	%
	% Create FGM L2 data in OMB
	%
	cmd = [ fileparts(mfilename('fullpath')) '/mms_fsm_fgm_l2.sh ' sc ' ' instr ' ' mode ' ' tstart];
	if tf_log
		cmd = [cmd ' ' logfile];
	end
	
	mrfprintf( 'logtext', '\n\n' );
	mrfprintf( 'logtext', '=====================================' );
	mrfprintf( 'logtext', '| Calling IDL to Create L2 Data     |' );
	mrfprintf( 'logtext', '=====================================' );
	mrfprintf( 'logtext', '%s', cmd );
	mrfprintf( 'logtext', '\n\n' );

	% Call IDL
	[status, cmdout] = system(cmd);

	% Parse results
	%   - Both the the log file and the data file are returned.
	l2_file = regexp(cmdout, '\n', 'split');
	il2     = find( ~cellfun(@isempty, l2_file), 1, 'last' );
	l2_file = l2_file{il2};
	
	% system will capture output from both stdout and stderr
	%   - If no log file is defined, redirect to stdout
	if ~tf_log
		mrfprintf( 'stdout', '%s', cmdout );
	end
	
	% Check status
	assert( status <= 100, 'Error processing L2 data.' );
	
	mrfprintf( 'logtext',  '\n\n' );
	mrfprintf( 'logtext', '========================================' );
	mrfprintf( 'logtext', '| Returning to MATLAB to Process FSM   |' );
	mrfprintf( 'logtext', '========================================' );
	mrfprintf( 'logtext',  '\n\n' );
	
	% Read the L2Pre files
	[b_omb, hirange, rate, t] = mms_fsm_fgm_read_l2temp(l2_file, trange);
	
	% Delete the temporary L2 file
	delete(l2_file);

%------------------------------------%
% Categorize Data                    %
%------------------------------------%
	
	% Categorize the data
	fgm = mms_fsm_fgm_read_categorize( t_packet, t, b_omb, rate, hirange, pfmode );
	
	% Clear all of the data
	clear t_packet t b_omb rate range pfmode

%------------------------------------%
% Isolate Modes                      %
%------------------------------------%
	
	% Flag
	%   - Bit 1 set:     Fast and Brst
	%   - Bit 1 not set: Slow 
	%   - iFast serves as iBrst when iSlow is empty
	iFast = find( bitget(fgm.flag, 1, 'uint8') );
	iSlow = find( bitget(fgm.flag, 1, 'uint8') == 0 );

%------------------------------------%
% High-Pass Filter                   %
%------------------------------------%
	if ~isempty(fc)
		% Slow
		if ~isempty(iSlow)
			% Create filter
			%   - filefilt requires doubles
			sr = double( fgm.sr( iSlow(1) ) );
			fN = sr / 2.0;
			d  = fdesign.highpass('N,F3dB', 9, fc/fN);
			h  = design(d, 'butter');

			% Apply filter
			fgm.b(1,iSlow) = filtfilt( h.sosMatrix, h.ScaleValues, double( fgm.b(1,iSlow) ));
			fgm.b(2,iSlow) = filtfilt( h.sosMatrix, h.ScaleValues, double( fgm.b(2,iSlow) ));
			fgm.b(3,iSlow) = filtfilt( h.sosMatrix, h.ScaleValues, double( fgm.b(3,iSlow) ));
		end

		% Fast (or burst)
		if ~isempty(iFast)
			% Create filter
			%   - butter() + filtfilt() was producing NaNs. This was the solution
			%   - https://www.mathworks.com/matlabcentral/answers/62553-filtfilt-function-returning-nan-at-certain-frequencies
			sr = double( fgm.sr( iFast(1) ) );
			fN = sr / 2.0;
			d  = fdesign.highpass('N,F3dB', 9, fc/fN);
			h  = design(d, 'butter');
		
			% Apply the filter
			fgm.b(1,iFast) = filtfilt( h.sosMatrix, h.ScaleValues, double( fgm.b(1,iFast) ));
			fgm.b(2,iFast) = filtfilt( h.sosMatrix, h.ScaleValues, double( fgm.b(2,iFast) ));
			fgm.b(3,iFast) = filtfilt( h.sosMatrix, h.ScaleValues, double( fgm.b(3,iFast) ));
		end
	end

%------------------------------------%
% Metadata                           %
%------------------------------------%
	
	% Add to data structure
	fgm.sc      = sc;
	fgm.instr   = instr;
	fgm.mode    = mode;
	fgm.optdesc = '';
end