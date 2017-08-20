%
% Name
%   mms_fsm_test_savgol
%
% Purpose
%   
%
% MATLAB release(s) MATLAB 7.14.0.739 (R2012a)
% Required Products None
%
% History:
%   2016-09-09      Written by Matthew Argall
%
%***************************************************************************

% Data file to plot
% file    = '/nfs/fsm/temp/mms4_fsm_srvy_l2plus_cal-dfg_20151101153733_v0.0.1.cdf';
% file    = '/nfs/fsm/temp/mms4_fsm_srvy_l2plus_cal-scm_20151103160519_v0.1.0.cdf';
% file    = '/nfs/fsm/temp/mms4_fsm_brst_l2plus_cal-scm-week_20151101000000_v0.0.0.cdf';
file    = '/nfs/fsm/temp/mms4_fsm_brst_l2plus_cal-dfg-month_20151101000000_v0.0.0.cdf';
theComp = 'z';
theFlag = 1;
theFreq = 1.0;
fc      = 0.5;             % Corner frequency of high-pass filter

%------------------------------------%
% Variable Names                     %
%------------------------------------%

% Parse file name to get variable name prefix and suffix
[sc, instr, mode, level, tstart, version, optdesc] = mms_dissect_filename(file);
prefix    = [sc '_' instr '_'];
suffix    = ['_' mode '_' level];
mag_instr = optdesc(5:7);

% Create Variable Names
amp_vname         = [prefix 'amp'   '_omb'   suffix];
phase_vname       = [prefix 'phase' '_omb'   suffix];
psd_vname         = [prefix 'psd'   '_omb'   suffix];
amp_hist_vname    = [prefix 'amp'   '_hist'  suffix];
phase_hist_vname  = [prefix 'phase' '_hist'  suffix];
psd_hist_vname    = [prefix 'psd'   '_hist'  suffix];
amp_floor_vname   = [prefix 'amp'   '_floor' suffix];
phase_floor_vname = [prefix 'phase' '_floor' suffix];
psd_floor_vname   = [prefix 'psd'   '_floor' suffix];

% MAT file
matfile = [prefix suffix(2:end) '_cal-floor-' optdesc(9:end) '_' tstart '_v' version '.mat'];

clear sc instr mode level optdesc prefix suffix % mag_instr

%------------------------------------%
% Read Data                          %
%------------------------------------%

% Spectra
% [amp, t, f] = MrCDF_Read(file, amp_vname,   'RowMajor', true);
% phase       = MrCDF_Read(file, phase_vname, 'RowMajor', true);
% [psd, t, f] = MrCDF_Read(file, psd_vname,   'RowMajor', true);

% Histogram
% [amp_hist,   ~, amp_bins, flag, comp] = MrCDF_Read(file, amp_hist_vname,   'RowMajor', true);
% [phase_hist, ~, phase_bins]           = MrCDF_Read(file, phase_hist_vname, 'RowMajor', true);
[psd_hist, f, psd_bins, flag, comp]     = MrCDF_Read(file, psd_hist_vname,   'RowMajor', true);

% Noise floor
% amp_floor   = MrCDF_Read(file, amp_floor_vname);
% phase_floor = MrCDF_Read(file, phase_floor_vname);
psd_floor   = MrCDF_Read(file, psd_floor_vname);

%------------------------------------%
% Pick Indices                       %
%------------------------------------%

% Component
switch theComp
	case 'x'
		iComp = 1;
	case 'y'
		iComp = 2;
	case 'z'
		iComp = 3;
	otherwise
		error( ['Invalid component: "' theComp '".'] );
end

% Flag
iFlag = find(flag == theFlag);
assert( ~isempty(iFlag), 'Invalid flag provided.' );

% Frequency
if0 = find(f >= theFreq, 1, 'first');

%------------------------------------%
% Create 2D Coordinate Grids         %
%------------------------------------%
%nTime     = length(t);
nFreq  = length(f);
nFlags = length(flag);
nBins  = length(psd_bins);

%psd      = squeeze(psd(:,iComp,:));
psd_hist  = squeeze(psd_hist(:,iFlag,iComp,:));
psd_floor = squeeze(psd_floor(:,iFlag,iComp));

clear iComp iFlag

%------------------------------------%
% Smooth Distribution                %
%------------------------------------%

% Probability distribution
sumdist = cumsum( double( psd_hist ), 1 );
sdist   = sumdist ./ repmat( sum( double( psd_hist ), 1 ), nBins, 1 );

% Smoothed step response
nCol = 5;
nRow = 11;
sdist_smooth = [ repmat(sdist(:,1),        1, (nCol-1)/2), sdist,        repmat(sdist(:,end),        1, (nCol-1)/2) ];
sdist_smooth = [ repmat(sdist_smooth(1,:), (nRow-1)/2, 1); sdist_smooth; repmat(sdist_smooth(end,:), (nRow-1)/2, 1) ];
kernel1      = ones(1,nRow) / nRow;
kernel2      = ones(1,nCol) / nCol;
sdist_smooth = conv2(kernel1, kernel2, sdist_smooth, 'valid');

clear nCol nRow kernel1 kernel2

%------------------------------------%
% Fit TanH                           %
%------------------------------------%

% Initial conditions
p = [0.5, 0.5, -0.5];

% Allocate memory
hfit = zeros(3, nFreq);

%
% FORWARD
%

% Fit for each frequency
for ii = if0 : nFreq

	% Create a new function handle
	%   - Must do so to update the Y-data
	x     = psd_bins;
	y     = sdist_smooth(:,ii)';
	ftanh = @(p) sum( ( y - ( 0.5*tanh(p(1)*x - p(2)) + p(3) ) ).^2 );

	% Get the fit parameters
	temp = fminsearch( ftanh, p );
	
	% Save results
	hfit(:,ii) = temp;
	
	% Next iteration
	p = temp;
end

%
% BACKWARD
%

% Initial conditions
p = [0.5, 0.5, -0.5];

ifc = find(f >= fc, 1, 'first');
for ii = if0-1 : -1 : ifc

	% Create a new function handle
	%   - Must do so to update the Y-data
	x     = double( psd_hist(:,ii)' );
	y     = cumsum(x) / sum(x);
	ftanh = @(p) sum( ( y - ( 0.5*tanh(p(1)*x - p(2)) + p(3) ) ).^2 );

	% Get the fit parameters
	temp = fminsearch( ftanh, p );
	
	% Save results
	hfit(:,ii) = temp(2);
	
	% Next iteration
	p = temp;
end

clear p temp x y ftanh

%------------------------------------%
% Resulting Distribution             %
%------------------------------------%

% Smoothed distribution
pdist_smooth = (sdist_smooth - [zeros(1,nFreq); sdist_smooth(1:end-1,:)]) .* repmat(sum(psd_hist, 1), nBins, 1);

% Fitted Results
p1        = repmat(hfit(1,:), nBins, 1);
p2        = repmat(hfit(2,:), nBins, 1);
p3        = repmat(hfit(3,:), nBins, 1);
sdist_fit = 0.5*tanh( p1 .* repmat(psd_bins', 1, nFreq) - p2 ) + p3;
pdist_fit = (sdist_fit - [zeros(1,nFreq); sdist_fit(1:end-1,:)]) .* repmat(sum(psd_hist, 1), nBins, 1);

clear p1 p2 p3


% The result
%   - The point of 50% probability at each frequency
[~,imin] = min( abs( sdist_fit - 0.5 ), [], 1 );
result   = psd_bins(imin);

%------------------------------------%
% Save as MAT File                   %
%------------------------------------%


if true
	% Output file name
	matfile = fullfile('/home/argall/', matfile);
	
	% Data structure
	data = struct( [ mag_instr '_f'           ], f,            ...
	               [ mag_instr '_bins'        ], psd_bins,     ...
	               [ mag_instr '_dist'        ], psd_hist,     ...
	               [ mag_instr '_dist_smooth' ], pdist_smooth, ...
	               [ mag_instr '_dist_fit'    ], pdist_fit,    ...
	               [ mag_instr '_fit'         ], hfit,         ...
	               [ mag_instr '_step'        ], sdist,        ...
	               [ mag_instr '_step_smooth' ], sdist_smooth, ...
	               [ mag_instr '_step_fit'    ], sdist_fit,    ...
	               [ mag_instr '_floor'       ], result        );
	
	% Append/overwrite data from an existing file
	if exist(matfile, 'file') == 2
		save( matfile, '-struct', 'data', '-append' );
	else
		save( matfile, '-struct', 'data' );
	end
	clear data
end


%------------------------------------%
% Scale Data                         %
%------------------------------------%

% Scale the PSD
%sclPSD = MrRescale( log10( psd ), 1, 64, ...
%                    'MinValue', -6, ...
%                    'MaxValue', -1, ...
%                    'Class',    'uint8' );

%
% Histogram
%
sclDist = MrRescale( psd_hist, 1, 64, ...
                     'Class', 'uint8' );

sclDist_smooth = MrRescale( pdist_smooth, 1, 64, ...
                            'Class', 'uint8' );

pdist_fit( pdist_fit > 1.2*max(pdist_smooth(:))  ) = 0;
pdist_fit( pdist_fit < 0 ) = 0;
sclDist_fit = MrRescale( pdist_fit, 1, 64, ...
                         'Class', 'uint8' );

%
% Step response
%
sclStep_hist = MrRescale( sdist, 1, 64, ...
                          'Class', 'uint8' );

sclStep_smooth = MrRescale( sdist_smooth, 1, 64, ...
                            'Class', 'uint8' );

sdist_fit( sdist_fit > 1 ) = 0;
sdist_fit( sdist_fit < 0 ) = 0;
sclStep_fit = MrRescale( sdist_fit, 1, 64, ...
                         'Class', 'uint8' );

%------------------------------------%
% Plot Distributions                 %
%------------------------------------%
% Create a figure
fig = figure();

% Original Histogram
subplot(3,3,1);
h = pcolor( f, psd_bins, sclDist );
% h = pcolor(f_hist_grid, bin_hist_grid, sclPSD_hist);
h.EdgeColor = 'none';
% colorbar
hold on
h = plot(f, result, '--m');
%line([theFreq, theFreq], [psd_bins(1) psd_bins(end)], 'LineStyle', '--', 'Color', [0,0,0])
hold off


% Smoothed Histogram
subplot(3,2,3);
h = pcolor( f, psd_bins, sclDist_smooth );
h.EdgeColor = 'none';
hold on
h = plot(f, result, '--m');
hold off


% Fitted Histogram
subplot(3,2,5);
h = pcolor( f, psd_bins, sclDist_fit );
h.EdgeColor = 'none';
hold on
h = plot(f, result, '--m');
hold off


% Original Step Response
subplot(3,2,2);
h = pcolor( f, psd_bins, sclStep_hist );
h.EdgeColor = 'none';
hold on
h = plot(f, result, '--m');
hold off


% Smoothed Step Response
subplot(3,2,4);
h = pcolor( f, psd_bins, sclStep_smooth );
h.EdgeColor = 'none';
hold on
h = plot(f, result, '--m');
hold off


% Fitted Step Response
subplot(3,2,6);
h = pcolor( f, psd_bins, sclStep_fit );
h.EdgeColor = 'none';
hold on
h = plot(f, result, '--m');
hold off

%------------------------------------%
% Plot A Cut                         %
%------------------------------------%

fig2 = figure();

% Distribution
subplot(2,1,1);
plot( psd_bins, psd_hist(:,if0), psd_bins, pdist_smooth(:,if0), psd_bins, pdist_fit(:,if0) );
legend('Signal', 'Smooth', 'TanH');


% Probability distribution
subplot(2,1,2);
plot( psd_bins, sdist(:,if0), psd_bins, sdist_smooth(:,if0), psd_bins, sdist_fit(:,if0) );
legend('Signal', 'Smooth', 'TanH');
ylim([-0.5,1.5])
