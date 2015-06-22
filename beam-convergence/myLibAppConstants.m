% myLibAppConstants % custom colors; set default axis colors
% http://www.rapidtables.com/web/color/RGB_Color.htm
% http://www.rapidtables.com/convert/color/rgb-to-hex.htm

myBlack         = [   0.0,   0.0,   0.0 ]        ; %
myBlue          = [   0.0, 114.0, 178.0 ] / 255.0; %
myDarkBlue      = [  23.0,  42.0, 135.0 ] / 255.0; %
myDarkGreen     = [   6.0, 150.0,   6.0 ] / 255.0; %
myDarkRed       = [ 237.0,  14.0,  51.0 ] / 255.0; %
myDarkTeal      = [   6.0, 150.0, 194.0 ] / 255.0; %
myFawn          = [ 229.0, 170.0, 112.0 ] / 255.0; %
myGrassGreen    = [   8.0, 209.0,  22.0 ] / 255.0; %
myLightFawn     = [ 248.0, 230.0, 212.0 ] / 255.0; %

myLightGrey1    = [  32.0,  32.0,  32.0 ] / 255.0; %
myLightGrey2    = [  64.0,  64.0,  64.0 ] / 255.0; %
myLightGrey3    = [  96.0,  96.0,  96.0 ] / 255.0; %
myLightGrey4    = [ 128.0, 128.0, 128.0 ] / 255.0; %
myLightGrey5    = [ 160.0, 160.0, 160.0 ] / 255.0; %
myLightGrey6    = [ 192.0, 192.0, 192.0 ] / 255.0; %

myGold          = [ 212.0, 175.0,  55.0 ] / 255.0; %
myGoldPowder    = [ 213.0,  94.0,   0.0 ] / 255.0; % D55E00
myLightTan      = [ 240.0, 232.0, 220.0 ] / 255.0; %
myLightWheat    = [ 255.0, 232.0, 200.0 ] / 255.0; %
myLightWheat    = [ 250.0, 240.0, 220.0 ] / 255.0; %
myOrange        = [ 247.0, 163.0,   7.0 ] / 255.0; %
myTan           = [ 210.0, 180.0, 140.0 ] / 255.0; %
myWheat         = [ 245.0, 222.0, 179.0 ] / 255.0; %

% color-blind compliant palette
% http://jfly.iam.u-tokyo.ac.jp/color/
myBluishGreen   = [   0.0, 158.0, 115.0 ] / 255.0; %
myReddishPurple = [ 204.0, 121.0, 167.0 ] / 255.0; %
mySkyBlue       = [  86.0, 180.0, 233.0 ] / 255.0; %
myVermilion     = [ 227.0,  66.0,  52.0 ] / 255.0; % E34234
myMedVermilion  = [ 217.0,  96.0,  59.0 ] / 255.0; %

% MMS-specific colors
MMS1_obsColor   = myBlack;
MMS2_obsColor   = myVermilion;
MMS3_obsColor   = myBluishGreen;
MMS4_obsColor   = mySkyBlue;
MMS_plotColorx  = mySkyBlue;
MMS_plotColory  = myBluishGreen;
MMS_plotColorz  = myVermilion;

DefaultAxesColors = [ myDarkBlue; myDarkRed; myDarkGreen; myOrange; myGrassGreen; myDarkTeal ];
DefaultAxesColors = [ MMS_plotColorx; MMS_plotColory; MMS_plotColorz; myOrange; myGrassGreen; myDarkTeal ];
set (0, 'DefaultAxesColorOrder', DefaultAxesColors);
% -~-~-~-~-~-~-~-~-~

ScreenSize      = get (0, 'ScreenSize'); % [ left bottom width height ]
ScreenSizeCells = num2cell (ScreenSize (1:4));
[ ScreenLeft, ScreenBottom, ScreenWidth, ScreenHeight ] = deal (ScreenSizeCells {:});

useAltDisplay = false;
monitorPositions = get (0, 'MonitorPositions');
% DisplayOffset (1:4) = monitorPositions (1, 1:4); % Xmin Ymin Xmax Ymax
DisplayOffset (1:4) = 0.0; % Xmin Ymin Xmax Ymax

if useAltDisplay
	monitorPositionsSize = size (monitorPositions);
	if (monitorPositionsSize (1) > 1)
		% sample monitorPositions = [[ 1 1 1920 1080 ]; [ 1921 1 3840 1200 ]]
		% The 2nd row is the extension of the display onto the 2nd display
		% These are now the offsets to the new display; at the simplest level,
		% these can simply be added to Display1 values.
		DisplayOffset (1:4) = monitorPositions (2, 1:4) - monitorPositions (1, 1:4) % extended display Xmin Ymin Xmax Ymax
	end
end
% -~-~-~-~-~-~-~-~-~

cPathSep = pathsep;
cFileSep = filesep;

EDI_presentation_beam_plot_style.Version         = '1';
EDI_presentation_beam_plot_style.Format          = 'png';
EDI_presentation_beam_plot_style.Preview         = 'none';
EDI_presentation_beam_plot_style.Width           = 'auto';
EDI_presentation_beam_plot_style.Height          = 'auto';
EDI_presentation_beam_plot_style.Units           = 'inches';
EDI_presentation_beam_plot_style.Color           = 'rgb';
EDI_presentation_beam_plot_style.Background      = 'w';
EDI_presentation_beam_plot_style.FixedFontSize   = '10';
EDI_presentation_beam_plot_style.ScaledFontSize  = 'auto';
EDI_presentation_beam_plot_style.FontMode        = 'scaled';
EDI_presentation_beam_plot_style.FontSizeMin     = '8';
EDI_presentation_beam_plot_style.FixedLineWidth  = '1';
EDI_presentation_beam_plot_style.ScaledLineWidth = 'auto';
EDI_presentation_beam_plot_style.LineMode        = 'fixed';
EDI_presentation_beam_plot_style.LineWidthMin    = '1';
EDI_presentation_beam_plot_style.FontName        = 'auto';
EDI_presentation_beam_plot_style.FontWeight      = 'auto';
EDI_presentation_beam_plot_style.FontAngle       = 'auto';
EDI_presentation_beam_plot_style.FontEncoding    = 'latin1';
EDI_presentation_beam_plot_style.PSLevel         = '2';
EDI_presentation_beam_plot_style.Renderer        = 'painters';
EDI_presentation_beam_plot_style.Resolution      = '300';
EDI_presentation_beam_plot_style.LineStyleMap    = 'none';
EDI_presentation_beam_plot_style.ApplyStyle      = '0';
EDI_presentation_beam_plot_style.Bounds          = 'loose';
EDI_presentation_beam_plot_style.LockAxes        = 'on';
EDI_presentation_beam_plot_style.LockAxesTicks   = 'on';
EDI_presentation_beam_plot_style.ShowUI          = 'on';
EDI_presentation_beam_plot_style.SeparateText    = 'off';
