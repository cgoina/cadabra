%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Declare Global variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



global intStartFrm;  
global NumPressed;
global Boxing;
global intSkipFrms;
global intNFrms;
global FigureHandle;
global MaxCnt;
global Cnt;
global FramesLeft;
global FramesBuffer;
global FrameNumbersBuffer;
global BufferSize;
global ActualBufferSize;
global InBufferIndex;
global StayOnCurrent;
global CurrentIndex;
global Rewinded;
global InPause;
global DimX;
global DimY;
global DimZ;
global Spacer;
global BlankEventsPanel;
global BlankInfoPanel;
global PanelDimX;
global EventsPanelDimY;
global InfoPanelDimY;
global ToBeUpdated;
global Events;
global NEvents;
global PanelSize;
global ReSizeFigure;
global ClosingEvent;
global CurrentEvents;
global MaxCurrentEvents;
global NCurrEvents;
global UniqueID;
global CaptionFileName;
global ESCDetected;
global BackupCaptionFID;
global BackupCaptionFileName;
global CaptionFileName;
global HImage;
global HMenu;
global HInfo;
global HEvents;
global slashstr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%strInVideoPath = ['/scratch'];
strInVideoPath = ['/common/fanti/PietrosOffice'];
%strOutVideoPath = ['/common/fanti/PietrosOffice'];
strOutVideoPath = [''];
strInCapPath = ['/common/fanti/PietrosOffice'];
strOutCapPath = ['/common/fanti/PietrosOffice'];
%strInCapPath = ['/scratch'];
%strOutCapPath = ['/scratch'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Video File Names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

strVideoFExt = ['avi'];
strInVideoFName = ['PietrosOfficeDownSampled'];
%strOutVideoFName = [strInVideoFName '_OUT'];
strOutVideoFName = [''];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Caption File Names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

strCapFExt = ['mat'];
strInCapFName = strInVideoFName;
strOutCapFName = strInCapFName;

CaptionFileName = [ strOutCapPath slashstr strOutCapFName '.' strCapFExt];
BackupCaptionFileName = [ strOutCapPath '/Backup_' strOutCapFName '.txt'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings for ProcessAVI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

strInFName = [strInVideoPath slashstr strInVideoFName '.' strVideoFExt];
strOutFName = [strOutVideoPath slashstr strOutVideoFName '.' strVideoFExt];
strProcFrmFcn = ['ProcessFrameReview'];          % Use this function to do the processing
bolMovOrFrame = 0;                                      % Get an AVI Movie out
bolMakeGray = 0;                                        % Leave RGB as is.
intStartFrm = 0;                                        % Start from this Frame
intNFrms = 100000;                                     % Number of Frames to Process 
intQuality = 10000;                                     % Highest Quality
intFps = 29.97;                                         % 30 Frame per Second


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check Existence of Input/Output Files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ResumeFromCaption = 0;
ResumeFromBackup = 0;
if ~exist(strInFName,'file'),
    disp(['ERROR: The Following Video-In file was not found']);
    disp(['NAME: ', strInVideoFName '.' strVideoFExt]);
    disp(['PATH: ', strInVideoPath slashstr]);
    return
end
if exist(strOutFName,'file') & bolMovOrFrame == 1,
    disp(['ERROR: The Following Video-Out file already exists']);
    disp(['         NAME: ', strOutVideoFName '.' strVideoFExt]);
    disp(['         PATH: ', strOutVideoPath slashstr]);
    return
end
CapExist = exist(CaptionFileName,'file');
if ~CapExist,
    disp(['ERROR: The Following Caption file DOES NOT already']);
    disp(['         NAME: ', strOutCapFName '.' strCapFExt]);
    disp(['         PATH: ', strOutCapPath slashstr]);
    return
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Events List
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Events = struct('Description',cell(NEvents,1));

NEvents = 9;
Events(1).Description = 'Cart Moving';
Events(2).Description = 'Person Walking';
Events(3).Description = 'Person Running';
Events(4).Description = 'Person Cycling';
Events(5).Description = 'Person Skating';
Events(6).Description = 'Person Other';
Events(7).Description = 'Insect Moving';
Events(8).Description = 'Car Moving';
Events(9).Description = 'Other';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remember to Initialize the Display Panels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ToBeUpdated.Info = 1;
ToBeUpdated.Image = 1;
ToBeUpdated.Menu = 1;
ToBeUpdated.Events = 1;
ReSizeFigure = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start Processing from a Known Coherent Status
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ClosingEvent = 0;
StayOnCurrent = 0;
Boxing = 0;
ESCDetected = 0;
InPause = 0;
Rewinded = 1;
CurrentIndex = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resume
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LoadCaptionFile;
MaxCurrentEvents = 9;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare a Buffer to FF/REW on already viewed frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BufferSize = 3;
ActualBufferSize = 1; % Buffer is smaller till we have seen "BufferSize" Frames
InBufferIndex = 1;
FramesBuffer = cell(BufferSize,1);
FrameNumbersBuffer = zeros(BufferSize,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Max and Initial Speed of PLAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cnt = 1;
MaxCnt = 16;
FramesLeft = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the Display and assign Handling Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FigureHandle = figure;
set(FigureHandle,'KeyPressFcn',@KeyboardHandling_Review);
set(FigureHandle,'DeleteFcn',@CloseImageHandling);
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call ProcessAVI to start processing Frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To Mex the File Use the following command :
%
% mex ProcessAVI.cc  -I/usr/include/avifile-0.7  /usr/lib/avifile-0.7/win32.so /usr/lib/avifile-0.7/xvid4.so -L/usr/lib -ljpeg -laviplay;


ProcessAVI(strInFName, strOutFName, strProcFrmFcn, bolMovOrFrame, ...
    bolMakeGray, intStartFrm, intNFrms, intQuality, intFps)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Close the Figure created
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close(FigureHandle);

disp(['Done.']);
