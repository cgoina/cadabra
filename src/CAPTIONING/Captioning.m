function Captioning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Declare Global variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global strVideoFExt;
global strInVideoFName;
global strInVideoPath;
global strInFName;
global findex;
global findex1;
global intStartFrm;  
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
global HImage;
global HMenu;
global HInfo;
global HEvents;
global firstfile;
global QuitPrg;
global filechg;
global slashstr;

firstfile = 0;
QuitPrg = 0;
ESCDetected = 0;

if (ispc), 
    slashstr = '\';
else
    slashstr = '/';
end

while ~QuitPrg && ~ESCDetected
findex = 0;
findex1 = 0;
filechg = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Events List
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NEvents = 9;
Events = struct('Description',cell(NEvents,1));
Events(1).Description = 'Left Wing';
Events(2).Description = 'Right Wing';
Events(3).Description = 'Both Wings';
Events(4).Description = 'Lunging';
Events(5).Description = 'Chasing';
Events(6).Description = 'Wing Threat';
Events(7).Description = 'Tussling';
Events(8).Description = 'Copulation';
Events(9).Description = 'Circling';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remember to Initialize the Display Panels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ToBeUpdated.Info = 0;
ToBeUpdated.Image = 0;
ToBeUpdated.Menu = 0;
ToBeUpdated.Events = 0;
ReSizeFigure = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start Processing from a Known Coherent Status
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ClosingEvent = 0;
StayOnCurrent = 0;
Boxing = 0;
InPause = 0;
Rewinded = 1;
CurrentIndex = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare a Buffer to FF/REW on already viewed frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BufferSize = 100;
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
%set(0,'DefaultFigureProperty',PropertyValue...)
%set(FigureHandle,'KeyPressFcn',@KeyboardHandling_Captioning);
%set(FigureHandle,'DeleteFcn',@CloseImageHandling_Captioning);
%drawnow;

    OldDimY = 576;
    OldDimX = 720;
    OldDimZ = 3;
    DimZ = OldDimZ;
    DimY = 480;
    DimX = round(OldDimX * DimY / OldDimY);
    Spacer =  0.02;
    PanelDimX = 160;

    % Resize the Figure
    rect = get(FigureHandle,'Position');
    if (ispc),
        rect(1) = 100;
        rect(2) = 100;
    else
        rect(1) = 20;
        rect(2) = 20;
    end
    rect(4) = DimY + (4 * DimY * Spacer);
    rect(3) = DimX + PanelDimX + (6 * DimY * Spacer);
    set(FigureHandle,'Position',rect);

    NormVec = rect([3,4,3,4]);

    PanelSize.Image.Left = (2 * DimY * Spacer);
    PanelSize.Image.Low = (2 * DimY * Spacer);
    PanelSize.Image.Height = DimY;
    PanelSize.Image.Width = DimX;
    PanelSize.Image.PVec = [PanelSize.Image.Left, PanelSize.Image.Low, PanelSize.Image.Width, PanelSize.Image.Height]./NormVec;

    PanelSize.Info.Left = (4 * DimY * Spacer) + DimX;
    PanelSize.Info.Low = (6 * DimY * Spacer) + (DimY - 4 * DimY * Spacer) * 14 / 15;
    PanelSize.Info.Height = (DimY - 2 * DimY * Spacer) * 1 / 15;
    PanelSize.Info.Width = PanelDimX;
    PanelSize.Info.PVec = [PanelSize.Info.Left, PanelSize.Info.Low, PanelSize.Info.Width, PanelSize.Info.Height]./NormVec;

    PanelSize.Menu.Left = (4 * DimY * Spacer) + DimX;
    PanelSize.Menu.Low = (4 * DimY * Spacer) + (DimY - 4 * DimY * Spacer) * 5 / 15;
    PanelSize.Menu.Height = (DimY - 4 * DimY * Spacer) * 9 / 15;
    PanelSize.Menu.Width = PanelDimX;
    PanelSize.Menu.PVec = [PanelSize.Menu.Left, PanelSize.Menu.Low, PanelSize.Menu.Width, PanelSize.Menu.Height]./NormVec;

    PanelSize.Events.Left = (4 * DimY * Spacer) + DimX;
    PanelSize.Events.Low = (2 * DimY * Spacer);
    PanelSize.Events.Height = (DimY - 4 * DimY * Spacer) * 5 / 15;
    PanelSize.Events.Width = PanelDimX;
    PanelSize.Events.PVec = [PanelSize.Events.Left, PanelSize.Events.Low, PanelSize.Events.Width, PanelSize.Events.Height]./NormVec;

    EventsPanelDimY = PanelSize.Events.Height;
    MenuPanelDimY = PanelSize.Menu.Height;
    InfoPanelDimY = PanelSize.Info.Height;
    BlankEventsPanel = ones(round(EventsPanelDimY),PanelDimX,DimZ)/4;
    BlankMenuPanel = ones(round(MenuPanelDimY),PanelDimX,DimZ)/4;
    BlankInfoPanel = ones(round(InfoPanelDimY),PanelDimX,DimZ)/4;
    HImage = subplot('position',PanelSize.Image.PVec);
    set(HImage,'Parent', FigureHandle,'XTick',[],'YTick',[]);
    HEvents = subplot('position',PanelSize.Events.PVec);
    set(HEvents, 'Parent', FigureHandle,'XTick',[],'YTick',[]);
    HMenu = subplot('position',PanelSize.Menu.PVec);
    set(HMenu, 'Parent', FigureHandle,'XTick',[],'YTick',[]);
    HInfo = subplot('position',PanelSize.Info.PVec);
    set(HInfo, 'Parent', FigureHandle,'XTick',[],'YTick',[]);

figure(FigureHandle);

PlotImage_Captioning(FramesBuffer{1},FigureHandle,FrameNumbersBuffer(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

strOutVideoPath = strInVideoPath;
strInCapPath = strInVideoPath;
strOutCapPath = strInVideoPath;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Video File Names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

strOutVideoFName = strInVideoFName;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Caption File Names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

strCapFExt = ['mat'];
strInCapFName = strInVideoFName;
strOutCapFName = strInCapFName;

CaptionFileName = [ strOutCapPath slashstr strOutCapFName '.' strCapFExt];
BackupCaptionFileName = [ strOutCapPath slashstr 'Backup_' strOutCapFName '.txt'];
strInFName = [strInVideoPath slashstr strInVideoFName '.' strVideoFExt];
strOutFName = [strOutVideoPath slashstr strOutVideoFName '.' strVideoFExt];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings for ProcessAVI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

strProcFrmFcn = ['ProcessFrameCaptioning'];          % Use this function to do the processing
bolMovOrFrame = 0;                                      % Get an AVI Movie out
bolMakeGray = 0;                                        % Leave RGB as is.
intStartFrm = 1;                                        % Start from this Frame
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
BackupCapExist = exist(BackupCaptionFileName,'file');
if CapExist & BackupCapExist,
    disp(['ERROR: The Following Caption and Backup files already exist']);
    disp(['         NAME: ', strOutCapFName '.' strCapFExt]);
    disp(['         PATH: ', strOutCapPath slashstr]);
    disp(['         NAME: Backup_', strOutCapFName '.txt']);
    disp(['         PATH: ', strOutCapPath slashstr]);
    ButtonName = [];
    while isempty(ButtonName),
        ButtonName = questdlg('Both the Caption and Backup Files Exist - What would you like to do ?', ...
                           'Inconsistent State Detected', ...
                           'Delete Backup','Delete Caption','Exit','Delete Caption');
    end
                   
    switch ButtonName,
        case 'Delete Backup',
            delete(BackupCaptionFileName);
            disp('WARNING: The Caption file has been deleted.');
        case 'Delete Caption',
            delete(CaptionFileName);
            disp('WARNING: The Backup file has been deleted.');
        case 'Exit',
            disp('WARNING: The user chose to EXIT.');
            return
    end
end
CapExist = exist(CaptionFileName,'file');
BackupCapExist = exist(BackupCaptionFileName,'file');

if CapExist,
    disp(['ERROR: The Following Caption file already exists']);
    disp(['         NAME: ', strOutCapFName '.' strCapFExt]);
    disp(['         PATH: ', strOutCapPath slashstr]);
    ButtonName = [];
    while isempty(ButtonName),
        ButtonName = questdlg('The Caption File Exists - What would you like to do ?', ...
                           'Partial Results Exist', ...
                           'Resume','Delete','Exit','Resume');
    end
    switch ButtonName,
        case 'Resume',
            disp('WARNING: Will RESUME working of the aforementioned file.');
            ResumeFromCaption = 1;
        case 'Delete',
            delete(CaptionFileName);
            disp('WARNING: The file has been deleted.');
        case 'Exit',
            disp('WARNING: The user chose to EXIT.');
            return
    end
elseif BackupCapExist,
    disp(['ERROR: The Following Backup file already exists']);
    disp(['         NAME: Backup_', strOutCapFName '.txt']);
    disp(['         PATH: ', strOutCapPath slashstr]);
    ButtonName = [];
    while isempty(ButtonName),
        ButtonName = questdlg('The Backup File Exists - What would you like to do ?', ...
                           'Partial Results Exist', ...
                           'Resume','Delete','Exit','Resume');
    end
    switch ButtonName,
        case 'Resume',
            disp('WARNING: Will RESUME working on the aforementioned file.');
            ResumeFromBackup = 1;
        case 'Delete',
            delete(BackupCaptionFileName);
            disp('WARNING: The file has been deleted.');
        case 'Exit',
            disp('WARNING: The user chose to EXIT.');
            return
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resume Options/Start anew
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ResumeFromBackup,
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Current Frame's Events List
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SaveFromBackupToCaptionFile;
    SaveFromCaptionToBackupFile;
    disp(['WARNING: The Following Backup file has been Converted to a Caption File and Deleted']);
    disp(['         NAME: Backup_', strOutCapFName '.txt']);
    disp(['         PATH: ', strOutCapPath slashstr]);
    disp(['WARNINIG: The Following Caption file has been Backed Up and Deleted']);
    disp(['         NAME: ', strOutCapFName '.' strCapFExt]);
    disp(['         PATH: ', strOutCapPath slashstr]);
    disp(['WARNINIG: The Following Backup file replaces it and is being used']);
    disp(['         NAME: Backup_', strOutCapFName '.txt']);
    disp(['         PATH: ', strOutCapPath slashstr]);
    delete(CaptionFileName);
    MaxCurrentEvents = 9;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Prepare Backup FILE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Note that we APPEND - not OVERWRITE
    [BackupCaptionFID] = fopen(BackupCaptionFileName,'a+');
    if BackupCaptionFID == -1,
        disp(['ERROR: The Following Backup file could not be Opened for Append']);
        disp(['         NAME: Backup_', strOutCapFName '.txt']);
        disp(['         PATH: ', strOutCapPath slashstr]);
        return
    end
    
elseif ResumeFromCaption

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Current Frame's Events List
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SaveFromCaptionToBackupFile;
    delete(CaptionFileName);
    disp(['WARNINIG: The Following Caption file has been Backed Up and Deleted']);
    disp(['         NAME: ', strOutCapFName '.' strCapFExt]);
    disp(['         PATH: ', strOutCapPath slashstr]);
    disp(['WARNINIG: The Following Backup file replaces it and is being used']);
    disp(['         NAME: Backup_', strOutCapFName '.txt']);
    disp(['         PATH: ', strOutCapPath slashstr]);
    MaxCurrentEvents = 9;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Prepare Backup FILE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Note that we APPEND - not OVERWRITE
    [BackupCaptionFID] = fopen(BackupCaptionFileName,'a+');
    if BackupCaptionFID == -1,
        disp(['ERROR: The Following Backup file could not be Opened for Append']);
        disp(['         NAME: Backup_', strOutCapFName '.txt']);
        disp(['         PATH: ', strOutCapPath slashstr]);
        return
    end

else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Current Frame's Events List
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    MaxCurrentEvents = 9;
    CurrentEvents = ...
        struct('StartFrame',cell(intNFrms - intStartFrm,MaxCurrentEvents),...
        'EventCode',cell(intNFrms - intStartFrm,MaxCurrentEvents),...
        'UniqueID',cell(intNFrms - intStartFrm,MaxCurrentEvents));
    NCurrEvents = zeros(intNFrms - intStartFrm,1);
    UniqueID = 0;
    intSkipFrms = 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Prepare Backup FILE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [BackupCaptionFID] = fopen(BackupCaptionFileName,'w+');
    if BackupCaptionFID == -1,
        disp(['ERROR: The Following Backup file could not be Created']);
        disp(['         NAME: Backup_', strOutCapFName '.txt']);
        disp(['         PATH: ', strOutCapPath slashstr]);
        return
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remember to Initialize the Display Panels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ToBeUpdated.Info = 1;
ToBeUpdated.Image = 1;
ToBeUpdated.Menu = 1;
ToBeUpdated.Events = 1;
ReSizeFigure = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Call ProcessAVI to start processing Frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To Mex the File Use the following command :
%
% mex ProcessAVI.cc  -I/usr/include/avifile-0.7  /usr/lib/avifile-0.7/win32.so /usr/lib/avifile-0.7/xvid4.so -L/usr/lib -ljpeg -laviplay;

if (ispc),
%     ProcessFrameCaptioning_win(intStartFrm);
    ProcessFrameCapWin;
    try
        mmread(strInFName,intStartFrm:intStartFrm+intNFrms,[],false,false,'ProcessFrameCapWin');
    catch
        err = lasterror;
        if ~strcmp(err.message(end-14:end),'STOP PROCESSING')
            rethrow(err); % if it was any other type of error, pass it on so that we can see it
        end
    end    
else
    if (strVideoFExt == 'avi') || (strVideoFExt == 'wmv'),
        ProcessAVI(strInFName, strOutFName, strProcFrmFcn, bolMovOrFrame, ...
            bolMakeGray, intStartFrm, intNFrms, intQuality, intFps);
    else
        ProcessFrameCaptioning_fmf(intStartFrm);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Close the Figure created
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close(FigureHandle);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save Files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['Saving Files. Please wait...']);
fclose(BackupCaptionFID);
SaveFromBackupToCaptionFile;
delete(BackupCaptionFileName);
disp(['Done.']);
QuitPrg = 1;

end
