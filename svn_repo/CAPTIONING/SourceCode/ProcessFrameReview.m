function [OutImage] = ProcessFrameReview(InImage,K);

%Declare Global variables
global intStartFrm;
global intSkipFrms;
global FigureHandle;
global Cnt;
global CurrentK;
global FramesLeft;
global FramesBuffer;
global FrameNumbersBuffer;
global ActualBufferSize;
global BufferSize;
global InBufferIndex;
global StayOnCurrent;
global CurrentIndex;
global InPause;
global DimX;
global DimY;
global DimZ;
global BlankEventsPanel;
global BlankMenuPanel;
global BlankInfoPanel;
global PanelDimX;
global EventsPanelDimY;
global MenuPanelDimY;
global InfoPanelDimY;
global Spacer;
global ToBeUpdated;
global PanelSize;
global ReSizeFigure;
global CurrentEvents;
global NCurrEvents;
global MaxCurrentEvents;
global CaptionFileName;
global Events;
global NEvents;
global ESCDetected;
global HImage;
global HMenu;
global HInfo;
global HEvents;



CurrentK = K;
ToBeUpdated.Image = 1;
ToBeUpdated.Info = 1;
ToBeUpdated.Events = 1;

CurrFrame = K-intStartFrm+1;
JustBuffering = (K < intSkipFrms);
if CurrFrame == 1 | ReSizeFigure,
    ReSizeFigure = 0;
    [OldDimY, OldDimX, OldDimZ] = size(InImage);
    DimZ = OldDimZ;
    DimY = 480;
    DimX = round(OldDimX * DimY / OldDimY);
    Spacer =  0.02;
    PanelDimX = 160;

    % Resize the Figure
    rect = get(FigureHandle,'Position');
    rect(1) = 10;
    rect(4) = DimY + (4 * DimY * Spacer);
    rect(3) = DimX + PanelDimX + (6 * DimY * Spacer);
    set(FigureHandle,'Position',rect);

    
    NormVec = rect([3,4,3,4]);

%    PanelSize.Image.Left = (2 * DimY * Spacer);
%    PanelSize.Image.Low = (2 * DimY * Spacer);
%    PanelSize.Image.Height = DimY;
%    PanelSize.Image.Width = DimX;
%    PanelSize.Image.PVec = [PanelSize.Image.Left, PanelSize.Image.Low, PanelSize.Image.Width, PanelSize.Image.Height]./NormVec;

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
    BlankEventsPanel = ones(round(EventsPanelDimY),PanelDimX,DimZ);
    BlankMenuPanel = ones(round(MenuPanelDimY),PanelDimX,DimZ);
    BlankInfoPanel = ones(round(InfoPanelDimY),PanelDimX,DimZ);
    HImage = subplot('position',PanelSize.Image.PVec);
    set(HImage,'Parent', FigureHandle,'XTick',[],'YTick',[]);
    HEvents = subplot('position',PanelSize.Events.PVec);
    set(HEvents, 'Parent', FigureHandle,'XTick',[],'YTick',[]);
    HMenu = subplot('position',PanelSize.Menu.PVec);
    set(HMenu, 'Parent', FigureHandle,'XTick',[],'YTick',[]);
    HInfo = subplot('position',PanelSize.Info.PVec);
    set(HInfo, 'Parent', FigureHandle,'XTick',[],'YTick',[]);
end
figure(FigureHandle);
set(FigureHandle,'KeyPressFcn',@KeyboardHandling_Review);
set(FigureHandle,'DeleteFcn',@CloseImageHandling);
% Prepare Output Image
% GrayInImage = rgb2gray(InImage); 
% OutImage = zeros(size(InImage));
% OutImage(:,:,1) = GrayInImage;
% OutImage(:,:,2) = GrayInImage;
% OutImage(:,:,3) = GrayInImage;
OutImage = InImage;



if ESCDetected & ~JustBuffering,
    OutImage = [];
elseif StayOnCurrent & ~JustBuffering,
    while StayOnCurrent,
        PlotImage_Review(FramesBuffer{CurrentIndex},FigureHandle,FrameNumbersBuffer(CurrentIndex));
        drawnow;
    end
elseif InPause & ~JustBuffering,
    while InPause,
        ToBeUpdated.Image = 0;
        ToBeUpdated.Info = 0;
        PlotImage_Review(InImage,FigureHandle,K);
        drawnow;
        pause(0.5);
    end
else
    if ActualBufferSize < BufferSize,
        ActualBufferSize = ActualBufferSize + 1;
    end
    FramesBuffer{InBufferIndex} = InImage;
    FrameNumbersBuffer(InBufferIndex) = K;
    
        
    if InBufferIndex < ActualBufferSize,
        InBufferIndex = InBufferIndex + 1;
    else
        InBufferIndex = 1;
    end
    
    if ~JustBuffering,
        if FramesLeft > 1,
            FramesLeft = FramesLeft - 1;
        else
            FramesLeft = Cnt;
            PlotImage_Review(OutImage,FigureHandle,K);
        end
    end
end
if K == 30,
    uuu = 0;
end
if ~JustBuffering,
%     CurrentEvents(K+1+1,:) = CurrentEvents(K+1,:);
%     NCurrEvents(K+1+1) = NCurrEvents(K+1);
end

