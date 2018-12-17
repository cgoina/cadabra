function ProcessFrameCapWin(data,width,height,frameNr,time)

%Declare Global variables
global strInFName;
global intStartFrm;
global intNFrms;
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
global QuitPrg;
global filechg;
global width1;
global height1;

persistent warned;

if nargin,
    if ESCDetected || QuitPrg || filechg, error('STOP PROCESSING'); end
try
    width1 = width; height1= height;
%     scanline = ceil(width*3/4)*4; % the scanline size must be a multiple of 4.
% 
%     % some times the amount of data doesn't match exactly what we expect...
%     if (numel(data) ~= scanline*height)
%         if (numel(data) > 3*width*height)
%             if (isemtpy(warned))
%                 warning('mmread:general','dimensions do not match data size. Guessing badly...');
%                 warned = true;
%             end
%             scanline = width*3;
%             data = data(1:3*width*height);
%         else
%             error('dimensions do not match data size. Too little data.');
%         end
%     end
    % if there is any extra scanline data, remove it
%     data = reshape(data,scanline,height);
%     data = data(1:3*width,:);

    % the data ordering is wrong for matlab images, so permute it
%     data = permute(reshape(data, 3, width, height),[3 2 1]);
    % the images are also upside down and colors were backwards.
%     InImage = double(flipdim(flipdim(data,1),3))/255;

%     InImage = double(reshape(data, 3, width, height))/255;
    
    CurrentK = frameNr;
%         %     if (intNFrms > abs(movie_arr.nrFramesTotal)), intNFrms = abs(movie_arr.nrFramesTotal); end;
    dt = time / CurrentK;

    ToBeUpdated.Image = 1;
    ToBeUpdated.Info = 0;

    JustBuffering = (CurrentK < intSkipFrms);
    if (CurrentK == intStartFrm) || ReSizeFigure,
        ReSizeFigure = 0;
%         [OldDimY, OldDimX, OldDimZ] = size(InImage);
        DimZ = 3; %OldDimZ;
        DimY = 480;
        DimX = 640; %round(OldDimX * DimY / OldDimY);
        Spacer =  0.02;
        PanelDimX = 160;

        % Resize the Figure
        rect = get(FigureHandle,'Position');
        rect(1) = 100;
        rect(2) = 100;
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
    end
    figure(FigureHandle);
    set(FigureHandle,'KeyPressFcn',@KeyboardHandling_Captioning);
    set(FigureHandle,'DeleteFcn',@CloseImageHandling);
    % Prepare Output Image
    % GrayInImage = rgb2gray(InImage);
    % OutImage = zeros(size(InImage));
    % OutImage(:,:,1) = GrayInImage;
    % OutImage(:,:,2) = GrayInImage;
    % OutImage(:,:,3) = GrayInImage;
%     OutImage = InImage;

    if StayOnCurrent && ~JustBuffering,
        while StayOnCurrent,
            ToBeUpdated.Info = 1;
            PlotImage_Captioning(FramesBuffer{CurrentIndex},FigureHandle,FrameNumbersBuffer(CurrentIndex));
            drawnow;
        end
    elseif InPause && ~JustBuffering,
        while InPause,
            if QuitPrg || filechg,
                error('STOP PROCESSING');
            end
            ToBeUpdated.Image = 0;
            ToBeUpdated.Info = 0;
            PlotImage_Captioning(data,FigureHandle,CurrentK);
            drawnow;
            pause(0.5);
        end
    else
        if ActualBufferSize < BufferSize,
            ActualBufferSize = ActualBufferSize + 1;
        end
        FramesBuffer{InBufferIndex} = data;
        FrameNumbersBuffer(InBufferIndex) = CurrentK;


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
                PlotImage_Captioning(data,FigureHandle,CurrentK);
            end
        end
    end
    if ~JustBuffering,
        CurrentEvents(CurrentK+1+1,:) = CurrentEvents(CurrentK+1,:);
        NCurrEvents(CurrentK+1+1) = NCurrEvents(CurrentK+1);
    end
catch
    % we will never see the error unless we catch and display it ourselves
    disp(lasterr);
end
end


function OpenAVI
global strVideoFExt;
global strInVideoFName;
global strInVideoPath;
global findex;
global filechg;
global InPause;

[strInVideoFName, strInVideoPath, findex] = uigetfile({'*.avi','AVI-file (*.avi)'});
if findex, 
    strInVideoPath = strInVideoPath(1:end-1);
    strVideoFExt = strInVideoFName(end-2:end);
    strInVideoFName = strInVideoFName(1:end-4);
    InPause = 0;
    filechg = 1; 
end
return

function exitprg
global QuitPrg;

QuitPrg = 1;
return
