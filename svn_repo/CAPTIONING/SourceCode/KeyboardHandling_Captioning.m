function KeyboardHandling_Captioning(src,evnt)
global Cnt;
global MaxCnt;
global FramesLeft;
global StayOnCurrent;
global CurrentIndex;
global FrameNumbersBuffer;
global FramesBuffer;
global ActualBufferSize;
global InBufferIndex;
global Rewinded;
global InPause;
global ToBeUpdated;
global NEvents;
global ClosingEvent;
global NCurrEvents;
global MaxCurrentEvents;
global CurrentEvents;
global CurrentK;
global UniqueID;
global ESCDetected;
global BackupCaptionFID;
global HImage;
global HMenu;
global HInfo;
global HEvents;
global EventNumber;
global Boxing;
global FigureHandle;

if StayOnCurrent,
    NumPressed = str2num(evnt.Character);
    IsNum = 1;
    if isempty(NumPressed),
        IsNum = 0;
    end
    if ClosingEvent,
        CF = FrameNumbersBuffer(CurrentIndex);
        if (evnt.Character == 30),
            ToBeUpdated.Menu = 1;
            ToBeUpdated.Events = 1;
            ClosingEvent = 0;
        elseif IsNum & ismember(NumPressed,[1:NCurrEvents(CF)]),
            % Close that event
            EndFrame = FrameNumbersBuffer(CurrentIndex);
            CF = EndFrame;
            ThisUniqueID = CurrentEvents(CF+1,NumPressed).UniqueID;
            ThisCode = CurrentEvents(CF+1,NumPressed).EventCode;
            % Select Image Area
            ToBeUpdated.Menu = 1;
            Boxing = 1;
            EventNumber = CurrentEvents(EndFrame+1,NumPressed).EventCode;
            PlotImage_Captioning(FramesBuffer{CurrentIndex},FigureHandle,FrameNumbersBuffer(CurrentIndex));
%             Boxing = 0;
%             rect = [];
%             while isempty(rect),
%                 rect = getrect(HImage);
%             end
%             rect = fix(rect);
            while  CF <= CurrentK,
                NCurrEvents(CF+1) = NCurrEvents(CF+1) - 1;
                Gap = 0;
                for i = 1:NCurrEvents(CF+1), 
                    if CurrentEvents(CF+1,i).UniqueID == ThisUniqueID,
                        Gap = 1;
                    end
                    CurrentEvents(CF+1,i) = CurrentEvents(CF+1,i+Gap);
                end
                CF = CF + 1;
            end
%             fprintf(BackupCaptionFID,'%010d E %010d %03d %d %d %d %d\n',ThisUniqueID, EndFrame, ThisCode,rect(1),rect(2),rect(3),rect(4));
            fprintf(BackupCaptionFID,'%010d E %010d %03d\n',ThisUniqueID, EndFrame, ThisCode);            ToBeUpdated.Menu = 1;
            ToBeUpdated.Events = 1;
            ClosingEvent = 0;
        else
            beep;
        end
    else
        if  (evnt.Character == 28) & (Rewinded < ActualBufferSize - 1),
            Rewinded = Rewinded + 1;
            if CurrentIndex > 1,
                CurrentIndex = CurrentIndex - 1;
            else
                CurrentIndex = ActualBufferSize;
            end
        elseif  (evnt.Character == 29) & (Rewinded > 1),
            Rewinded = Rewinded - 1;
            if CurrentIndex < ActualBufferSize,
                CurrentIndex = CurrentIndex + 1;
            else
                CurrentIndex = 1;
            end
        elseif  evnt.Character == 31,
            ClosingEvent = 1;
        elseif  evnt.Character == 30,
            StayOnCurrent = 0;
        elseif IsNum & ismember(NumPressed,[1:NEvents]),
            % Select Image Area
            ToBeUpdated.Menu = 1;
            Boxing = 1;
            EventNumber = NumPressed;
            PlotImage_Captioning(FramesBuffer{CurrentIndex},FigureHandle,FrameNumbersBuffer(CurrentIndex));
%             Boxing = 0;
%             rect = [];
%             while isempty(rect),
%                 rect = getrect(HImage);
%             end
%             rect = fix(rect);
            % Open that event
            StartFrame = FrameNumbersBuffer(CurrentIndex);
            CF = StartFrame;
            UniqueID = UniqueID + 1;
            while CF <= CurrentK & NCurrEvents(CF+1) < MaxCurrentEvents,
                NCurrEvents(CF+1) = NCurrEvents(CF+1) + 1;
                CurrentEvents(CF+1,NCurrEvents(CF+1)).StartFrame = StartFrame;
                CurrentEvents(CF+1,NCurrEvents(CF+1)).EventCode = NumPressed;
                CurrentEvents(CF+1,NCurrEvents(CF+1)).UniqueID = UniqueID;
                CF = CF + 1;
            end
%             fprintf(BackupCaptionFID,'%010d S %010d %03d %d %d %d %d\n',UniqueID, StartFrame, NumPressed,rect(1),rect(2),rect(3),rect(4));
            fprintf(BackupCaptionFID,'%010d S %010d %03d\n',UniqueID, StartFrame, NumPressed);
        end
        ToBeUpdated.Image = 4;
        ToBeUpdated.Info = 4;
        ToBeUpdated.Menu = 4;
        ToBeUpdated.Events = 4;
    end
elseif InPause,
    ToBeUpdated.Image = 0;
    if evnt.Character == 30,
        InPause = 0;
        ToBeUpdated.Image = 1;
        ToBeUpdated.Info = 1;
        ToBeUpdated.Menu = 1;
        ToBeUpdated.Events = 1;
    end
else
    if evnt.Character == 27,
        ESCDetected = 1;
    elseif evnt.Character == 29,
        if Cnt < MaxCnt,
            Cnt = Cnt * 2;
            FramesLeft = Cnt;
            ToBeUpdated.Info = 1;
        else
            beep;
        end
    elseif  evnt.Character == 28,
        if Cnt > 1,
            Cnt = Cnt / 2;
            FramesLeft = Cnt;
            ToBeUpdated.Info = 1;
        else
            beep;
        end
    elseif  evnt.Character == 31,
        if InBufferIndex == 1,
            CurrentIndex = ActualBufferSize;
        else
            CurrentIndex = InBufferIndex - 1;
        end
        Rewinded = 1;
        Cnt = 1;
        FramesLeft = 1;
        StayOnCurrent = 1;
        ToBeUpdated.Info = 1;
        ToBeUpdated.Menu = 1;
        ToBeUpdated.Events = 1;
    elseif  evnt.Character == 30,
        ToBeUpdated.Menu = 1;
        InPause = 1;
    end
end

