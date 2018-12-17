function KeyboardHandling_Review(src,evnt)
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
global HImage;
global HMenu;
global HInfo;
global HEvents;
global EventNumber;
global Boxing;
global FigureHandle;

if InPause,
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
    elseif  evnt.Character == 30,
        ToBeUpdated.Menu = 1;
        InPause = 1;
    end
end

