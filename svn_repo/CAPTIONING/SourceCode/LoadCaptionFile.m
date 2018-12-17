function [] = LoadCaptionFile();
global Events;
global NEvents;
global NCurrEvents;
global CurrentEvents;
global MaxCurrentEvents;
global UniqueID;
global intStartFrm;          
global intSkipFrms;
global intNFrms;
global CaptionFileName;
global BackupCaptionFileName;
global BufferSize;


MaxCurrentEvents = 9;
CurrentEvents = ...
    struct('StartFrame',cell(intNFrms - intStartFrm,MaxCurrentEvents),...
           'EventCode',cell(intNFrms - intStartFrm,MaxCurrentEvents),...
           'UniqueID',cell(intNFrms - intStartFrm,MaxCurrentEvents));
NCurrEvents = zeros(intNFrms - intStartFrm,1);

MaxUniqueID = 0;    
MaxFrameNum = 0;
load(CaptionFileName);
ListLength = size(EventList,1);
intStartFrm = 0; 
intSkipFrms = 0;



NActiveEvents = 0;
ActiveEvents = [];
NextActive = 1;
for K = 0:intNFrms-1,
    % Active is what we had before
    NewActiveEvents = [];

    % Minus whatever Ends now
    j = 1;
    for i = 1:NActiveEvents,
        if EventList(ActiveEvents(i)).EndFrame > K,
            NewActiveEvents = [NewActiveEvents ActiveEvents(i)];
        end
    end
    % Plus whatever Starts now
    if NextActive <= NumOfEvents;
        while EventList(NextActive).StartFrame == K,
            NewActiveEvents = [NewActiveEvents NextActive];
            NextActive = NextActive + 1;
            if NextActive > NumOfEvents
                break;
            end
        end
    end
    ActiveEvents = NewActiveEvents;
    NActiveEvents = length(ActiveEvents);
    
    NCurrEvents(K+1) = NActiveEvents;
    for i = 1:NActiveEvents,
        CurrentEvents(K+1,i).StartFrame = EventList(ActiveEvents(i)).StartFrame;
        CurrentEvents(K+1,i).EventCode = EventList(ActiveEvents(i)).EventCode;
        CurrentEvents(K+1,i).UniqueID = EventList(ActiveEvents(i)).UniqueID;
    end
end
