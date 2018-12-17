function [] = SaveFromBackupToCaptionFile
global Events;
global NEvents;
global NCurrEvents;
global CurrentEvents;
global MaxCurrentEvents;
global CaptionFileName;
global BackupCaptionFileName;

% Determine How many events are there in 
BackupCaptionFID = fopen(BackupCaptionFileName,'r');
% [Line] = fscanf(BackupCaptionFID,'%010d %c %010d %03d %d %d %d %d',8);
[Line] = fscanf(BackupCaptionFID,'%010d %c %010d %03d',8);
NumOfEvents = 0;
while ~feof(BackupCaptionFID),
    ThisUniqueID = Line(1);
    if ThisUniqueID > NumOfEvents,
        NumOfEvents = ThisUniqueID;
    end
%     [Line] = fscanf(BackupCaptionFID,'%010d %c %010d %03d %d %d %d %d',8);
    [Line] = fscanf(BackupCaptionFID,'%010d %c %010d %03d',8);
end
fclose(BackupCaptionFID);

clear CurrentEvents;
clear NCurrEvents;
EventList = struct('UniqueID',cell(NumOfEvents,1),...
                   'EventCode',cell(NumOfEvents,1),...
                   'StartFrame',cell(NumOfEvents,1),...
                   'EndFrame',cell(NumOfEvents,1),...
                   'StartRect',cell(NumOfEvents,1),...
                   'EndRect',cell(NumOfEvents,1),...
                   'Description',cell(NumOfEvents,1));



% Read the entries of the Backup file
BackupCaptionFID = fopen(BackupCaptionFileName,'r');
CaptionFID = fopen([CaptionFileName(1:end-3) 'txt'],'w');
OrderingIndex = zeros(NumOfEvents,1);          
% [Line] = fscanf(BackupCaptionFID,'%010d %c %010d %03d %d %d %d %d',8);
% fprintf(CaptionFID,'%010d %c %010d %03d %d %d %d %d\n',[Line]);
[Line] = fscanf(BackupCaptionFID,'%010d %c %010d %03d',8);
fprintf(CaptionFID,'%010d %c %010d %03d\n',[Line]);
while ~feof(BackupCaptionFID),
    EventList(Line(1)).UniqueID = Line(1);
    EventList(Line(1)).EventCode = Line(4);
    if Line(2) == 'S',
        EventList(Line(1)).StartFrame = Line(3);
%         EventList(Line(1)).StartRect = Line(4:8);
        OrderingIndex(Line(1)) = Line(3);
        if isempty(EventList(Line(1)).EndFrame), 
            EventList(Line(1)).EndFrame = Inf;
        end
    else
        EventList(Line(1)).EndFrame = Line(3);
%         EventList(Line(1)).EndRect = Line(4:8);
    end
    EventList(Line(1)).Description = Events(Line(4)).Description;
%     [Line] = fscanf(BackupCaptionFID,'%010d %c %010d %03d %d %d %d %d',8);
%     fprintf(CaptionFID,'%010d %c %010d %03d %d %d %d %d\n',[Line]);
    [Line] = fscanf(BackupCaptionFID,'%010d %c %010d %03d',8);
    fprintf(CaptionFID,'%010d %c %010d %03d\n',[Line]);
end
fclose(BackupCaptionFID);
fclose(CaptionFID);

% Reorder the Even List
[delme, Order] = sort(OrderingIndex, 'ascend');
EventList = EventList(Order);
save(CaptionFileName,'Events','NEvents','EventList','NumOfEvents');

