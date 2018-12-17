function PlotImage_Review(Image, FigureHandle, FrameNumber);
global Cnt;
global MaxCnt;
global DimX;
global DimY;
global DimZ;
global BlankEventsPanel;
global BlankMenuPanel;
global BlankInfoPanel;
global PanelDimX;
global Spacer;
global ToBeUpdated;
global Events;
global NEvents;
global PanelSize;
global StayOnCurrent;
global InPause;
global ClosingEvent;
global MaxCurrentEvents;
global CurrentEvents;
global NCurrEvents;
global HImage;
global HMenu;
global HInfo;
global HEvents;
global EventNumber;
global Boxing;
global slashstr;


set(FigureHandle,'DoubleBuffer','on');
set(FigureHandle,'Resize','off');
set(FigureHandle,'Color','k');
set(FigureHandle,'ToolBar','none');
set(FigureHandle,'MenuBar','none');
set(FigureHandle,'Name','Video Captioning');
set(FigureHandle,'NumberTitle','off');
rect = get(FigureHandle,'Position');
rect(4) = DimY + (4 * DimY * Spacer);
rect(3) = DimX + PanelDimX + (6 * DimY * Spacer);
set(FigureHandle,'Position',rect);



if ToBeUpdated.Info,
    % Plot The Info
    figure(FigureHandle);
    axes(HInfo);
    cla(HInfo);
    image(BlankInfoPanel);
    set(HInfo,'XTick',[],'YTick',[]);
    text(5,15,['Frame : ' num2str(FrameNumber,'%06d')],'FontSize',8);
    text(90,15,['Speed : ' num2str(Cnt,'%02d') slashstr num2str(MaxCnt,'%02d')],'FontSize',8);
    ToBeUpdated.Info = ToBeUpdated.Info - 1;
end

if ToBeUpdated.Menu,
    % Plot The Info
    figure(FigureHandle);
    axes(HMenu);
    cla(HMenu);
    image(BlankMenuPanel);
    set(HMenu,'XTick',[],'YTick',[]);
    if InPause,
        h = text(round(PanelDimX/2),15,['MAIN MENU']);
        set(h,'HorizontalAlignment','center','FontSize',8);
        h = text(round(PanelDimX/2),40,['RESUME']);
        set(h,'HorizontalAlignment','center','FontSize',8);
        Arr.X = [round(PanelDimX/2)-5,round(PanelDimX/2),round(PanelDimX/2),round(PanelDimX/2)+5];
        Arr.Y = [56,56-2,56+2,56];
        Arr.VX = [-20,0, 0, 20];
        Arr.VY = [0,-12, 12, 0];
        hold on;
        h = quiver(Arr.X,Arr.Y,Arr.VX,Arr.VY, 0);
        set(h,'Color','r','LineWidth',2,'AutoScale','off','ShowArrowHead','off','Visible','on');
        set(h,'MaxHeadSize',0.2)
        set(h,'ShowArrowHead','on')
    else
        h = text(round(PanelDimX/2),15,['MAIN MENU']);
        set(h,'HorizontalAlignment','center','FontSize',8);
        h = text(round(PanelDimX/2),40,['PAUSE']);
        set(h,'HorizontalAlignment','center','FontSize',8);
        h = text(round(PanelDimX/5),60,['SLOW']);
        set(h,'HorizontalAlignment','center','FontSize',8);
        h = text(round(PanelDimX*4/5),60,['FAST']);
        set(h,'HorizontalAlignment','center','FontSize',8);
        Arr.X = [round(PanelDimX/2)-5,round(PanelDimX/2),round(PanelDimX/2),round(PanelDimX/2)+5];
        Arr.Y = [56,56-2,56+2,56];
        Arr.VX = [-20,0, 0, 20];
        Arr.VY = [0,-12, 12, 0];
        hold on;
        h = quiver(Arr.X,Arr.Y,Arr.VX,Arr.VY, 0);
        set(h,'Color','r','LineWidth',2,'AutoScale','off','ShowArrowHead','off','Visible','on');
        set(h,'MaxHeadSize',0.2)
        set(h,'ShowArrowHead','on')
        h = text(round(PanelDimX/2),120,['ESC - Exit']);
        set(h,'HorizontalAlignment','center','FontSize',8);
    end
    ToBeUpdated.Menu = ToBeUpdated.Menu - 1;
end

if ToBeUpdated.Events,
    figure(FigureHandle);
    % Plot The Info
    axes(HEvents);
    cla(HEvents);
    image(BlankEventsPanel);
    set(HEvents, 'XTick',[],'YTick',[]);
    h = text(round(PanelDimX/2),15,['CURRENT EVENTS']);
    set(h,'HorizontalAlignment','center','FontSize',8);
    for i = 1:NCurrEvents(FrameNumber+1),
        Description = Events(CurrentEvents(FrameNumber+1,i).EventCode).Description;
        text(15,30+15*i,[num2str(CurrentEvents(FrameNumber+1, i).StartFrame,'%06d'),' - ', Description],'FontSize',8);
    end
    ToBeUpdated.Events = ToBeUpdated.Events - 1;
end


if ToBeUpdated.Image,
    figure(FigureHandle);
    % Plot The Image
    axes(HImage);
    cla(HImage);
    image(Image);
    set(HImage,'XTick',[],'YTick',[]);
    drawnow;
    ToBeUpdated.Image = ToBeUpdated.Image - 1;
end
    
return