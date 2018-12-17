function PlotImage_Captioning(Image,FigureHandle,FrameNumber)
global strVideoFExt;
global strInVideoFName;
global strInVideoPath;
global findex;
global findex1;
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
global firstfile;
global QuitPrg;
global filechg;
global slashstr;
global width1;
global height1;

set(FigureHandle,'DoubleBuffer','on');
set(FigureHandle,'Resize','off');
set(FigureHandle,'Color','k');
set(FigureHandle,'ToolBar','none');
set(FigureHandle,'MenuBar','none');
set(FigureHandle,'Name',['Drosophila Behavior Captioning - ' strInVideoPath slashstr strInVideoFName '.' strVideoFExt]);
set(FigureHandle,'NumberTitle','off');
rect = get(FigureHandle,'Position');
rect(4) = DimY + (4 * DimY * Spacer);
rect(3) = DimX + PanelDimX + (6 * DimY * Spacer);
set(FigureHandle,'Position',rect);

if ~findex1,
    f = uimenu('Label','File');
    uimenu(f,'Label','Open','enable','on','Callback','OpenAVI');
    uimenu(f,'Label','Quit [Esc]','enable','on','Callback','exitprg');
    if ~firstfile,
        OpenAVI;
        firstfile = 1;
    end
    QuitPrg = 0;
    filechg = 0;
    findex1 = 1;
    if ~findex,
        Quitprg = 1;
    else
        Quitprg = 0;
    end
end

if ToBeUpdated.Info,
    % Plot The Info
    figure(FigureHandle);
    axes(HInfo);
    cla(HInfo);
    image(BlankInfoPanel);
    set(HInfo,'XTick',[],'YTick',[]);
    text(5,15,['Frame : ' num2str(FrameNumber,'%06d')],'FontSize',8,'color','w');
    text(90,15,['Speed : ' num2str(Cnt,'%02d') '/' num2str(MaxCnt,'%02d')],'FontSize',8,'color','w');
    ToBeUpdated.Info = ToBeUpdated.Info - 1;
end

if ToBeUpdated.Menu,
    % Plot The Info
    figure(FigureHandle);
    axes(HMenu);
    cla(HMenu);
    image(BlankMenuPanel);
    set(HMenu,'XTick',[],'YTick',[]);
    if StayOnCurrent,
%         if Boxing,
%             Description = Events(EventNumber).Description;
%             h = text(round(PanelDimX/2),105,['DRAW BOX AROUND']);
%             set(h,'HorizontalAlignment','center','FontSize',8,'color','w');
%             h = text(round(PanelDimX/2),120,['AROUND']);
%             set(h,'HorizontalAlignment','center','FontSize',8,'color','w');
%             h = text(round(PanelDimX/2),150,[Description]);
%             set(h,'HorizontalAlignment','center','FontSize',10,'color','w');
%         else
            if ClosingEvent,
                h = text(round(PanelDimX/2),15,['CLOSE EVENTS']);
                set(h,'HorizontalAlignment','center','FontSize',8,'color','w');
                h = text(round(PanelDimX/2),40,['RESUME']);
                set(h,'HorizontalAlignment','center','FontSize',8,'color','w');
                Arr.X = [round(PanelDimX/2)-5,round(PanelDimX/2),round(PanelDimX/2),round(PanelDimX/2)+5];
                Arr.Y = [56,56-2,56+2,56];
                Arr.VX = [-20,0, 0, 20];
                Arr.VY = [0,-12, 12, 0];
                hold on;
                h = quiver(Arr.X,Arr.Y,Arr.VX,Arr.VY, 0);
                set(h,'Color','r','LineWidth',2,'AutoScale','off','ShowArrowHead','off','Visible','on');
                set(h,'MaxHeadSize',0.2)
                set(h,'ShowArrowHead','on')
                for i = 1:NCurrEvents(FrameNumber+1),
                    Description = Events(CurrentEvents(FrameNumber+1,i).EventCode).Description;
                    text(15,90+15*i,[num2str(i),' - Close ', Description],'color','w');
                end
            else
                h = text(round(PanelDimX/2),15,['OPEN/CLOSE EVENTS']);
                set(h,'HorizontalAlignment','center','FontSize',8, 'Color', 'w');
                h = text(round(PanelDimX/2),40,['RESUME']);
                set(h,'HorizontalAlignment','center','FontSize',8, 'Color', 'w');
                h = text(round(PanelDimX/4),60,['REW']);
                set(h,'HorizontalAlignment','center','FontSize',8, 'Color', 'w');
                h = text(round(PanelDimX*3/4),60,['FF']);
                set(h,'HorizontalAlignment','center','FontSize',8, 'Color', 'w');
                h = text(round(PanelDimX/2),80,['CLOSE']);
                set(h,'HorizontalAlignment','center','FontSize',8, 'Color', 'w');
                Arr.X = [round(PanelDimX/2)-5,round(PanelDimX/2),round(PanelDimX/2),round(PanelDimX/2)+5];
                Arr.Y = [56,56-2,56+2,56];
                Arr.VX = [-20,0, 0, 20];
                Arr.VY = [0,-12, 12, 0];
                hold on;
                h = quiver(Arr.X,Arr.Y,Arr.VX,Arr.VY, 0);
                set(h,'Color','r','LineWidth',2,'AutoScale','off','ShowArrowHead','off','Visible','on');
                set(h,'MaxHeadSize',0.2)
                set(h,'ShowArrowHead','on')
                for i = 1:NEvents,
                    text(35,90+15*i,[num2str(i),' - ', Events(i).Description],'color','w');
                end
            end
%         end
    elseif InPause,
        h = text(round(PanelDimX/2),15,['MAIN MENU']);
        set(h,'HorizontalAlignment','center','FontSize',8,'color','w');
        h = text(round(PanelDimX/2),40,['RESUME']);
        set(h,'HorizontalAlignment','center','FontSize',8,'color','w');
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
        set(h,'HorizontalAlignment','center','FontSize',8,'Color','w');
        h = text(round(PanelDimX/2),40,['PAUSE']);
        set(h,'HorizontalAlignment','center','FontSize',8,'Color','w');
        h = text(round(PanelDimX/5),60,['SLOW']);
        set(h,'HorizontalAlignment','center','FontSize',8,'Color','w');
        h = text(round(PanelDimX*4/5),60,['FAST']);
        set(h,'HorizontalAlignment','center','FontSize',8,'Color','w');
        h = text(round(PanelDimX/2),80,['EVENT']);
        set(h,'HorizontalAlignment','center','FontSize',8,'Color','w');
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
        set(h,'HorizontalAlignment','center','FontSize',8,'color','w');
    end
    ToBeUpdated.Menu = ToBeUpdated.Menu - 1;
end

if ToBeUpdated.Events,
    figure(FigureHandle);
    % Plot The Info
    axes(HEvents);
    cla(HEvents);
    image(BlankEventsPanel);
    set(HEvents, 'XTick',[],'YTick',[],'Color','w');
    h = text(round(PanelDimX/2),15,['CURRENT EVENTS']);
    set(h,'HorizontalAlignment','center','FontSize',8,'color','w');
    for i = 1:NCurrEvents(FrameNumber+1),
        Description = Events(CurrentEvents(FrameNumber+1,i).EventCode).Description;
        text(15,30+15*i,[num2str(CurrentEvents(FrameNumber+1, i).StartFrame,'%06d'),' - ', Description],'FontSize',8,'color','w');
    end
    ToBeUpdated.Events = ToBeUpdated.Events - 1;
end


if ToBeUpdated.Image,
    figure(FigureHandle);
    % Plot The Image
    axes(HImage);
    cla(HImage);
    image(flipdim(flipdim(permute(reshape(Image, 3, width1, height1),[3 2 1]),1),3));
    set(HImage,'XTick',[],'YTick',[]);
    drawnow;
    ToBeUpdated.Image = ToBeUpdated.Image - 1;
end
return


function OpenAVI
global strVideoFExt;
global strInVideoFName;
global strInVideoPath;
global findex;
global filechg;
global InPause;

[strInVideoFName, strInVideoPath, findex] = uigetfile({'*.avi;*.wmv','AVI/WMV file (*.avi,*.wmv)';'*.fmf','fmf-file (*.fmf)'}, ...
                                                        'Open movie file','MultiSelect','off');
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
