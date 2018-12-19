function gVision
% gVision - Video Capture Tools for Ethology 
%   and other Machine Vision Applications
%
% Version 1.1b - Updated for Matlab 2010b - 9/14/2010
%
% gVision is distributed under the GNU Public License
% http://www.gnu.org/licenses/gpl.txt
%
% Gus K. Lott, PhD (c)2010
% lottg@janelia.hhmi.org
% HHMI Janelia Farm Research Campus
% Neurobiological Instrumentation Engineer
%
% TODO:
%
% - Thoroughly comment code
% - Custom DCAM support
% - State Saving currently only activates the gvBasic mode; save method?
% - Populate About Menu
%
% - Synchronize Multiple gVisions through a file (or otherwise) socket



warning off imaq:peekdata:tooManyFramesRequested
warning off MATLAB:aviclose:noAssociatedFrames
warning off MATLAB:timer:deleterunning

makegui('1.1b - Gus K Lott III, PhD - HHMI Janelia Farm 2010')



%Capture a snapshot of the video frame as a tiff/jpeg/whatever
function snapShot(~,~,fig) 
gui = get(fig,'userdata');

im = getsnapshot(gui.vi);
tfig = figure; imshow(im)
drawnow

[fname, pname] = uiputfile({'*.jpg','JPEG Image (*.jpg)';'*.tif','TIFF Image (*.tif)';'*.png','PNG Image (*.png)'},'Save Snap Shot',fullfile(get(gui.menu(1),'userdata'),'snapshot.jpg'));
if isequal(fname,0) || isequal(pname,0)
    delete(tfig)
    return
end
set(gui.menu(1),'userdata',pname)
imwrite(im,fullfile(pname,fname),fname(end-2:end))
delete(tfig)


%Connect to a camera and configure the UI - 
% Called when a camera mode is selected from the the Camera menu
function addCam(obj,event,fig,flag)
gui=get(fig,'userdata');
try
    
imaqreset
if flag == 0 %pulls camera connect info from menu items
    mode = get(obj,'label');
    id = get(get(obj,'parent'),'userdata');
    adaptor = get(get(get(obj,'parent'),'parent'),'label');
else %Called when loading a saved state
    mode = obj.format;
    id = obj.id;
    adaptor = obj.adaptor;
end
gui.adaptor = adaptor; %i.e. dcam, winvideo, etc

%Create video object
gui.vi = videoinput(adaptor,id,mode);
gui.vi.framespertrigger = inf;
gui.vi.triggerrepeat = 0;

%Update UI
set(gui.FPT,'string','inf','userdata',inf)
set(gui.TR,'string',0,'userdata',0)

%Create the preview figure
delete(findobj('tag','gVIprev'))
ROI = gui.vi.ROIPosition;
gui.vidFig = figure('units','pixels','position',[0 0 ROI(3) ROI(4)],'numbertitle','off','menubar','none',...
    'name','gVision Video Preview - Gus K Lott III, PhD - HHMI Janelia Farm','tag','gVIprev');
centerfig(gui.vidFig)
set(gui.vidFig,'closerequestfcn',{@prevFigClose,gui.fig});

gui.prevAx=axes('parent',gui.vidFig,'position',[0 0 1 1],'xtick',[],'ytick',[],'box','on','color','w');
vidRes = get(gui.vi, 'VideoResolution');
nBands = get(gui.vi, 'NumberOfBands');
gui.pImage = image( zeros(vidRes(2), vidRes(1), nBands) );
set(gui.prevAx,'xtick',[],'ytick',[],'xlim',[0 vidRes(1)+1],'ylim',[0 vidRes(2)+1])
preview(gui.vi,gui.pImage) %start preview into it

%Enable UI 
set(gui.gROI,'string',num2str(get(gui.vi,'ROIPosition')),'enable','on')
set([gui.fileSelect gui.addSROI gui.remSROI gui.dispSROI gui.StartStop gui.snapShot gui.sROInameRoot gui.sROInameNum gui.sROIcolor],'enable','on')
set(gui.previewONOFF,'enable','on','value',1)
set([gui.FPT gui.TR],'enable','on')
set([gui.FPStxt,gui.FPS],'visible','off')

%Configure Alignment Crosshairs in the preview window
gui.pxH = line([-100000 100000],[nan nan],'color','r','linewidth',2,'linestyle','--','buttondownfcn',{@dragCross,gui.fig});
gui.pxHt = text(0,0,'','parent',gui.prevAx,'color','r','fontweight','bold');
gui.pxV = line([nan nan],[-100000 100000],'color','r','linewidth',2,'linestyle','--','buttondownfcn',{@dragCross,gui.fig});
gui.pxVt = text(0,0,'','parent',gui.prevAx,'color','r','fontweight','bold');
set([gui.xHair gui.histActivate gui.MetaData gui.gROI gui.gROIrst gui.gROIdraw],'enable','on')

%draw histogram bars
delete(get(gui.histAx,'children'));
if gui.vi.NumberOfBands == 1
    for j = 0:255
        gui.histBars(j+1)=patch([-.5 .5 .5 -.5]+j, [0 0 0 0],'b','edgecolor','none','parent',gui.histAx);
    end
    set(gui.histBars(1),'xdata',[-5 .5 .5 -5],'facecolor','r')
    set(gui.histBars(256),'xdata',[-.5 5 5 -.5]+255,'facecolor','r')
    set(gui.histAx,'xtick',[],'ytick',[],'xlim',[-5 260],'ylim',[0 10^3],'box','on')
else %if RGB image
    c = 'rgb';
    for i = 1:3
        for j = 0:255
            gui.histBars(i,j+1)=patch([-.5 .5 .5 -.5]+j, [0 0 0 0],c(i),'edgecolor','none','parent',gui.histAx);
        end
        switch i
            case 1; modc = [1 .8 .8];
            case 2; modc = [.8 1 .8];
            case 3; modc = [.8 .8 1];
        end
        set(gui.histBars(1),'xdata',[-5 .5 .5 -5],'facecolor',modc)
        set(gui.histBars(256),'xdata',[-.5 5 5 -.5]+255,'facecolor',modc)
        set(gui.histAx,'nextplot','add')
    end
    set(gui.histAx,'xtick',[],'ytick',[],'xlim',[-5 260],'ylim',[0 10^3],'box','on') 
end

%hist box for drawing a limited region to calculate hist from
ROI = gui.vi.ROIPosition;
hist.ROI = [1 1 ROI(3)-1 ROI(4)-1];
ROI = hist.ROI;
    hist.lines(1) = line([0 ROI(3)]+ROI(1),[0 0]+ROI(2),'parent',gui.prevAx,'buttondownfcn',{@dragBox,gui.fig,1});
    hist.lines(2) = line([ROI(3) ROI(3)]+ROI(1),[0 ROI(4)]+ROI(2),'parent',gui.prevAx,'buttondownfcn',{@dragBox,gui.fig,2});
    hist.lines(3) = line([ROI(3) 0]+ROI(1),[ROI(4) ROI(4)]+ROI(2),'parent',gui.prevAx,'buttondownfcn',{@dragBox,gui.fig,3});
    hist.lines(4) = line([0 0]+ROI(1),[ROI(4) 0]+ROI(2),'parent',gui.prevAx,'buttondownfcn',{@dragBox,gui.fig,4});
set(hist.lines,'linewidth',3,'color',[0 1 1],'userdata',hist,'visible','off','marker','.','tag','gHist')
gui.histLines = hist.lines;


%Line scan plot line
gui.lineScanTrace = line(nan,nan,'color','b','linewidth',1,'parent',gui.histAx,'visible','off');


% Detect 8 vs 16 bit images for tif only logging
info = imaqhwinfo(gui.vi);
if ~strcmpi(info.NativeDataType,'uint8')
    file = get(gui.fName,'userdata');
    file.fname = 'temp.tif';
    set(gui.fName,'string',fullfile(file.pname,file.fname),'userdata',file)
end

%Enable colormap changes for grayscale streams
if get(gui.vi,'NumberOfBands') == 1
    set(gui.streamMen(1),'enable','on')
    set(gui.cmapMen,'checked','off')
    set(gui.cmapMen(1),'checked','on')
else
    set(gui.streamMen(1),'enable','off')
end
set(gui.fig,'userdata',gui)

updatePropStat(fig); %Read properties from camera
parsePlugs(fig); %Install plugins
ModeChange([],[],fig); %Enable basic mode
catch e
    save(sprintf('Exception_%s.mat',datestr(now)),'e')
    errordlg(e.message)
end

%Create and apply a colormap to a grayscale preview window
function changePreviewCMap(obj,event,fig)
gui = get(fig,'userdata');
switch get(obj,'Label')
    case 'Default'
        cmap = 'default';
    case 'Hot'
        cmap = hot(256);
    case 'Cool'
        cmap = cool(256);
    case 'HSV'
        cmap = hsv(256);
    case 'Jet'
        cmap = jet(256);
end
colormap(gui.prevAx,cmap)

%Validate Recording Property string entries
function setRec(obj,event,fig)
gui = get(fig,'userdata');
val = str2double(get(obj,'string'));
if isnan(val)
    set(obj,'string',num2str(get(obj,'userdata')))
    return
end
switch(obj)
    case gui.FPT; val = round(val);
        if val<=0
            set(obj,'string',num2str(get(obj,'userdata')))
            return
        end
    case gui.TR; val = round(val);
    case gui.FPS
        if val<=0
            set(obj,'string',num2str(get(obj,'userdata')))
            return
        end
end
if val<0
    set(obj,'string',num2str(get(obj,'userdata')))
    return
end
set(obj,'userdata',val,'string',num2str(val))

switch(obj)
    case gui.FPT; gui.vi.FramesPerTrigger = val;
    case gui.TR; gui.vi.TriggerRepeat = val;
end

%Toggle visualization of the video stream
function switchPreview(obj,event,fig)
gui = get(fig,'userdata');
switch get(obj,'value')
    case 0
        stoppreview(gui.vi)
        set(gui.vidFig,'visible','off')
    case 1
        preview(gui.vi,gui.pImage)
        set(gui.vidFig,'visible','on')
end

%Check for trigger selection change
function ModeChange(~,event,fig)
gui = get(fig,'userdata');
data = get(gui.ModeOptions,'data');

if isempty(event)
    for i=1:size(data,1)
        if data{i,1}
            break;
        end
    end
    ind = i;
else
    ind = event.Indices(1);
end

%Delete an existing class
if ~isempty(gui.mode); delete(gui.mode); gui.mode = [];end

% updatePropStat(gui.fig)
for i = 1:size(data,1)
    data{i,1} = false;
end
data{ind,1}=true;
set(gui.ModeOptions,'data',data)

classnames = get(gui.ModeOptions,'userdata');

%Create a class of the name specified
eval(sprintf('gui.mode = plugins.%s(fig);',classnames{ind}))
set(fig,'userdata',gui);



%handling global ROI
function setGROI(obj,event,fig)
gui = get(fig,'userdata');
ROI = str2num(get(obj,'string'));
if length(ROI)~=4; 
    set(obj,'string',num2str(gui.vi.ROIPosition))
    return; 
end

res = gui.vi.VideoResolution;

%check for illogical values
errf = 0;
if ROI(1)+ROI(3) > res(1); errf = 1; end
if ROI(2)+ROI(4) > res(2); errf = 1; end
if ROI(1)<0 | ROI(2)<0 | ROI(3)<1 | ROI(4)<1; errf = 1; end
if errf
    set(obj,'string',num2str(gui.vi.ROIPosition))
    return
end
gui.vi.ROIPosition = ROI;
ROI = gui.vi.ROIPosition;
set(obj,'string',num2str(ROI))
set(gui.prevAx,'xlim',[0 res(1)+1]-ROI(1),'ylim',[0 res(2)+1]-ROI(2))

%Recenter crosshairs
set(gui.pxH,'ydata',[1 1]*round(ROI(4)/2));
set(gui.pxHt,'string',num2str(round(ROI(4)/2)))
set(gui.pxHt,'position',[1 ROI(4)/2+13])
set(gui.pxV,'xdata',[1 1]*round(ROI(3)/2));
set(gui.pxVt,'string',num2str(round(ROI(3)/2)))
set(gui.pxVt,'position',[ROI(3)/2+13 13])
if get(gui.xHair,'value')
    set([gui.pxH gui.pxHt gui.pxV gui.pxVt],'visible','on')
else
    set([gui.pxH gui.pxHt gui.pxV gui.pxVt],'visible','off')
end

%make sure that histogram window still fits
ROI = gui.vi.ROIPosition;
ROI = [1 1 ROI(3)-1 ROI(4)-1];
set(gui.histLines(1),'xdata',[0 ROI(3)]+ROI(1),'ydata',[0 0]+ROI(2))
set(gui.histLines(2),'xdata',[ROI(3) ROI(3)]+ROI(1),'ydata',[0 ROI(4)]+ROI(2))
set(gui.histLines(3),'xdata',[ROI(3) 0]+ROI(1),'ydata',[ROI(4) ROI(4)]+ROI(2))
set(gui.histLines(4),'xdata',[0 0]+ROI(1),'ydata',[ROI(4) 0]+ROI(2))
if get(gui.histActivate,'value')
    set(gui.histLines,'visible','on')
else
    set(gui.histLines,'visible','off')
end

hist = get(gui.histLines(1),'userdata');
hist.ROI = ROI;
set(gui.histLines,'userdata',hist)

function drawGROI(obj,event,fig)
gui = get(fig,'userdata');
res = gui.vi.VideoResolution;
ROI = gui.vi.ROIPosition;

switch get(obj,'value')
    case 0 %apply the drawn boxes value
        set([gui.gROI gui.gROIrst],'enable','on')
        delete(findobj('tag','gROIdraw'))
        setGROI(gui.gROI,'',fig)
    case 1 %draw a box and let the user manually specify the ROI
        set([gui.gROI gui.gROIrst],'enable','off')
        box.lines(1) = line([1 ROI(3)],[1 1],'parent',gui.prevAx,'buttondownfcn',{@dragGROI,gui.fig,1},'tag','gROIdraw','color','r','linewidth',3);
        box.lines(2) = line([ROI(3) ROI(3)],[1 ROI(4)],'parent',gui.prevAx,'buttondownfcn',{@dragGROI,gui.fig,2},'tag','gROIdraw','color','r','linewidth',3);
        box.lines(3) = line([ROI(3) 1],[ROI(4) ROI(4)],'parent',gui.prevAx,'buttondownfcn',{@dragGROI,gui.fig,3},'tag','gROIdraw','color','r','linewidth',3);
        box.lines(4) = line([1 1],[ROI(4) 1],'parent',gui.prevAx,'buttondownfcn',{@dragGROI,gui.fig,4},'tag','gROIdraw','color','r','linewidth',3);
        box.ROI = ROI;
        set(box.lines,'userdata',box)
end
function dragGROI(obj,event,fig,side)
gui = get(fig,'userdata');

set(gui.vidFig,'windowbuttonmotionfcn',{@moveGROI,fig,obj,side})
set(gui.vidFig,'windowbuttonupfcn','set(gcbf,''windowbuttonmotionfcn'','''')');
function moveGROI(obj,event,fig,hline,side)
gui = get(fig,'userdata');
pos = get(gca,'currentpoint');
gROI = gui.vi.ROIPosition;
box = get(hline,'userdata');
ROI = box.ROI;

switch side
    case 3 %Bottom
        box.ROI(4) = round(pos(1,2))-ROI(2);
        if box.ROI(4)<0; box.ROI(4) = 0; end
    case 2 %Right
        box.ROI(3) = round(pos(1,1))-ROI(1);
        if box.ROI(3)<0; box.ROI(3) = 0; end
    case {1,4} %Top & Left
        box.ROI(2) = round(pos(1,2));
        box.ROI(1) = round(pos(1,1));
end
ROI = box.ROI;
if ROI(1)+ROI(3) > gROI(3); return; end
if ROI(2)+ROI(4) > gROI(4); return; end
if pos(1,2)<1; return; end
if pos(1,1)<1; return; end


set(box.lines(1),'xdata',[1 ROI(3)]+ROI(1),'ydata',[1 1]+ROI(2))
set(box.lines(2),'xdata',[ROI(3) ROI(3)]+ROI(1),'ydata',[1 ROI(4)]+ROI(2))
set(box.lines(3),'xdata',[ROI(3) 1]+ROI(1),'ydata',[ROI(4) ROI(4)]+ROI(2))
set(box.lines(4),'xdata',[1 1]+ROI(1),'ydata',[ROI(4) 1]+ROI(2))
try set(box.txt,'position',[ROI(1) ,ROI(2)+ROI(4)/2]); end

set(gui.gROI,'string',num2str(ROI))

set(box.lines,'userdata',box)
        


function rstGROI(obj,event,fig)
gui = get(fig,'userdata');
res = gui.vi.VideoResolution;
gui.vi.ROIPosition = [0 0 res(1) res(2)];
set(gui.gROI,'string',num2str(gui.vi.ROIPosition))

%Recenter alignment crosshairs
ROI = gui.vi.ROIPosition;
    set(gui.pxH,'ydata',[1 1]*round(ROI(4)/2));
    set(gui.pxHt,'string',num2str(round(ROI(4)/2)))
    set(gui.pxHt,'position',[1 ROI(4)/2+13])
    set(gui.pxV,'xdata',[1 1]*round(ROI(3)/2));
    set(gui.pxVt,'string',num2str(round(ROI(3)/2)))
    set(gui.pxVt,'position',[ROI(3)/2+13 13])
set(gui.prevAx,'xlim',[0 ROI(3)+1],'ylim',[0 ROI(4)+1])
if get(gui.xHair,'value')
    set([gui.pxH gui.pxHt gui.pxV gui.pxVt],'visible','on')
else
    set([gui.pxH gui.pxHt gui.pxV gui.pxVt],'visible','off')
end

%make sure that histogram window still fits
ROI = gui.vi.ROIPosition;
ROI = [1 1 ROI(3)-1 ROI(4)-1];
set(gui.histLines(1),'xdata',[0 ROI(3)]+ROI(1),'ydata',[0 0]+ROI(2))
set(gui.histLines(2),'xdata',[ROI(3) ROI(3)]+ROI(1),'ydata',[0 ROI(4)]+ROI(2))
set(gui.histLines(3),'xdata',[ROI(3) 0]+ROI(1),'ydata',[ROI(4) ROI(4)]+ROI(2))
set(gui.histLines(4),'xdata',[0 0]+ROI(1),'ydata',[ROI(4) 0]+ROI(2))

if get(gui.histActivate,'value')
    set(gui.histLines,'visible','on')
else
    set(gui.histLines,'visible','off')
end

hist = get(gui.histLines(1),'userdata');
hist.ROI = ROI;
set(gui.histLines,'userdata',hist)


%For handling software-subrois
function addSROI(~,event,fig)
gui = get(fig,'userdata');
ROIlist = get(gui.listSROI,'userdata');
gROI = gui.vi.ROIPosition;

if strcmpi(event,'loading')
    %When reloading the program
    ROI = [10 10 10 10];
else
    s=inputdlg('Left, Top, Width, Height','New Sub-ROI',1,{'10 10 100 100'});
    if isempty(s); return; end
    ROI = str2num(s{1});
end


if isempty(ROI); return; end
if length(ROI)~=4 | min(ROI)<0; return; end

if ROI(1)+ROI(3) > gROI(3)
    return
end
if ROI(2)+ROI(4) > gROI(4)
    return
end

ROIlist(end+1).ROI = ROI;

set(gui.sROIControls,'visible','off')

%Draw UI for control of this specific sub-roi
%Allow color selection on a per channel basis?
nROIs = length(gui.sROIControls)+1;

%Find next available name
sROIname = 1;
for i = 1:length(gui.sROIControls)
    ui = get(gui.sROIControls(i),'userdata');
    thisname = get(ui.name,'string');
    id = sscanf(thisname,'sub%02d');
    if id>=sROIname
        sROIname=id+1;
    end
end

%Draw channel details for a subROI
gui.sROIControls(end+1) = uipanel(gui.ROIPanel,'units','normalized','position',[.54 .01 .45 .75]);
pnt = gui.sROIControls(end);
   
uicontrol(pnt,'style','text','string','Name:','units','normalized','horizontalalignment','left',...
    'position',[.01 .9 .28 .09],'value',1,'backgroundcolor',get(pnt,'backgroundcolor'),'fontweight','bold');
ui.name = uicontrol(pnt,'style','edit','string',sprintf('sub%02d',sROIname),'units','normalized',...
    'position',[.25 .9 .59 .09],'value',1,'backgroundcolor','w','horizontalalignment','left','callback',{@updateSROI,gui.fig});
ui.color = uicontrol(pnt,'style','pushbutton','string','','units','normalized','position',[.85 .9 .14 .09],...
    'callback',{@colorSROI,gui.fig},'backgroundcolor',[0 .7 0]);

uicontrol(pnt,'style','text','string','Loc:','units','normalized','horizontalalignment','left',...
    'position',[.01 .8 .28 .09],'value',1,'backgroundcolor',get(pnt,'backgroundcolor'),'fontweight','bold');
ui.x = uicontrol(pnt,'style','edit','string',num2str(ROI(1)),'units','normalized',...
    'position',[.2 .8 .19 .09],'backgroundcolor','w','tooltipstring','x','callback',{@updateSROI,gui.fig});
ui.y = uicontrol(pnt,'style','edit','string',num2str(ROI(2)),'units','normalized',...
    'position',[.4 .8 .19 .09],'backgroundcolor','w','tooltipstring','y','callback',{@updateSROI,gui.fig});
ui.w = uicontrol(pnt,'style','edit','string',num2str(ROI(3)),'units','normalized',...
    'position',[.6 .8 .19 .09],'backgroundcolor','w','tooltipstring','w','callback',{@updateSROI,gui.fig});
ui.h = uicontrol(pnt,'style','edit','string',num2str(ROI(4)),'units','normalized',...
    'position',[.8 .8 .19 .09],'backgroundcolor','w','tooltipstring','h','callback',{@updateSROI,gui.fig});

ui.capture = uicontrol(pnt,'style','checkbox','string','Capture','units','normalized',...
    'position',[.5 .7 .5 .09],'value',1,'backgroundcolor',get(pnt,'backgroundcolor'),'callback',{@updateSROI,gui.fig});
ui.dispname = uicontrol(pnt,'style','checkbox','string','Show Name','units','normalized',...
    'position',[0 .7 .5 .09],'value',1,'backgroundcolor',get(pnt,'backgroundcolor'),'callback',{@updateSROI,gui.fig});

uicontrol(pnt,'style','text','string','Notes:','units','normalized','horizontalalignment','left',...
    'position',[.01 .62 .28 .08],'value',1,'backgroundcolor',get(pnt,'backgroundcolor'),'fontweight','bold');
ui.notes = uicontrol(pnt,'style','edit','max',2,'backgroundcolor','w','units','normalized',...
    'position',[.01 .01 .98 .61],'horizontalalignment','left');


set(pnt,'userdata',ui)
set([gui.gROI gui.gROIdraw gui.gROIrst],'enable','off')
set(gui.sROIpopDown,'enable','on')
set(gui.fig,'userdata',gui)

%Update lines in overlay
ind = length(ROIlist);
drawROI(ROIlist,fig);

ROIlist = get(gui.listSROI,'userdata');
set(ROIlist(ind).lines,'linewidth',2)
set(gui.listSROI,'value',ind)
dispSROI(gui.dispSROI,'',fig)

%Called when a channel name is changed, update name in display window and
%string list
function updateSROI(~,~,fig)
gui = get(fig,'userdata');
ROIlist = get(gui.listSROI,'userdata');

for i = 1:length(gui.sROIControls)
    ui = get(gui.sROIControls(i),'userdata');
    a = str2double(get(ui.x,'string'));
    if ~isnan(a); ROIlist(i).ROI(1) = a;
    else set(ui.x,'string',num2str(ROIlist(i).ROI(1))); end
    a = str2double(get(ui.y,'string'));
    if ~isnan(a); ROIlist(i).ROI(2) = a; 
    else set(ui.y,'string',num2str(ROIlist(i).ROI(2))); end
    a = str2double(get(ui.w,'string'));
    if ~isnan(a); ROIlist(i).ROI(3) = a;
    else set(ui.w,'string',num2str(ROIlist(i).ROI(3))); end
    a = str2double(get(ui.h,'string'));
    if ~isnan(a); ROIlist(i).ROI(4) = a; 
    else set(ui.h,'string',num2str(ROIlist(i).ROI(4))); end
end

drawROI(ROIlist,fig)


%Main function called to re-draw the sub-ROI in the main window
function drawROI(ROIlist,fig)
gui = get(fig,'userdata');
delete(findobj('parent',gui.prevAx,'tag','gSROI'))

s={};
for i = 1:length(ROIlist)
    ui = get(gui.sROIControls(i),'userdata');
    s{i} = get(ui.name,'string');
    borderColor = get(ui.color,'backgroundcolor');
    
    ROI = ROIlist(i).ROI;
    ROIlist(i).lines(1) = line([0 ROI(3)]+ROI(1),[0 0]+ROI(2),'parent',gui.prevAx,'buttondownfcn',{@dragSROI,gui.fig,i,1});
    ROIlist(i).lines(2) = line([ROI(3) ROI(3)]+ROI(1),[0 ROI(4)]+ROI(2),'parent',gui.prevAx,'buttondownfcn',{@dragSROI,gui.fig,i,2});
    ROIlist(i).lines(3) = line([ROI(3) 0]+ROI(1),[ROI(4) ROI(4)]+ROI(2),'parent',gui.prevAx,'buttondownfcn',{@dragSROI,gui.fig,i,3});
    ROIlist(i).lines(4) = line([0 0]+ROI(1),[ROI(4) 0]+ROI(2),'parent',gui.prevAx,'buttondownfcn',{@dragSROI,gui.fig,i,4});
    ROIlist(i).txt = text(ROI(1),ROI(2),s{i},'parent',gui.prevAx,'fontweight','bold','color',borderColor,'tag','gSROI','verticalalignment','top');
    ROIlist(i).name = s{i};
    ROIlist(i).color = borderColor;
    ROIlist(i).notes = get(ui.notes,'string');
    ROIlist(i).dispname = get(ui.dispname,'value');
    ROIlist(i).capture = get(ui.capture,'value');
    set(ROIlist(i).lines,'linewidth',1,'userdata',i,'color',borderColor,'tag','gSROI')
    
    if get(ui.dispname,'value')==0
        set(ROIlist(i).txt,'visible','off')
    end
    
    if get(ui.capture,'value')==0 | get(gui.dispSROI,'value')==0
        set(ROIlist(i).lines,'visible','off');
        set(ROIlist(i).txt,'visible','off');
    end
   
end
set(gui.listSROI,'string',s,'userdata',ROIlist)
if isempty(ROIlist); return; end
val = get(gui.listSROI,'value');
if (val>length(ROIlist)); return; end
set(ROIlist(val).lines,'linewidth',2)

function dragSROI(obj,event,fig,ind,side)
gui = get(fig,'userdata');
ROIlist = get(gui.listSROI,'userdata');
for i = 1:length(ROIlist)
    set(ROIlist(i).lines ,'linewidth',1)
end
set(gui.listSROI,'value',ind)
set(ROIlist(ind).lines,'linewidth',2)
set(gui.sROIControls,'visible','off')
set(gui.sROIControls(ind),'visible','on')
set(gui.vidFig,'windowbuttonmotionfcn',{@moveSROI,fig,obj,ind,side})
set(gui.vidFig,'windowbuttonupfcn','set(gcbf,''windowbuttonmotionfcn'','''')');
function moveSROI(obj,event,fig,hline,ind,side)
gui = get(fig,'userdata');
pos = get(gca,'currentpoint');
gROI = gui.vi.ROIPosition;
ROIlist = get(gui.listSROI,'userdata');
ROI = ROIlist(ind).ROI;

switch side
    case 3 %Bottom
        ROIlist(ind).ROI(4) = round(pos(1,2))-ROI(2);
        if ROIlist(ind).ROI(4)<1; ROIlist(ind).ROI(4) = 1; end
    case 2 %Right
        ROIlist(ind).ROI(3) = round(pos(1,1))-ROI(1);
        if ROIlist(ind).ROI(3)<1; ROIlist(ind).ROI(3) = 1; end
    case {1,4} %Top or left
        ROIlist(ind).ROI(2) = round(pos(1,2));
        ROIlist(ind).ROI(1) = round(pos(1,1));
end
ROI = ROIlist(ind).ROI;
if ROI(1)+ROI(3) > gROI(3); return; end
if ROI(2)+ROI(4) > gROI(4); return; end
if round(pos(1,2))<1; return; end
if round(pos(1,1))<1; return; end


set(ROIlist(ind).lines(1),'xdata',[0 ROI(3)]+ROI(1),'ydata',[0 0]+ROI(2))
set(ROIlist(ind).lines(2),'xdata',[ROI(3) ROI(3)]+ROI(1),'ydata',[0 ROI(4)]+ROI(2))
set(ROIlist(ind).lines(3),'xdata',[ROI(3) 0]+ROI(1),'ydata',[ROI(4) ROI(4)]+ROI(2))
set(ROIlist(ind).lines(4),'xdata',[0 0]+ROI(1),'ydata',[ROI(4) 0]+ROI(2))
set(ROIlist(ind).txt,'position',[ROI(1) ,ROI(2)])

ui = get(gui.sROIControls(ind),'userdata');
set(ui.x,'string',num2str(ROI(1)))
set(ui.y,'string',num2str(ROI(2)))
set(ui.w,'string',num2str(ROI(3)))
set(ui.h,'string',num2str(ROI(4)))


set(gui.listSROI,'userdata',ROIlist)


function remSROI(obj,event,fig)
gui = get(fig,'userdata');
ROIlist = get(gui.listSROI,'userdata');
if isempty(ROIlist); return; end

ind = get(gui.listSROI,'value');
delete(gui.sROIControls(ind))
gui.sROIControls(ind)=[];
set(gui.fig,'userdata',gui)

ROIlist(ind) = [];
drawROI(ROIlist,fig);
if length(gui.sROIControls)<ind
    set(gui.listSROI,'value',ind-1)
else
    set(gui.listSROI,'value',ind)
end
    
dispSROI(gui.dispSROI,'',fig)

if ~isempty(ROIlist)
    try set(ROIlist(ind).lines,'linewidth',2); end
    set([gui.gROI gui.gROIdraw gui.gROIrst],'enable','off')
    ind = get(gui.listSROI,'value');
    set(gui.sROIControls(ind),'visible','on');
else
    set([gui.gROI gui.gROIdraw gui.gROIrst],'enable','on')
    set(gui.listSROI,'value',1)
    set(gui.sROIpopDown,'enable','off')
end


function copyNameSROI(obj,event,fig)
gui = get(fig,'userdata');
val = get(gui.listSROI,'value');
ROIlist = get(gui.listSROI,'userdata');

srNameRoot = get(gui.sROInameRoot,'string');
srID = str2double(get(gui.sROInameNum,'string'));
if isnan(srID); return; end

%Update sROI name
ui = get(gui.sROIControls(val),'userdata');
set(ui.name,'string',sprintf('%s_%02d',srNameRoot,srID))

%Update Color
set(ui.color,'backgroundcolor',get(gui.sROIcolor,'backgroundcolor'));

%itterate index
set(gui.sROInameNum,'string',num2str(srID+1))

%Step to next entry in the list
if val~=length(ROIlist)
    val = val+1;
    set(gui.listSROI,'value',val)
end

%Display current control box
set(gui.sROIControls,'visible','off')
set(gui.sROIControls(val),'visible','on')

drawROI(ROIlist,fig)


function dispSROI(obj,event,fig)
gui = get(fig,'userdata');
ROIlist = get(gui.listSROI,'userdata');

drawROI(ROIlist,fig)
   

function listSROI(obj,event,fig)
gui = get(fig,'userdata');
ROIlist = get(gui.listSROI,'userdata');
if isempty(ROIlist); return; end
ind = get(gui.listSROI,'value');

drawROI(ROIlist,fig)

set(gui.sROIControls,'visible','off')
set(gui.sROIControls(ind),'visible','on')

function colorSROI(obj,event,fig)
gui = get(fig,'userdata');
ROIlist = get(gui.listSROI,'userdata');
c=uisetcolor('Select a ROI Color');
if length(c)~=3; return; end

set(obj,'backgroundcolor',c)

drawROI(ROIlist,fig)







function setMetric(obj,event,fig)
gui = get(fig,'userdata');

delete(timerfind('tag','gHist'))
delete(timerfind('tag','gLineScan'))

set(gui.histLines,'visible','off')
set([gui.pxH gui.pxHt gui.pxV gui.pxVt],'visible','off')
set([gui.pxH gui.pxV],'linewidth',2)
set(gui.lineScanTrace,'visible','off')

set(gui.HistPanel,'visible','on')
set(gui.MetaDataPanel,'visible','off')

if get(obj,'value')==0
    return
end
set(gui.histBars,'visible','off')
set([gui.histActivate gui.xHair gui.MetaData],'value',0)
switch obj
    %Histogram
    case gui.histActivate
        set(gui.HistPanel,'visible','on')
        set(gui.histActivate,'value',1)
        set(gui.histLines,'visible','on')
        set(gui.histBars,'visible','on')
        set(gui.histAx,'xlim',[-5 260])
        
        hist = get(gui.histLines(1),'userdata');        
        t = timer('tag','gHist','busymode','drop','period',.1,'TasksToExecute',inf,...
            'ExecutionMode','FixedRate','timerfcn',{@histUpdate,fig,hist.lines(1)});
        start(t)
    %Metadata
    case gui.MetaData
        set(gui.HistPanel,'visible','off')
        set(gui.MetaDataPanel,'visible','on')
        set(gui.MetaData,'value',1)
        %Display metadata Text panel
        
	%Cross Hairs and line scan
    case gui.xHair
        set(gui.HistPanel,'visible','on')
        set(gui.xHair,'value',1)
        set(gui.lineScanTrace,'visible','on')
        set([gui.pxH gui.pxHt gui.pxV gui.pxVt],'visible','on')
        ROI = get(gui.vi,'ROIPosition');
        set([gui.pxH gui.pxHt gui.pxV gui.pxVt],'visible','on')
        set(gui.pxH,'ydata',[1 1]*round(ROI(4)/2));
        set(gui.pxHt,'string',num2str(round(ROI(4)/2)))
        set(gui.pxHt,'position',[1 ROI(4)/2+13])
        set(gui.pxV,'xdata',[1 1]*round(ROI(3)/2));
        set(gui.pxVt,'string',num2str(round(ROI(3)/2)))
        set(gui.pxVt,'position',[ROI(3)/2+13 13])
        
        set(gui.pxH,'linewidth',3)
       
        t = timer('tag','gLineScan','busymode','drop','period',.1,'TasksToExecute',inf,...
            'ExecutionMode','FixedRate','timerfcn',{@lineScanUpdate,fig});
        start(t)
end

%Functions for managing the line scan and alignment crosshairs
function lineScanUpdate(obj,event,fig)
gui = get(fig,'userdata');
im = get(gui.pImage,'cdata');
switch findobj([gui.pxH,gui.pxV],'linewidth',3)
    case gui.pxH
        ind = min(get(gui.pxH,'ydata'));
        ydat = mean(im(ind,:,:),3);
    case gui.pxV
        ind = min(get(gui.pxV,'xdata'));
        ydat = mean(im(:,ind,:),3);
end

xdat = 1:length(ydat);
set(gui.lineScanTrace,'xdata',xdat,'ydata',ydat)
set(gui.histAx,'xlim',[1 length(ydat)],'ylim',[-10 270])
function dragCross(obj,event,fig)
gui=get(fig,'userdata');
gROI = gui.vi.ROIPosition;
set([gui.pxH gui.pxV],'linewidth',2)
switch obj
    case gui.pxH
        set(gui.vidFig,'windowbuttonmotionfcn',{@moveCross,obj,0,gui.pxHt,gROI})
        set(gui.pxH,'linewidth',3)
    case gui.pxV
        set(gui.vidFig,'windowbuttonmotionfcn',{@moveCross,obj,1,gui.pxVt,gROI})
        set(gui.pxV,'linewidth',3)
end
set(gui.vidFig,'windowbuttonupfcn','set(gcbf,''windowbuttonmotionfcn'','''')');
function moveCross(obj,event,pTrace,ID,txt,gROI)
pos = get(gca,'currentpoint');
if pos(1,2)<=0.5||pos(1,1)<=0.5; return; end
if round(pos(1,2))>gROI(4)||round(pos(1,1))>gROI(3); return; end

switch ID
    case 0
        set(pTrace,'ydata',[1 1]*round(pos(1,2)));
        set(txt,'position',[1 pos(1,2)+13],'string',num2str(round(pos(1,2))))
    case 1
        set(pTrace,'xdata',[1 1]*round(pos(1,1)));
        set(txt,'position',[pos(1,1)+5 13],'string',num2str(round(pos(1,1))))
end

%histogram functions
function histUpdate(obj,event,fig,hline)
gui = get(fig,'userdata');
ghist = get(hline,'userdata');
ROI = ghist.ROI;
im = get(gui.pImage,'cdata');
im = im(ROI(2):(ROI(2)+ROI(4)),ROI(1):(ROI(1)+ROI(3)),:);
if gui.vi.NumberOfBands == 1
    n = hist(im(:),0:255);
    for i = 1:256
        set(gui.histBars(i),'ydata',[0 0 1 1]*n(i))
    end
    set(gui.histAx,'ylim',[0 1.05*max(n)])
else %3 chans
    maxval = 0;
    for i = 1:3
        temp = im(:,:,i);
        n = hist(temp(:),0:255);
        for j = 1:256
            set(gui.histBars(i,j),'ydata',[0 0 1 1]*n(j))
        end   
        maxval = max([maxval,max(n)]);
    end
    set(gui.histAx,'ylim',[0 1.05*maxval])
end
set(gui.histAx,'xscale','linear','xlim',[-5 260])
drawnow

%pop out for hist or line scan
function histPop(obj,event,fig)
gui = get(fig,'userdata');
if get(gui.histActivate,'value')==1
% if strcmpi(get(gui.histBars(1),'visible'),'on')
    tempFig = figure('numbertitle','off','name','Camera Histogram');
    ax = copyobj(gui.histAx,tempFig);
    set(ax,'position',[.1 .1 .8 .8],'xtick',0:20:255,'ytickmode','auto','xlim',[-5 260],'xscale','linear')
    xlabel('Pixel Intensity')
    ylabel('Number of Pixels')
else
    tempFig = figure('numbertitle','off','name','Camera Histogram');
    ax = copyobj(gui.histAx,tempFig);
    set(ax,'position',[.1 .1 .8 .8],'xtick',0:20:255,'ytickmode','auto','xlim',[-5 260],'xscale','linear')
    ylabel('Pixel Intensity')
    xlabel('Pixel Coordinate')
end



function dragBox(obj,event,fig,side)
gui = get(fig,'userdata');

set(gui.vidFig,'windowbuttonmotionfcn',{@moveBox,fig,obj,side})
set(gui.vidFig,'windowbuttonupfcn','set(gcbf,''windowbuttonmotionfcn'','''')');
function moveBox(obj,event,fig,hline,side)
gui = get(fig,'userdata');
pos = get(gca,'currentpoint');
gROI = gui.vi.ROIPosition;
box = get(hline,'userdata');
ROI = box.ROI;

switch side
    case 3 %Bottom
        box.ROI(4) = round(pos(1,2))-ROI(2);
        if box.ROI(4)<1; box.ROI(4) = 1; end
    case 2 %Right
        box.ROI(3) = round(pos(1,1))-ROI(1);
        if box.ROI(3)<1; box.ROI(3) = 1; end
    case {1,4} %Top & Left
        box.ROI(2) = round(pos(1,2));
        box.ROI(1) = round(pos(1,1));
end
ROI = box.ROI;
if ROI(1)+ROI(3) > gROI(3); return; end
if ROI(2)+ROI(4) > gROI(4); return; end
if pos(1,2)<1; return; end
if pos(1,1)<1; return; end


set(box.lines(1),'xdata',[0 ROI(3)]+ROI(1),'ydata',[0 0]+ROI(2))
set(box.lines(2),'xdata',[ROI(3) ROI(3)]+ROI(1),'ydata',[0 ROI(4)]+ROI(2))
set(box.lines(3),'xdata',[ROI(3) 0]+ROI(1),'ydata',[ROI(4) ROI(4)]+ROI(2))
set(box.lines(4),'xdata',[0 0]+ROI(1),'ydata',[ROI(4) 0]+ROI(2))
try set(box.txt,'position',[ROI(1) ,ROI(2)+ROI(4)/2]); end

set(box.lines,'userdata',box)


function setCameraProp(obj,event,fig)
gui = get(fig,'userdata');
data = get(gui.Params,'data');
ss = getselectedsource(gui.vi);

for i = 1:size(data,1)
    try
        set(ss,data{i,1},data{i,2})
    end
end
updatePropStat(fig)

function updatePropStat(fig)
gui = get(fig,'userdata');
ss = getselectedsource(gui.vi);
a = get(ss);
c = fieldnames(a);
b=imaqhwinfo(gui.vi);
    

%fill status window
%DeviceID, FramesAcquired, videoformat, returnedcolorspace, ss-UniqueID,
data{1,1}='Device ID'; data{1,2}=gui.vi.DeviceID;
data{2,1}='Video Format'; data{2,2}=gui.vi.VideoFormat;
data{3,1}='Data Type'; data{3,2}=b.NativeDataType;
data{4,1}='Color Space'; data{4,2}=gui.vi.ReturnedColorSpace;
try data{5,1}='Unique ID'; data{5,2}=ss.UniqueID; end
data{end+1,1}='Adaptor'; data{end,2}=gui.adaptor;
set(gui.Status,'data',data)

if ~isfield(a,'FrameRate')
    set([gui.FPStxt gui.FPS],'visible','on');
end

data={}; j = 1;
%ignore properties: parent, selected, tag, type, frameTimeout,
for i = 1:length(c)
    if strcmpi(c{i},'parent')|strcmpi(c{i},'selected')|strcmpi(c{i},'tag')|strcmpi(c{i},'type')|strcmpi(c{i},'NormalizedBytesPerPacket')|strcmpi(c{i},'FrameTimeout')
        continue
    end
    data{j,1} = c{i}; data{j,2} = eval(['a.',c{i}]); 
    
    pinfo=propinfo(ss,c{i});
    
    switch pinfo.Constraint
        case 'none'
            data{j,3} = pinfo.ConstraintValue;
        case 'bounded'
            tmp =  pinfo.ConstraintValue;
            data{j,3} = ['  [',num2str(tmp(1)),' ',num2str(tmp(2)),']'];
        case 'enum'
            str = '';
            for i = 1:length(pinfo.ConstraintValue)
                str = [str,', ',pinfo.ConstraintValue{i}];
            end
            str(1:2) = ' ';
            data{j,3} = str;
    end
    j = j+1;
end
set(gui.Params,'data',data)


function parsePlugs(fig)
gui = get(fig,'userdata');
%Set trigger properties
data={};
classnames={};

plugs = what('plugins');
plugs = plugs.m;

%Scan for available plugins
for i=1:length(plugs)
    %These files are not gVision Modes
    if (strcmpi(plugs{i},'text2im.m')); continue; end
    if (strcmpi(plugs{i},'gvPlugin.m')); continue; end    
    
    switch plugs{i}
        case 'gvBasic.m' %The default mode
            [cName dispName description] = plugins.gvBasic.gvPlugName;
            data{end+1,1} = true;
            data{end,2} = dispName;
            data{end,3} = description;
            classnames{end+1} = cName;
        otherwise
            %Other Modes, including user generated modes
            [p cl ext] = fileparts( plugs{i} );
            eval(sprintf('[cName dispName description] = plugins.%s.gvPlugName;',cl))
            data{end+1,1} = false;
            data{end,2} = dispName;
            data{end,3} = description;
            classnames{end+1} = cName;
            
    end
end
set(gui.ModeOptions,'data',data,'userdata',classnames)


    

function makegui(ver)

imaqreset

delete(findobj('tag','gfh08'))
delete(findobj('tag','gVIprev'))
delete(timerfind('tag','gHist'))
gui.ver = ver;
gui.mode = [];
gui.prevAx = [];

%SPLASH SCREEN
tempFig=figure('position',[0 0 720 217],'menubar','none','numbertitle','off','name',...
    ['gVision ',ver],'tag','gfh08','resize','off','color',[.6 .8 .6]);
centerfig(tempFig)
axes('position',[0 0 1 1],'color',[.8 .9 .8],'yticklabel',[]...
    ,'xticklabel',[],'xlim',[0 1],'ylim',[0 1],'xcolor',[.4 .4 .4],'ycolor',[.4 .4 .4],'visible','off');
aTemp=imread('gVisionSplash.png');
image(aTemp)
axis  off
pause(0.1); 
%END SPLASH SCREEN

gui.fig=figure('position',[0 0 800 600],'numbertitle','off','menubar','none',...
    'name',['gVision ',ver],'tag','gfh08','resize','off');
centerfig(gui.fig)
delete(tempFig)
set(gui.fig,'closerequestfcn',{@programQuit,gui.fig});


%Create Pull Down Menus
gui.menu(1)=uimenu('Label','&File','userdata',pwd);
gui.fileMen(1) = uimenu(gui.menu(1),'Label','L&oad State','callback',{@stateLoad,gui.fig});
gui.fileMen(2) = uimenu(gui.menu(1),'Label','&Save State','callback',{@stateSave,gui.fig});


gui.menu(2)=uimenu('Label','&Camera');
gui.menu(3)=uimenu('Label','Stream &Options');
gui.streamMen(1) = uimenu(gui.menu(3),'Label','Preview Colormap','enable','off');
    gui.cmapMen(1) = uimenu(gui.streamMen(1),'Label','Default','checked','on','callback',{@changePreviewCMap,gui.fig});
    gui.cmapMen(2) = uimenu(gui.streamMen(1),'Label','Hot','checked','off','callback',{@changePreviewCMap,gui.fig});
    gui.cmapMen(3) = uimenu(gui.streamMen(1),'Label','Cool','checked','off','callback',{@changePreviewCMap,gui.fig});
    gui.cmapMen(4) = uimenu(gui.streamMen(1),'Label','HSV','checked','off','callback',{@changePreviewCMap,gui.fig});
    gui.cmapMen(5) = uimenu(gui.streamMen(1),'Label','Jet','checked','off','callback',{@changePreviewCMap,gui.fig});
gui.streamMen(2) = uimenu(gui.menu(3),'Label','Iterate File Names','checked','on','callback',@toggleChecked);
gui.streamMen(3) = uimenu(gui.menu(3),'Label','Time Stamp Overlay','checked','off','callback',@toggleChecked);
gui.streamMen(4) = uimenu(gui.menu(3),'Label','Force Color to Gray','checked','off','callback',@toggleChecked);
gui.streamMen(5) = uimenu(gui.menu(3),'Label','Compression');
    profiles = VideoWriter.getProfiles();
    for i=1:length(profiles)
        menu(i) = uimenu(gui.streamMen(5),'Label',profiles(i).Name,'checked','off');
        if strcmpi(profiles(i).Name,'Uncompressed AVI')
            set(menu(i),'checked','on')
        end
    end
    set(menu,'callback',{@menuMutex,gui.streamMen(5)})
gui.menu(4)=uimenu('Label','&About','enable','off');
delete(imaqfind)

%Scan Available Cameras
gTEMP=imaqhwinfo;
for i=1:length(gTEMP.InstalledAdaptors)
    
    caminfo=imaqhwinfo(gTEMP.InstalledAdaptors{i});
    if ~isempty(caminfo.DeviceIDs); 
        gT2=uimenu('parent',gui.menu(2),'label',caminfo.AdaptorName);
    end
    for j=1:length(caminfo.DeviceIDs)
        camformat=imaqhwinfo(gTEMP.InstalledAdaptors{i},caminfo.DeviceIDs{j});
        str = [num2str(caminfo.DeviceIDs{j}),': ',camformat.DeviceName];
        
        try
            vi = videoinput(gTEMP.InstalledAdaptors{i},caminfo.DeviceIDs{j});
            str = [str,' ',' - ID:',get(vi.Source,'uniqueID')];
            delete(vi)
        end
        
        gT3=uimenu('parent',gT2,'label',str,'userdata',caminfo.DeviceIDs{j});
        
        for jj=length(camformat.SupportedFormats):-1:1
            uimenu('parent',gT3,'callback',{@addCam,gui.fig,0},'label',...
                camformat.SupportedFormats{jj})
        end
    end
end



%Camera Controls and Status Info
gui.Status=uitable('units','pixels','position',[10 470 400 120]);
cnames = {'Parameter','Value'}; ceditable = [false false];
set(gui.Status,'columnName',cnames,'columnEditable',ceditable)
gui.Params=uitable('units','pixels','position',[10 300 400 160]);
cnames={'Parameter','Value','Range'}; ceditable=[false true false];
set(gui.Params,'ColumnName',cnames,'ColumnEditable',ceditable,'CellEditCallback',{@setCameraProp,gui.fig})

%ROI Controls for Global and sub-ROI
gui.ROIPanel=uipanel('Title','Regions of Interest','units','pixels','position',[10 5 400 280],'backgroundcolor',get(gui.fig,'color'),'fontweight','bold','foregroundcolor',[.5 .3 .3]);
%Global
uicontrol(gui.ROIPanel,'style','text','string','Global ROI:','units','normalized','position',[.05 .85 .25 .07],'horizontalalignment','left','backgroundcolor',get(gui.ROIPanel,'backgroundcolor'),'fontweight','bold')
gui.gROI=uicontrol(gui.ROIPanel,'style','edit','string','','units','normalized','position',[.25 .85 .4 .1],'backgroundcolor','w',...
    'callback',{@setGROI,gui.fig},'enable','off');
gui.gROIdraw = uicontrol(gui.ROIPanel,'style','toggle','string','Draw','units','normalized','position',[.68 .86 .12 .1],...
    'callback',{@drawGROI,gui.fig},'enable','off');
gui.gROIrst = uicontrol(gui.ROIPanel,'style','pushbutton','string','Reset','units','normalized','position',[.83 .86 .12 .1],...
    'callback',{@rstGROI,gui.fig},'enable','off');
uicontrol(gui.ROIPanel,'style','text','backgroundcolor',[.5 .3 .3],'units','normalized','position',[.05 .8 .9 .02]);

%software sub-roi
uicontrol(gui.ROIPanel,'style','text','string','Sub ROI:','units','normalized','position',[.01 .65 .25 .1],'horizontalalignment','left','backgroundcolor',get(gui.ROIPanel,'backgroundcolor'),'fontweight','bold')
gui.addSROI=uicontrol(gui.ROIPanel,'style','pushbutton','string','Add','units','normalized','position',[.01 .6 .2 .07],...
    'callback',{@addSROI,gui.fig},'enable','off');
gui.remSROI=uicontrol(gui.ROIPanel,'style','pushbutton','string','Remove','units','normalized','position',[.01 .52 .2 .07],...
    'callback',{@remSROI,gui.fig},'enable','off');

gui.dispSROI=uicontrol(gui.ROIPanel,'style','checkbox','string','Display All','units','normalized','position',[.01 .01 .2 .07],...
    'callback',{@dispSROI,gui.fig},'enable','off','backgroundcolor',get(gui.ROIPanel,'backgroundcolor'),'value',1);

%Naming support, name root, itterator.  Will allow user to specify a root
%and an itteration number (%02d) and copies color, name, and increments
%itterator by 1 and steps down to next existing sROI
gui.sROInameRoot = uicontrol(gui.ROIPanel,'style','Edit','string','NameRoot','units','normalized','position',[.01 .31 .2 .07],...
    'callback','','enable','off','backgroundcolor','w','tooltipstring','Name Root to Propagate in sub-ROI List');
gui.sROInameNum = uicontrol(gui.ROIPanel,'style','Edit','string','0','units','normalized','position',[.01 .22 .14 .07],...
    'callback','','enable','off','backgroundcolor','w','tooltipstring','Next Number to Label sub-ROI name with');
gui.sROIcolor=uicontrol(gui.ROIPanel,'style','pushbutton','string','','units','normalized','position',[.16 .22 .05 .07],...
    'callback',{@colorSROI,gui.fig},'enable','off','backgroundcolor',[.7 0 0],'tooltipstring','Color to Propogate in sub-ROI list');
gui.sROIpopDown=uicontrol(gui.ROIPanel,'style','pushbutton','string','Pop Down \/','units','normalized','position',[.01 .13 .2 .07],...
    'callback',{@copyNameSROI,gui.fig},'enable','off');


%Main list for sub-regions of interest (i.e. channels) in the global video
gui.listSROI=uicontrol(gui.ROIPanel,'style','listbox','string','','units','normalized','position',[.22 .01 .3 .75],...
    'callback',{@listSROI,gui.fig},'backgroundcolor','w');

gui.sROIControls = [];

%Histogram, Focus, and Alignment Controls
gui.HistPanel=uipanel('Title','Histogram & Image Alignment','units','pixels','position',[420 32 370 168],'backgroundcolor',get(gui.fig,'color'),'fontweight','bold','visible','off');
gui.histAx = axes('parent',gui.HistPanel,'position',[0.01 .01 .98 .98]);
set(gui.histAx,'xtick',[],'ytick',[],'xlim',[0 255],'ylim',[0 10^3],'box','on')
gui.histPop = uicontrol(gui.HistPanel,'style','pushbutton','units','normalized','position',[0.01 .9 .04 .10],'backgroundcolor',[.7 0 0],'callback',{@histPop,gui.fig});

%MetaData Panel
gui.MetaDataPanel=uipanel('Title','Metadata','units','pixels','position',[420 32 370 168],'backgroundcolor',get(gui.fig,'color'),'fontweight','bold','visible','on');
gui.MetaDataChk = uicontrol(gui.MetaDataPanel,'units','normalized','position',[0.01 0.9 0.49 0.1],...
                    'style','check','value',1,'string','Log MetaData with Recording','backgroundcolor',get(gui.MetaDataPanel,'backgroundcolor'));
gui.MetaDataClear = uicontrol(gui.MetaDataPanel,'units','normalized','position',[0.5 0.9 0.49 0.1],...
                    'style','check','value',0,'string','Clear Metadata/Notes on Stop','backgroundcolor',get(gui.MetaDataPanel,'backgroundcolor'));
gui.MetaDataTxt = uicontrol(gui.MetaDataPanel,'units','normalized','position',[0.01 0.01 0.98 0.89],...
                    'style','edit','max',2,'backgroundcolor','w','horizontalalignment','left');

gui.histActivate = uicontrol(gui.fig,'style','toggle','string','Histogram','units','pixels',...
    'position',[420 5 70 25],'backgroundcolor',[.7 .7 .8],'fontweight','bold','enable','off','callback',{@setMetric,gui.fig});
gui.xHair = uicontrol(gui.fig,'style','toggle','string','Alignment Crosshairs','units','pixels',...
    'position',[492 5 150 25],'backgroundcolor',[.8 .7 .7],'fontweight','bold','enable','off','callback',{@setMetric,gui.fig});
gui.MetaData = uicontrol(gui.fig,'style','toggle','string','Video Metadata','units','pixels',...
    'position',[644 5 146 25],'backgroundcolor',[.7 .8 .7],'fontweight','bold','enable','off','value',1,'callback',{@setMetric,gui.fig});
    
%Mode Controls
gui.ModePanel=uipanel('Title','gVision Modes','units','pixels','position',[420 210 370 120],'backgroundcolor',get(gui.fig,'color'),'fontweight','bold','foregroundcolor','b');
gui.ModeOptions = uitable('units','normalized','parent',gui.ModePanel,'position',[0.01 .01 .98 .98]);
cnames = {' ','Name','Description'}; ceditable = [true false false];
set(gui.ModeOptions,'columnName',cnames,'columnEditable',ceditable,'CellEditCallback',{@ModeChange,gui.fig})

%Recording Controls
gui.RecPanel=uipanel('Title','Recording Controls','units','pixels','position',[420 330 370 270],'backgroundcolor',get(gui.fig,'color'),'fontweight','bold','foregroundcolor',[.8 0 0]);
gui.previewONOFF = uicontrol(gui.RecPanel,'style','toggle','string','Preview','units','normalized',...
    'position',[0.03 0.89 0.3 0.1],'backgroundcolor',[.8 .7 .7],'fontweight','bold','enable','off','callback',{@switchPreview,gui.fig});
gui.snapShot = uicontrol(gui.RecPanel,'style','pushbutton','string','Snap Shot','units','normalized',...
    'position',[0.35 0.89 0.3 0.1],'backgroundcolor',[.7 .8 .7],'fontweight','bold','enable','off','callback',{@snapShot,gui.fig});
gui.simpleAdvanced = uicontrol(gui.RecPanel,'style','checkbox','string','Advanced','units','normalized',...
    'position',[0.67 0.89 0.3 0.1],'backgroundcolor',get(gui.RecPanel,'backgroundcolor'),'fontweight','bold','value',0,'enable','on','callback',{@simpleAdvanced,gui.fig});

%Default File Name
file.fname = 'temp.avi';
file.pname = pwd;
gui.fName = uicontrol(gui.RecPanel,'style','text','string','temp.avi','units','normalized',...
    'position',[0.1 0.8 0.8 0.06],'backgroundcolor',[.8 .8 .8]-.1,'fontweight','bold','userdata',file);
gui.fileSelect = uicontrol(gui.RecPanel,'style','pushbutton','string','Browse...','units','normalized',...
    'position',[0.1 0.67 0.35 0.08],'fontweight','bold','callback',{@setRecFile,gui.fig},'enable','off');
gui.OpenDir = uicontrol(gui.RecPanel,'style','pushbutton','string','Open Directory','units','normalized',...
    'position',[0.55 0.67 0.35 0.08],'fontweight','bold','callback',{@openDir,gui.fig});

uicontrol(gui.RecPanel,'style','text','backgroundcolor',[.5 .3 .3],'units','normalized','position',[.05 .64 .9 .01]);

%FramesPerTrigger
uicontrol(gui.RecPanel,'style','text','string','Frames/Trigger:','units','normalized','position',[.05 .29 .3 .06],'horizontalalignment','left',...
    'backgroundcolor',get(gui.RecPanel,'backgroundcolor'),'fontweight','bold')
gui.FPT = uicontrol(gui.RecPanel,'style','edit','string','inf','units','normalized','backgroundcolor','w',...
    'position',[0.32 0.28 0.2 0.08],'callback',{@setRec,gui.fig},'enable','off','userdata',inf);
%TriggerRepeat
uicontrol(gui.RecPanel,'style','text','string','Trigger Repeat:','units','normalized','position',[.05 .19 .3 .06],'horizontalalignment','left',...
    'backgroundcolor',get(gui.RecPanel,'backgroundcolor'),'fontweight','bold')
gui.TR = uicontrol(gui.RecPanel,'style','edit','string','0','units','normalized','backgroundcolor','w',...
    'position',[0.32 0.18 0.2 0.08],'callback',{@setRec,gui.fig},'enable','off','userdata',0);

uicontrol(gui.RecPanel,'style','text','backgroundcolor',[.5 .3 .3],'units','normalized','position',[.05 .15 .9 .01]);

%Frame Rate for file write(if needed)
gui.FPStxt = uicontrol(gui.RecPanel,'style','text','string','File Frame Rate:','units','normalized','position',[.05 .05 .3 .06],'horizontalalignment','left',...
    'backgroundcolor',get(gui.RecPanel,'backgroundcolor'),'fontweight','bold','visible','off');
gui.FPS = uicontrol(gui.RecPanel,'style','edit','string','30','units','normalized','backgroundcolor','w',...
    'position',[0.32 0.05 0.2 0.08],'callback',{@setRec,gui.fig},'visible','off','userdata',30);



%Status Text Window
gui.RecStatus = uicontrol(gui.RecPanel,'style','text','string',sprintf('nFrames: 0\neTime: \nfps: '),...
    'units','normalized','position',[.55 .19 .4 .26],'horizontalalignment','left','backgroundcolor',[1 .7 .7],'fontweight','bold','userdata',0);
gui.StartStop = uicontrol(gui.RecPanel,'style','toggle','string','Start','units','normalized',...
    'position',[0.6 0.01 0.35 0.1],'backgroundcolor',[.9 .9 .2],'fontweight','bold','fontsize',12,'callback',{@StartStop,gui.fig},'enable','off');

set(gui.fig,'userdata',gui)
simpleAdvanced([],[],gui.fig) %enter into simple mode to start with

%Seek out and parse a defaultGvState.mat 
a = dir;
for i = 1:length(a)
    if strcmpi('defaultGvState.mat',a(i).name)
        fprintf('Loading Default gVision State saved in defaultGvState.mat\n')
        stateLoad([],fullfile(pwd,'defaultGvState.mat'),gui.fig)
        break
    end
end

%Switch between simple and advanced mode
% This displays more or less user interface options (like camera property
% controls) and collapses the UI
function simpleAdvanced(obj,event,fig)
gui = get(fig,'userdata');

switch get(gui.simpleAdvanced,'value');
    case 1 %Advanced Mode
        %Standard UI mode with all displays
        set(gui.RecPanel,'position',[420 330 370 270])
        set(gui.MetaDataPanel,'position',[420 32 370 168]);
        set([gui.histActivate gui.xHair],'value',0); set(gui.MetaData,'value',1); try setMetric(gui.MetaData,'',fig); end
        set([gui.ModePanel gui.ROIPanel gui.Status gui.Params],'visible','on')
        set(gui.fig,'position',[0 0 800 600]);
        centerfig(gui.fig)
    case 0 %Simple Mode,
        %RecPanel & MetaDataPanel only, shrink window
        set(gui.RecPanel,'position',[1 170 370 270])
        set(gui.MetaDataPanel,'position',[1 1 370 168],'visible','on')
        set([gui.ModePanel gui.ROIPanel gui.Status gui.Params gui.HistPanel],'visible','off')
        set(gui.fig,'position',[0 0 370 270+170]);
        centerfig(gui.fig)
end


%If simply closing the preview figure, stop previewing
function prevFigClose(obj,event,fig)
gui = get(fig,'userdata');
set(gui.previewONOFF,'value',0)
switchPreview(gui.previewONOFF,'',fig);


function programQuit(~,~,fig)
gui = get(fig,'userdata');
imaqreset
delete(findobj('tag','gfh08'))
delete(findobj('tag','gVIprev'))
delete(timerfind('tag','gHist'))


function openDir(obj,event,fig)
gui = get(fig,'userdata');
dos(['explorer ',get(gui.menu(1),'userdata')]);

%Present the user with the file selection dialog for file name and format
function setRecFile(obj,event,fig)
gui = get(fig,'userdata');

if get(gui.vi, 'NumberOfBands')==1
    foptions = {'*.avi','AVI files (*.avi)';...
                '*.fmf','Fly Movie Format (*.fmf)';...
                '*.tif','TIFF Stack 16-bit (*.tif)'};
else %Options for color video
    foptions = {'*.avi','AVI files (*.avi)';...
        '*.tif','TIFF Stack 8/16-bit (*.tif)'};
end

info = imaqhwinfo(gui.vi);
if strcmpi(info.NativeDataType,'uint16')
    foptions = {'*.tif','TIFF Stack 16-bit (*.tif)'};
    [fname, pname] = uiputfile(foptions,'Target Video File',fullfile(get(gui.menu(1),'userdata'),[datestr(now,30),'.tif']));
    if isequal(fname,0) || isequal(pname,0)
       return
    end
    set(gui.menu(1),'userdata',pname)

    file.fname = fname;
    file.pname = pname;
    set(gui.fName,'string',fullfile(pname,fname),'userdata',file);

else
    [fname, pname] = uiputfile(foptions,'Target Video File',fullfile(get(gui.menu(1),'userdata'),[datestr(now,30),'.avi']));
    if isequal(fname,0) || isequal(pname,0)
       return
    end
    set(gui.menu(1),'userdata',pname)

    file.fname = fname;
    file.pname = pname;
    set(gui.fName,'string',fname,'userdata',file);
end

%Load compression options for selected file type
[p n ext] = fileparts(fullfile(pname,fname));
delete(get(gui.streamMen(5),'children'))
switch ext
    case '.avi'
        profiles = VideoWriter.getProfiles();
        for i=1:length(profiles)
            menu(i) = uimenu(gui.streamMen(5),'Label',profiles(i).Name,'checked','off');
            if strcmpi(profiles(i).Name,'Uncompressed AVI')
                set(menu(i),'checked','on')
            end
        end
    case '.tif'
        profiles = Tiff.Compression;
        s = fieldnames(profiles);
        for i = 1:length(s)
            menu(i) = uimenu(gui.streamMen(5),'Label',s{i},'checked','off',...
                'userdata',getfield(profiles,s{i}));
            if strcmpi(s{i},'None')
                set(menu(i),'checked','on')
            end
        end       
    case '.fmf'
        menu(1) = uimenu(gui.streamMen(5),'Label','None','checked','on');
end

set(menu,'callback',{@menuMutex,gui.streamMen(5)})

%Mutual Exclusion for Compression Options
function menuMutex(obj,event,parent)
set(get(get(obj,'parent'),'children'),'checked','off')
set(obj,'checked','on');

function toggleChecked(obj,event)
if strcmpi(get(obj,'checked'),'on')
    set(obj,'checked','off')
else
    set(obj,'checked','on')
end



%State Saving/Loading
function stateSave(~,~,fig)
gui = get(fig,'userdata');

state.MetaDataTxt = get(gui.MetaDataTxt,'string');
state.MetaDataClear = get(gui.MetaDataClear,'value');
state.MetaDataChk = get(gui.MetaDataChk,'value');

data = get(gui.Status,'data');
state.adaptor = gui.adaptor;%
state.id = data{1,2};%
state.format = data{2,2};%

data = get(gui.Params,'data');%
state.params = data;%


state.sROIlist = get(gui.listSROI,'userdata');  %contains all data for sROI recreation

data = get(gui.ModeOptions,'data');
state.trig = data;

state.pwd = get(gui.menu(1),'userdata');
state.fname = get(gui.fName,'string');
state.prev = get(gui.previewONOFF,'value');
state.gROI = str2num(get(gui.gROI,'string'));%
state.FPT = get(gui.FPT,'userdata');%
state.TR = get(gui.TR,'userdata');%
state.FPS = get(gui.FPS,'userdata');%

state.pos = get(gui.fig,'position');
state.vpos = get(gui.vidFig,'position');

[fname, pname] = uiputfile('*.mat','Save Program State',fullfile(get(gui.menu(1),'userdata'),'defaultGvState.mat'));
if isequal(fname,0) || isequal(pname,0)
   return
end
set(gui.menu(1),'userdata',pname)

save(fullfile(pname,fname),'state')


%%%% STATE LOADING
function stateLoad(~,passfile,fig)
gui = get(fig,'userdata');


if isempty(passfile)
    [fname, pname] = uigetfile('*.mat','Load Program State',fullfile(get(gui.menu(1),'userdata'),'default.mat'));
    if isequal(fname,0) || isequal(pname,0)
        return
    end
else
    [pname,fname,ext] = fileparts(passfile);
    fname = [fname,ext];
end
set(gui.menu(1),'userdata',pname)
load(fullfile(pname,fname))
%state struct loaded

addCam(state,'',gui.fig,1)
gui = get(fig,'userdata');
ss = getselectedsource(gui.vi);
for i = 1:size(state.params,1);
    try set(ss,state.params{i,1},state.params{i,2}); end
end
set(gui.Params,'data',state.params);

gui.vi.ROIPosition = state.gROI;
ROI = gui.vi.ROIPosition;
res = gui.vi.VideoResolution;
set(gui.prevAx,'xlim',[1 res(1)]-ROI(1),'ylim',[1 res(2)]-ROI(2))
set(gui.gROI,'string',num2str(state.gROI),'userdata',state.gROI)

gui.vi.FramesPerTrigger = state.FPT;
set(gui.FPT,'string',num2str(state.FPT),'userdata',state.FPT)
set(gui.FPS,'string',num2str(state.FPS),'userdata',state.FPS)
gui.vi.TriggerRepeat = state.TR;
set(gui.TR,'string',num2str(state.TR),'userdata',state.TR)

for i = 1:size(state.trig,1)
    if state.trig{i,1}
        break
    end
end

%Create UIs for channels, update names and checkboxes and positions
for i = 1:length(state.sROIlist)
    addSROI([],'loading',gui.fig)  %Add default UIs
end
gui = get(gui.fig,'userdata'); %Grab pointers to new sROI uipanels
for i = 1:length(state.sROIlist)
    ui = get(gui.sROIControls(i),'userdata');
    set(ui.name,'string',state.sROIlist(i).name)
    set(ui.color,'backgroundcolor',state.sROIlist(i).color)
    set(ui.x,'string',num2str(state.sROIlist(i).ROI(1)))
    set(ui.y,'string',num2str(state.sROIlist(i).ROI(2)))
    set(ui.w,'string',num2str(state.sROIlist(i).ROI(3)))
    set(ui.h,'string',num2str(state.sROIlist(i).ROI(4)))
    set(ui.capture,'value',state.sROIlist(i).capture)
    set(ui.dispname,'value',state.sROIlist(i).dispname)
    set(ui.notes,'string',state.sROIlist(i).notes)
end

%Draw channels and populate list
if ~isempty(state.sROIlist)
    set(gui.listSROI,'value',1)
    drawROI(state.sROIlist,gui.fig);
end

%Reload metadata info
set(gui.MetaDataTxt,'string',state.MetaDataTxt);
set(gui.MetaDataClear,'value',state.MetaDataClear);
set(gui.MetaDataChk,'value',state.MetaDataChk);

% set(gui.fig,'position',state.pos)
% set(gui.vidFig,'position',state.vpos)



%EOF