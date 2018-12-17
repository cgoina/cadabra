function varargout = border(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @border_OpeningFcn, ...
                   'gui_OutputFcn',  @border_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

function border_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for border
global Files slashstr;
handles.output = hObject;

border_width = 4;
path=Files.strInVideoPath;
file = Files.strInVideoFName;
%path = varargin{1}.file(1).PathName;
%file = varargin{1}.file(1).FileName;
    
% Load background image
BackgroundFileName = [path slashstr file '_meanimg.dat_win'];
[BackgroundFID] = fopen(BackgroundFileName,'r');
nrows = fread(BackgroundFID,1,'double');
ncols = fread(BackgroundFID,1,'double');
handles.vars.mean_image = fread(BackgroundFID,[nrows ncols],'double');
fclose(BackgroundFID);

% Load ROI
numchambers = length(dir([path file '*_roi.mat']));

%JL092809 Changed from round(numChambers/2)) to 1 because it caused a
%problem when using GView to acquire movies.
%load([path file '_' num2str(round(numchambers/2)) '_roi.mat']);
load([path slashstr file '_1_roi.mat']);
handles.vars.ROI = ROI;
handles.vars.scale = scale;
handles.vars.r0 = numel(ROI.cols)/2; % minus 4-pixel border from qtrak
%handles.vars.r = handles.vars.r0-...
%    varargin{1}.params.border_width/handles.vars.scale.x;

handles.vars.r = handles.vars.r0-...
    border_width/handles.vars.scale.x;

min_r = 2.5; % mm

d = 0.1/handles.vars.scale.x/(handles.vars.r0-min_r/handles.vars.scale.x);
set(handles.slider,'SliderStep',[d d]);
set(handles.slider,'Min',min_r/handles.vars.scale.x/handles.vars.r0);
set(handles.slider,'Value',handles.vars.r/handles.vars.r0);
set(handles.edit_radius,'Value',handles.vars.r*handles.vars.scale.x);
set(handles.edit_radius,'String',...
    num2str(handles.vars.r*handles.vars.scale.x,'%4.1f'));

% Display one arena
figure(handles.border); cla;
imagesc(handles.vars.mean_image(handles.vars.ROI.rows,handles.vars.ROI.cols));
axis equal tight off;
set(gca,'XTick',[],'YTick',[]);
colormap(gray); hold on; drawnow;

plot(handles.vars.r0, handles.vars.r0, 'r+');
rectangle('Position',[handles.vars.r0-handles.vars.r, ...
    handles.vars.r0-handles.vars.r, ...
    2*handles.vars.r, 2*handles.vars.r],...
    'Curvature', [1,1], 'edgecolor', 'b', 'linewidth', 1.5);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes border wait for user response (see UIRESUME)
uiwait(handles.border);


function varargout = border_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.vars.r * handles.vars.scale.x;
varargout{2} = handles.vars.r0 * handles.vars.scale.x;
% Hint: delete(hObject) closes the figure
delete(hObject);


function image_CreateFcn(hObject, eventdata, handles)


function edit_radius_Callback(hObject, eventdata, handles)
handles.vars.r = str2double(get(handles.edit_radius,'String'))/handles.vars.scale.x;
if handles.vars.r*handles.vars.scale.x < 1,
    handles.vars.r = 1/handles.vars.scale.x;
end
set(handles.edit_radius,'Value',handles.vars.r*handles.vars.scale.x);
set(handles.edit_radius,'String',...
    num2str(handles.vars.r*handles.vars.scale.x,'%4.1f'));
plot_radius(hObject, handles)

function slider_Callback(hObject, eventdata, handles)
handles.vars.r = get(handles.slider,'Value')*handles.vars.r0;
set(handles.edit_radius,'Value',handles.vars.r*handles.vars.scale.x);
set(handles.edit_radius,'String',...
    num2str(handles.vars.r*handles.vars.scale.x,'%4.1f'));
plot_radius(hObject, handles)

function slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function edit_radius_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function border_CloseRequestFcn(hObject, eventdata, handles)
button = questdlg('Correct?','select analysis area','Yes','No','Yes');
if strcmp(button,'Yes') || ~numel(button),
    uiresume(handles.border);
    guidata(hObject, handles);
end

function plot_radius(hObject, handles)
figure(handles.border); cla;
imagesc(handles.vars.mean_image(handles.vars.ROI.rows,handles.vars.ROI.cols));
axis equal tight off;
set(gca,'XTick',[],'YTick',[]);
colormap(gray); hold on; drawnow;

plot(handles.vars.r0, handles.vars.r0, 'r+');
rectangle('Position',[handles.vars.r0-handles.vars.r, ...
    handles.vars.r0-handles.vars.r, ...
    2*handles.vars.r, 2*handles.vars.r],...
    'Curvature', [1,1], 'edgecolor', 'b', 'linewidth', 1.5);

guidata(hObject, handles);
