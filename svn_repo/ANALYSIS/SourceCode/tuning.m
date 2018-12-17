function varargout = tuning(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tuning_OpeningFcn, ...
                   'gui_OutputFcn',  @tuning_OutputFcn, ...
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


% --- Executes just before tuning is made visible.
function tuning_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for tuning
handles.output = hObject;

% TAKE OVER PROVIDED PARAMETER SETTINGS
handles.vars.tune = varargin{1};
handles.vars.tune.reanalyze = 0;

% SET HANDLES
setall(hObject, handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tuning wait for user response (see UIRESUME)
uiwait(handles.tuning);


% --- Outputs from this function are returned to the command line.
function varargout = tuning_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.vars.tune;

% Hint: delete(hObject) closes the figure
delete(hObject);

function setall(hObject, handles)
% SET ALL HANDLES IN GUI TO PRESENT VALUES
set(handles.OpPointValue, 'String', num2str(handles.vars.tune.lunging.thresh,'%f4.2'));
set(handles.OpPointValue, 'Value', handles.vars.tune.lunging.thresh);
set(handles.OperatingPoint, 'Value', handles.vars.tune.lunging.thresh);
set(handles.checkbox_correct_orient, 'Value', handles.vars.tune.correct_orient);
set(handles.checkbox_correct_positions, 'Value', handles.vars.tune.correct_positions);
guidata(hObject, handles);

function DefaultButton_Callback(hObject, eventdata, handles)
% SET DEFAULT DETECTION PARAMETERS
handles.vars.tune = default_parameters;
setall(hObject, handles);
guidata(hObject, handles);

function LoadButton_Callback(hObject, eventdata, handles)
% LOAD DIFFERENT PARAMETER SET
[FileName,PathName] = uigetfile({'*.par','Parameter Files (*.par)';'*.*','All Files (*.*)'},...
                                'Select Parameter File');
if FileName,
    load('-mat',[PathName FileName]);
    handles.vars.tune = default;
    setall(hObject, handles);
end
guidata(hObject, handles);

function SaveButton_Callback(hObject, eventdata, handles)
% SAVE PRESENT PARAMETER SET
[FileName,PathName] = uiputfile({'*.par','Parameter Files (*.par)';'*.*','All Files (*.*)'},...
                                     'Save Parameter File');
if FileName,
    default = handles.vars.tune; %#ok<NASGU>
    save([PathName FileName],'default');
end
guidata(hObject, handles);

function OperatingPoint_Callback(hObject, eventdata, handles)
% SLIDER OPERATING POINT
old = handles.vars.tune.lunging.thresh;
handles.vars.tune.lunging.thresh = get(hObject,'Value');
if old ~= handles.vars.tune.lunging.thresh,
    handles.vars.tune.reanalyze = 1;
end
setall(hObject, handles);
guidata(hObject, handles);

function OperatingPoint_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function OpPointValue_Callback(hObject, eventdata, handles)
% VALUE OPERATING POINT
handles.vars.tune.lunging.thresh = get(hObject,'Value');
setall(hObject, handles);
guidata(hObject, handles);

function checkbox_correct_orient_Callback(hObject, eventdata, handles)
% CORRECT ORIENTATIONS (most likely path)
old = handles.vars.tune.correct_orient;
handles.vars.tune.correct_orient = get(hObject,'Value');
if old ~= handles.vars.tune.correct_orient,
    handles.vars.tune.reanalyze = 1;
end
setall(hObject, handles);
guidata(hObject, handles);

function checkbox_correct_positions_Callback(hObject, eventdata, handles)
% CORRECT FLY POSITIONS (most likely path)
old = handles.vars.tune.correct_positions;
handles.vars.tune.correct_positions = get(hObject,'Value');
if old ~= handles.vars.tune.correct_positions,
    handles.vars.tune.reanalyze = 1;
end
setall(hObject, handles);
guidata(hObject, handles);


function OpPointValue_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function slider1_Callback(hObject, eventdata, handles)

function slider1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function tuning_CloseRequestFcn(hObject, eventdata, handles)
% 'CLOSE' BUTTON (red 'X') HANDLING - DO NOT CHANGE.
uiresume(handles.tuning);
guidata(hObject, handles);


function edit2_Callback(hObject, eventdata, handles)

function edit2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
