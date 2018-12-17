%% Analysis
% Copyright (C) 2008 Heiko Dankert, California Institute of Technology

% This file is part of ANALYSIS
% and the "Caltech Automated Drosophila Aggression-Courtship 
% Behavioral Repertoire Analysis (CADABRA)".

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%
% * Implementation by Heiko Dankert
%
%% Analysis GUI (MAIN PROGRAM)

function varargout = analysis(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @analysis_OpeningFcn, ...
                   'gui_OutputFcn',  @analysis_OutputFcn, ...
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


function analysis_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

% WIDTH OF BORDER THAT MAY BE EXCLUDED FROM ANALYSIS
handles.vars.params.border_width0 = 4; % mm
handles.vars.chamber.width0 = 40;
handles.vars.chamber.height0 = 50;
handles.vars.params.radius = -99; % default (no circular arena)
handles.vars.params.tune.correct_orient = 0; % correct fly orientation
handles.vars.params.tune.correct_positions = 0; % find most likely fly path
handles.vars.params.tune.reanalyze = 0;

% START VALUES FOR GUI and PLOT RANGES (maximum y-axis)
handles.reset = 0;
handles.vars.file_index = 0;
handles.vars.SavePathName = [];
% handles.vars.new = 1;
handles.vars.auto = 1;
handles.vars.params.fid = 0;
handles.vars.editgen = 0;
handles.vars.editchamber = 0;
handles.vars.oldgen = '';
if ispc, handles.vars.params.slash = '\'; else handles.vars.params.slash = '/'; end

handles.vars.file = []; 
% handles.vars.file.Chambers = 1;
% handles.vars.file.Genotypes = 1;
handles.vars.chamber.width = handles.vars.chamber.width0;
handles.vars.chamber.height = handles.vars.chamber.height0;
handles.vars.params.border_width = handles.vars.params.border_width0;
set(handles.edit_chamberwidth,'String',num2str(handles.vars.chamber.width,'%3.0f'));
set(handles.edit_chamberheight,'String',num2str(handles.vars.chamber.height,'%3.0f'));
set(handles.edit_border,'String',num2str(handles.vars.params.border_width,'%4.1f'));
handles.vars.chamber_index = []; handles.vars.file_list = [];
handles.vars.params.tim = get(handles.cumul_time,'Value');
handles.vars.params.range.heatmaps.occ.chase = 10;
handles.vars.params.range.heatmaps.occ.pos = 100;
handles.vars.params.range.heatmaps.occ.vel = 10;
handles.vars.params.range.heatmaps.occ.wing = 10;
handles.vars.params.range.heatmaps.occ.circl = 50;
handles.vars.params.range.heatmaps.occ.lunge = 10;
handles.vars.params.range.heatmaps.occ.wingthreat = 10;
handles.vars.params.range.heatmaps.occ.tussl = 10;
handles.vars.params.range.stat.occ.walk = 10;
handles.vars.params.range.stat.occ.flydist = 20;
handles.vars.params.range.stat.occ.vel = 20;
handles.vars.params.range.stat.occ.wing = 100;
handles.vars.params.range.stat.occ.chase = 100;
handles.vars.params.range.stat.occ.circl = 100;
handles.vars.params.range.stat.occ.lunge = 50;
handles.vars.params.range.stat.occ.wingthreat = 10;
handles.vars.params.range.stat.occ.tussl = 50;
handles.vars.params.range.heatmaps.tim.chase = 1;
handles.vars.params.range.heatmaps.tim.pos = 10;
handles.vars.params.range.heatmaps.tim.vel = 10;
handles.vars.params.range.heatmaps.tim.wing = 1;
handles.vars.params.range.heatmaps.tim.circl = 5;
handles.vars.params.range.heatmaps.tim.lunge = 1;
handles.vars.params.range.heatmaps.tim.wingthreat = 1;
handles.vars.params.range.heatmaps.tim.tussl = 1;
handles.vars.params.range.stat.tim.walk = 10;
handles.vars.params.range.stat.tim.flydist = 20;
handles.vars.params.range.stat.tim.vel = 20;
handles.vars.params.range.stat.tim.wing = 1;
handles.vars.params.range.stat.tim.chase = 20;
handles.vars.params.range.stat.tim.circl = 10;
handles.vars.params.range.stat.tim.lunge = 1;
handles.vars.params.range.stat.tim.wingthreat = 1;
handles.vars.params.range.stat.tim.tussl = 10;
handles = set_ranges(hObject, handles);

handles.vars.params.plots.heatmaps.position_vec = 0;
handles.vars.params.plots.heatmaps.velocity_vec = 0;
handles.vars.params.plots.heatmaps.wing_vec = 0;
handles.vars.params.plots.heatmaps.wingrel_vec = 0;
handles.vars.params.plots.heatmaps.chase_vec = 0;
handles.vars.params.plots.heatmaps.circl_vec = 0;
handles.vars.params.plots.heatmaps.lunge_vec = 0;
handles.vars.params.plots.heatmaps.wingthreat_vec = 0;
handles.vars.params.plots.heatmaps.tussl_vec = 0;
handles.vars.params.plots.stat.flydist_vec = 0;
handles.vars.params.plots.stat.body_vec = 0;
handles.vars.params.plots.stat.vel_vec = 0;
handles.vars.params.plots.stat.wing_vec = 0;
handles.vars.params.plots.stat.lunge_vec = 0;
handles.vars.params.plots.stat.wingthreat_vec = 0;
handles.vars.params.plots.stat.tussl_vec = 0;
handles.vars.params.plots.stat.scat_vec = 0;
handles.vars.params.plots.stat.copul_vec = 0;
handles.vars.params.plots.timeseries_vec = 0;
handles.vars.params.plots.stat.actionrel_vec = 0;
handles.vars.params.plots.ethogram_vec = 0;

% Set default detection parameters
handles.vars.params.tune = default_parameters;

guidata(hObject, handles);

function varargout = analysis_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;


function handles = get_fileindex(handles,hObject)
% GET THE ACTUAL FILE INDEX
for i=1:length(handles.vars.file),
    for ii=1:length(handles.vars.file(i).Genotypes)
        if strcmp(handles.vars.file(i).Genotypes{ii},handles.vars.oldgen),
            handles.vars.file_index = i;
            guidata(hObject, handles);  % Update handles structure
            return;
        end
    end
end
    
function handles = addfiles(FileName,handles,hObject)
% ADD FILES TO PRESENT STRUCTURE
if isstruct(FileName),
    FileName = struct2cell(FileName);
    FileName = FileName(1,:);
end
if isempty(handles.vars.SavePathName),
    % Set the config file path and default file name, in case it has 
    % no been set yet
    handles.vars.SavePathName = handles.vars.PathName;
    handles.vars.saveas = 'analysis.cfg';
    set(handles.edit_save,'String',[handles.vars.SavePathName handles.vars.saveas]);
end
% if handles.vars.new,
%     set(handles.edit_save, 'String', [handles.vars.PathName handles.vars.saveas]);
%     handles.vars.new = 0;
% end

% New file name or path?
if ~(isequal(FileName, 0) || isequal(handles.vars.PathName, 0)),
    % File list or single file?
    if iscell(FileName),
        [sorted_names,sorted_index] = sortrows(FileName');
    else
        sorted_names = cell(1,1);
        sorted_names{1} = FileName;
    end
    % Extract file names without arena number
    fnames = cell(1,length(sorted_names));
    for j=1:length(sorted_names),
        fn = sorted_names{j}; in = strfind(fn,'_');
        fnames{j} = fn(1:in(end)-1);
    end
    [b,m,n] = unique(fnames); nfiles = numel(b);

    % Scan for all arenas per file
    chamb = cell(1,nfiles);
    for j=1:nfiles,
        FileName = dir([handles.vars.PathName b{j} '_*.feat']);
        chamb{j} = zeros(1,numel(FileName));
        for jj=1:numel(FileName),
            fn = FileName(jj).name; in = strfind(fn,'_'); in = in(end);
            chamb{j}(jj) = str2num(fn(in+1:end-5));
        end
        [chamb{j},si] = sort(chamb{j}); fname = FileName;
        for jj=1:numel(FileName),
            FileName(jj).name = fname(si(jj)).name;
        end
    end

    % Update variables and lists in the GUI
    for j=1:nfiles,
        handles.vars.file_index = length(handles.vars.file) + 1;
        handles.vars.file(handles.vars.file_index).FileName = b{j}; 
        handles.vars.file(handles.vars.file_index).PathName = handles.vars.PathName;
        handles.vars.file(handles.vars.file_index).Chambers = chamb{j};
        for jj=1:numel(chamb{j}),
            handles.vars.file(handles.vars.file_index).Genotypes{jj} = b{j};
        end
        handles.vars.chamber_index(handles.vars.file_index) = 1;
        in = strfind(handles.vars.PathName,handles.vars.params.slash); in = in(end-1:end);
        handles.vars.file_list{handles.vars.file_index} = [handles.vars.PathName(in(1)+1:in(2)) b{j}];
        set(handles.listbox_files, 'String', handles.vars.file_list);
        set(handles.listbox_chambers, 'String', handles.vars.file(handles.vars.file_index).Chambers);
        set(handles.edit_genotype, 'String', handles.vars.file(handles.vars.file_index).Genotypes{handles.vars.chamber_index(handles.vars.file_index)});
        guidata(hObject, handles);  % Update handles structure
    end
    handles = refresh_genlist(handles,hObject,1);
    guidata(hObject, handles);  % Update handles structure
end

function handles = savedat(hObject, handles)
% AUTOMATICALLY SAVE GUI STATUS TO CONFIG FILE
handles = setall(hObject, handles);
guidata(hObject, handles);  % Update handles structure

if ~handles.reset,
    in = strfind(handles.vars.saveas,'.'); in = in(end)-1;
    savevars = handles.vars;
    save('-mat',[handles.vars.SavePathName handles.vars.saveas(1:in) '.cfg'],'savevars');
end


function handles = refresh_genlist(handles,hObject,renew)
% REFRESH GENOTYPE LIST
if renew,
    % Complete actualization of genotype list
    fnames = []; cnt = 0; 
    for j=1:numel(handles.vars.file),
        for jj=1:numel(handles.vars.file(j).Genotypes),
            cnt = cnt + 1;
            fnames{cnt} = handles.vars.file(j).Genotypes{jj};
        end
    end
    % Sort genotype list
    [a,b] = unique(fnames); handles.vars.gens = fnames(sort(b));
    ngens = numel(handles.vars.gens); handles.vars.n = zeros(1,ngens);
    % Update n listbox entries
    for j=1:ngens,
        handles.vars.n(j) = numel(find(strcmp(fnames,handles.vars.gens(j))));
    end
    % Currently selected genotype?
    if handles.vars.file_index && numel(handles.vars.gens), 
        in = find(strcmp(handles.vars.gens,handles.vars.file(handles.vars.file_index).Genotypes{handles.vars.chamber_index(handles.vars.file_index)})); 
    end
    % Refresh number of rows and columns in GUI and for the
    % plot routines (currently all genotypes onto one sheet of paper)
    if ngens,
        handles.vars.params.k(1) = floor(sqrt(ngens));
        handles.vars.params.k(2) = ceil(ngens/handles.vars.params.k(1));
        set(handles.edit_rows,'String',num2str(handles.vars.params.k(1)));
        set(handles.edit_cols,'String',num2str(handles.vars.params.k(2)));
    end
else
    % Refresh label index of genotype and n listbox in GUI with the currently selected genotype
    ngens = numel(handles.vars.gens);
    in = get(handles.listbox_genotypes,'Value');
    oldg = handles.vars.gens; oldn = handles.vars.n;
    handles.vars.gens(handles.vars.slider_val) = oldg(in); 
    handles.vars.gens(in) = oldg(handles.vars.slider_val);
    handles.vars.n(handles.vars.slider_val) = oldn(in); 
    handles.vars.n(in) = oldn(handles.vars.slider_val);
    % Currently selected genotype?
    in = handles.vars.slider_val;
end
guidata(hObject, handles);

if handles.vars.file_index,
    % Update genotype and n listbox
    if numel(handles.vars.gens),
        set(handles.listbox_genotypes, 'String', handles.vars.gens, 'Value', in);
        set(handles.listbox_n, 'String', handles.vars.n, 'Value', in);
        if ngens>1, set(handles.slider, 'SliderStep', [1/(ngens-1) 1/(ngens-1)], 'Value', (ngens-in)/(ngens-1)); end
    else
        set(handles.listbox_genotypes, 'String', '', 'Value', 0);
        set(handles.listbox_n, 'String', '', 'Value', 0);
    end
    
    % REFRESH ALL VARIABLES by reading out current GUI parameters and selections
    handles.vars.params.PSFileN = [handles.vars.SavePathName handles.vars.saveas(1:end-3) 'ps'];
    handles.vars.params.max_frames = str2double(get(handles.edit_timelimit,'String')) * 60; %time limit [s]
    handles.vars.params.cross_to = '';
    handles.vars.params.oneobj = get(handles.checkbox_oneobj,'Value');
    handles.vars.params.autoscale = get(handles.checkbox_autoscale,'Value');
    handles.vars.params.winlos = get(handles.checkbox_winlos,'Value');
    handles.vars.params.courtship = get(handles.checkbox_courtship,'Value');
    handles.vars.params.k = [str2double(get(handles.edit_rows,'String')) str2double(get(handles.edit_cols,'String'))];
    handles.vars.params.pdf = get(handles.checkbox_pdf,'Value');
    handles.vars.params.bool_xls = get(handles.checkbox_text,'Value');
    handles.vars.params.boxplot = get(handles.checkbox_boxplot,'Value');
    handles.vars.params.bonf = get(handles.checkbox_bonf,'Value');
%     handles.vars.params.feat_read_new = get(handles.checkbox_readnew,'Value');
    handles.vars.params.analyze_new = get(handles.checkbox_analyzenew,'Value');
    handles.vars.params.feat_read_new = handles.vars.params.analyze_new;
    handles.vars.params.border = get(handles.checkbox_border,'Value');
    handles.vars.params.axisfontsize = str2double(get(handles.edit_axisfontsize,'String'));
    handles.vars.chamber.width = str2double(get(handles.edit_chamberwidth,'String'));
    handles.vars.chamber.height = str2double(get(handles.edit_chamberheight,'String'));
    handles.vars.params.border_width = str2double(get(handles.edit_border,'String'));
    handles.vars.params.oneobj = get(handles.checkbox_oneobj,'Value');
    handles.vars.params.circular = get(handles.checkbox_circular,'Value');
    
    
    handles.vars.params.tim = get(handles.cumul_time,'Value');

    handles.vars.params.plots.heatmaps.position = get(handles.heatmap_position,'Value');
    handles.vars.params.plots.heatmaps.velocity = get(handles.heatmap_velocity,'Value');
    handles.vars.params.plots.heatmaps.ext_leftwing = get(handles.heatmap_ext_leftwing,'Value');
    handles.vars.params.plots.heatmaps.ext_rightwing = get(handles.heatmap_ext_rightwing,'Value');
    handles.vars.params.plots.heatmaps.ext_diffleftrightwing = get(handles.heatmap_ext_diffleftrightwing,'Value');
    handles.vars.params.plots.heatmaps.ext_onewing = get(handles.heatmap_ext_onewing,'Value');
    handles.vars.params.plots.heatmaps.chasing = get(handles.heatmap_chasing,'Value');
    handles.vars.params.plots.heatmaps.circling = get(handles.heatmap_circling,'Value');
    handles.vars.params.plots.heatmaps.tussling = get(handles.heatmap_tussling,'Value');
    handles.vars.params.plots.heatmaps.lunging = get(handles.heatmap_lunging,'Value');
    handles.vars.params.plots.heatmaps.wingthreat = get(handles.heatmap_wingthreat,'Value');
    handles.vars.params.plots.heatmaps.wingrel = get(handles.heatmap_wingrel,'Value');

    handles.vars.params.plots.stat.velocity = get(handles.stat_velocity,'Value');
    handles.vars.params.plots.stat.chasing = get(handles.stat_chasing,'Value');
    handles.vars.params.plots.stat.traveldist = get(handles.stat_walkingdist,'Value');
    handles.vars.params.plots.stat.dist = get(handles.stat_dist,'Value');
    handles.vars.params.plots.stat.movestop = get(handles.stat_movestop,'Value');
    handles.vars.params.plots.stat.body = get(handles.stat_bodysize,'Value');
    handles.vars.params.plots.stat.dendro = get(handles.stat_dendro,'Value');
    handles.vars.params.plots.stat.actionrel = get(handles.stat_actionrel,'Value');
    handles.vars.params.plots.stat.tussl = get(handles.stat_tussl,'Value');
    handles.vars.params.plots.stat.lunging = get(handles.stat_lunging,'Value');
    handles.vars.params.plots.stat.wingthreat = get(handles.stat_wingthreat,'Value');
    handles.vars.params.plots.stat.scat = get(handles.stat_scat,'Value');
    handles.vars.params.plots.stat.leftwing = get(handles.stat_leftwing,'Value');
    handles.vars.params.plots.stat.rightwing = get(handles.stat_rightwing,'Value');
    handles.vars.params.plots.stat.onewing = get(handles.stat_onewing,'Value');
    handles.vars.params.plots.stat.circling = get(handles.stat_circling,'Value');
    handles.vars.params.plots.stat.wingextflydist = get(handles.stat_wingextflydist,'Value');
    handles.vars.params.plots.stat.copulation = get(handles.stat_copulation,'Value');

    handles.vars.params.plots.timeseries = get(handles.timeseries,'Value');
    handles.vars.params.plots.ethogram = get(handles.ethogram,'Value');
    handles.vars.params.plots.movieclips = get(handles.movieclips,'Value');
    handles.vars.params.plots.showhtml = get(handles.showhtml,'Value');

    if handles.vars.params.tim,
        handles.vars.params.range.heatmaps.tim.chase = (abs(str2double(get(handles.edit_chaseheat,'String'))));
        handles.vars.params.range.heatmaps.tim.pos = (abs(str2double(get(handles.edit_posheat,'String'))));
        handles.vars.params.range.heatmaps.tim.vel = (abs(str2double(get(handles.edit_velheat,'String'))));
        handles.vars.params.range.heatmaps.tim.wing = (abs(str2double(get(handles.edit_wingheat,'String'))));
        handles.vars.params.range.heatmaps.tim.circl = (abs(str2double(get(handles.edit_circlheat,'String'))));
        handles.vars.params.range.heatmaps.tim.lunge = (abs(str2double(get(handles.edit_lungeheat,'String'))));
        handles.vars.params.range.heatmaps.tim.wingthreat = (abs(str2double(get(handles.edit_wingthreatheat,'String'))));
        handles.vars.params.range.heatmaps.tim.tussl = (abs(str2double(get(handles.edit_tusslheat,'String'))));
        handles.vars.params.range.stat.tim.walk = (abs(str2double(get(handles.edit_walkstat,'String'))));
        handles.vars.params.range.stat.tim.flydist = (abs(str2double(get(handles.edit_flydiststat,'String'))));
        handles.vars.params.range.stat.tim.vel = (abs(str2double(get(handles.edit_velstat,'String'))));
        handles.vars.params.range.stat.tim.wing = (abs(str2double(get(handles.edit_wingstat,'String'))));
        handles.vars.params.range.stat.tim.chase = (abs(str2double(get(handles.edit_chasestat,'String'))));
        handles.vars.params.range.stat.tim.circl = (abs(str2double(get(handles.edit_circlstat,'String'))));
        handles.vars.params.range.stat.tim.lunge = (abs(str2double(get(handles.edit_lungestat,'String'))));
        handles.vars.params.range.stat.tim.tussl = (abs(str2double(get(handles.edit_tusslstat,'String'))));
        handles.vars.params.range.stat.tim.wingthreat = (abs(str2double(get(handles.edit_wingthreatstat,'String'))));
    else
        handles.vars.params.range.heatmaps.occ.chase = (abs(str2double(get(handles.edit_chaseheat,'String'))));
        handles.vars.params.range.heatmaps.occ.pos = (abs(str2double(get(handles.edit_posheat,'String'))));
        handles.vars.params.range.heatmaps.occ.vel = (abs(str2double(get(handles.edit_velheat,'String'))));
        handles.vars.params.range.heatmaps.occ.wing = (abs(str2double(get(handles.edit_wingheat,'String'))));
        handles.vars.params.range.heatmaps.occ.circl = (abs(str2double(get(handles.edit_circlheat,'String'))));
        handles.vars.params.range.heatmaps.occ.lunge = (abs(str2double(get(handles.edit_lungeheat,'String'))));
        handles.vars.params.range.heatmaps.occ.wingthreat = (abs(str2double(get(handles.edit_wingthreatheat,'String'))));
        handles.vars.params.range.heatmaps.occ.tussl = (abs(str2double(get(handles.edit_tusslheat,'String'))));
        handles.vars.params.range.stat.occ.walk = (abs(str2double(get(handles.edit_walkstat,'String'))));
        handles.vars.params.range.stat.occ.flydist = (abs(str2double(get(handles.edit_flydiststat,'String'))));
        handles.vars.params.range.stat.occ.vel = (abs(str2double(get(handles.edit_velstat,'String'))));
        handles.vars.params.range.stat.occ.wing = (abs(str2double(get(handles.edit_wingstat,'String'))));
        handles.vars.params.range.stat.occ.chase = (abs(str2double(get(handles.edit_chasestat,'String'))));
        handles.vars.params.range.stat.occ.circl = (abs(str2double(get(handles.edit_circlstat,'String'))));
        handles.vars.params.range.stat.occ.lunge = (abs(str2double(get(handles.edit_lungestat,'String'))));
        handles.vars.params.range.stat.occ.wingthreat = (abs(str2double(get(handles.edit_wingthreatstat,'String'))));
        handles.vars.params.range.stat.occ.tussl = (abs(str2double(get(handles.edit_tusslstat,'String'))));
    end
else
    % If filelist is empty deactivate everything
    enable_all('off', hObject, handles);
end
handles = savedat(hObject, handles);


function handles = setall(hObject, handles)
% REFRESH GUI WITH CURRENT VALUES AND ENTRIES OF VARIABLES
ngens = numel(handles.vars.gens);
if prod(handles.vars.params.k) > ngens && numel(handles.vars.gens),
    handles.vars.params.k(1) = floor(sqrt(ngens));
    handles.vars.params.k(2) = ceil(ngens/handles.vars.params.k(1));
end
set(handles.listbox_files, 'String', handles.vars.file_list,'Value',handles.vars.file_index);
if handles.vars.file_index && numel(handles.vars.gens),
    set(handles.listbox_chambers, 'String', handles.vars.file(handles.vars.file_index).Chambers);
    set(handles.edit_genotype, 'String', handles.vars.file(handles.vars.file_index).Genotypes{handles.vars.chamber_index(handles.vars.file_index)});
else
    set(handles.listbox_chambers, 'String', '');
    set(handles.edit_genotype, 'String', '');
end
% set(handles.edit_save,'String',[handles.vars.SavePathName handles.vars.saveas]);
handles.vars.params.PSFileN = [handles.vars.SavePathName handles.vars.saveas(1:end-3) 'ps'];
set(handles.edit_timelimit,'String',num2str(round(handles.vars.params.max_frames/60)));
handles.vars.params.cross_to = '';
set(handles.checkbox_oneobj,'Value',handles.vars.params.oneobj);
set(handles.checkbox_autoscale,'Value',handles.vars.params.autoscale);
set(handles.checkbox_winlos,'Value',handles.vars.params.winlos);
set(handles.checkbox_courtship,'Value',handles.vars.params.courtship);
set(handles.edit_rows,'String',num2str(handles.vars.params.k(1)));
set(handles.edit_cols,'String',num2str(handles.vars.params.k(2)));
set(handles.checkbox_pdf,'Value',handles.vars.params.pdf);
set(handles.checkbox_text,'Value',handles.vars.params.bool_xls);
set(handles.checkbox_boxplot,'Value',handles.vars.params.boxplot);
set(handles.checkbox_bonf,'Value',handles.vars.params.bonf);
set(handles.checkbox_meansem,'Value',1-handles.vars.params.boxplot);
% set(handles.checkbox_readnew,'Value',handles.vars.params.feat_read_new);
set(handles.checkbox_analyzenew,'Value',handles.vars.params.analyze_new);
set(handles.checkbox_border,'Value',handles.vars.params.border);
set(handles.edit_axisfontsize,'String',num2str(handles.vars.params.axisfontsize));
set(handles.edit_chamberwidth,'String',num2str(handles.vars.chamber.width,'%3.0f'));
set(handles.edit_chamberheight,'String',num2str(handles.vars.chamber.height,'%3.0f'));
set(handles.edit_border,'String',num2str(handles.vars.params.border_width,'%4.1f'));
set(handles.checkbox_circular,'Value',handles.vars.params.circular);
set(handles.cumul_time,'Value',handles.vars.params.tim);

set(handles.heatmap_position,'Value',handles.vars.params.plots.heatmaps.position);
set(handles.heatmap_velocity,'Value',handles.vars.params.plots.heatmaps.velocity);
set(handles.heatmap_ext_leftwing,'Value',handles.vars.params.plots.heatmaps.ext_leftwing);
set(handles.heatmap_ext_rightwing,'Value',handles.vars.params.plots.heatmaps.ext_rightwing);
set(handles.heatmap_ext_diffleftrightwing,'Value',handles.vars.params.plots.heatmaps.ext_diffleftrightwing);
set(handles.heatmap_ext_onewing,'Value',handles.vars.params.plots.heatmaps.ext_onewing);
set(handles.heatmap_chasing,'Value',handles.vars.params.plots.heatmaps.chasing);
set(handles.heatmap_circling,'Value',handles.vars.params.plots.heatmaps.circling);
set(handles.heatmap_tussling,'Value',handles.vars.params.plots.heatmaps.tussling);
set(handles.heatmap_lunging,'Value',handles.vars.params.plots.heatmaps.lunging);
set(handles.heatmap_wingthreat,'Value',handles.vars.params.plots.heatmaps.wingthreat);
set(handles.heatmap_wingrel,'Value',handles.vars.params.plots.heatmaps.wingrel);

set(handles.stat_velocity,'Value',handles.vars.params.plots.stat.velocity);
set(handles.stat_chasing,'Value',handles.vars.params.plots.stat.chasing);
set(handles.stat_walkingdist,'Value',handles.vars.params.plots.stat.traveldist);
set(handles.stat_dist,'Value',handles.vars.params.plots.stat.dist);
set(handles.stat_movestop,'Value',handles.vars.params.plots.stat.movestop);
set(handles.stat_bodysize,'Value',handles.vars.params.plots.stat.body);
set(handles.stat_dendro,'Value',handles.vars.params.plots.stat.dendro);
set(handles.stat_actionrel,'Value',handles.vars.params.plots.stat.actionrel);
set(handles.stat_tussl,'Value',handles.vars.params.plots.stat.tussl);
set(handles.stat_lunging,'Value',handles.vars.params.plots.stat.lunging);
set(handles.stat_wingthreat,'Value',handles.vars.params.plots.stat.wingthreat);
set(handles.stat_scat,'Value',handles.vars.params.plots.stat.scat);
set(handles.stat_leftwing,'Value',handles.vars.params.plots.stat.leftwing);
set(handles.stat_rightwing,'Value',handles.vars.params.plots.stat.rightwing);
set(handles.stat_onewing,'Value',handles.vars.params.plots.stat.onewing);
set(handles.stat_circling,'Value',handles.vars.params.plots.stat.circling);
set(handles.stat_wingextflydist,'Value',handles.vars.params.plots.stat.wingextflydist);
set(handles.stat_copulation,'Value',handles.vars.params.plots.stat.copulation);

set(handles.timeseries,'Value',handles.vars.params.plots.timeseries);
set(handles.ethogram,'Value',handles.vars.params.plots.ethogram);
set(handles.movieclips,'Value',handles.vars.params.plots.movieclips);
set(handles.showhtml,'Value',handles.vars.params.plots.showhtml);

handles = set_ranges(hObject, handles);


function handles = set_ranges(hObject, handles)
if handles.vars.params.tim,
    set(handles.edit_chaseheat,'String',num2str((handles.vars.params.range.heatmaps.tim.chase)));
    set(handles.edit_posheat,'String',num2str((handles.vars.params.range.heatmaps.tim.pos)));
    set(handles.edit_velheat,'String',num2str((handles.vars.params.range.heatmaps.tim.vel)));
    set(handles.edit_wingheat,'String',num2str((handles.vars.params.range.heatmaps.tim.wing)));
    set(handles.edit_circlheat,'String',num2str((handles.vars.params.range.heatmaps.tim.circl)));
    set(handles.edit_lungeheat,'String',num2str((handles.vars.params.range.heatmaps.tim.lunge)));
    set(handles.edit_wingthreatheat,'String',num2str((handles.vars.params.range.heatmaps.tim.wingthreat)));
    set(handles.edit_tusslheat,'String',num2str((handles.vars.params.range.heatmaps.tim.tussl)));
    set(handles.edit_walkstat,'String',num2str((handles.vars.params.range.stat.tim.walk)));
    set(handles.edit_flydiststat,'String',num2str((handles.vars.params.range.stat.tim.flydist)));
    set(handles.edit_velstat,'String',num2str((handles.vars.params.range.stat.tim.vel)));
    set(handles.edit_wingstat,'String',num2str((handles.vars.params.range.stat.tim.wing)));
    set(handles.edit_chasestat,'String',num2str((handles.vars.params.range.stat.tim.chase)));
    set(handles.edit_circlstat,'String',num2str((handles.vars.params.range.stat.tim.circl)));
    set(handles.edit_lungestat,'String',num2str((handles.vars.params.range.stat.tim.lunge)));
    set(handles.edit_tusslstat,'String',num2str((handles.vars.params.range.stat.tim.tussl)));
    set(handles.edit_wingthreatstat,'String',num2str((handles.vars.params.range.stat.tim.wingthreat)));
else
    set(handles.edit_chaseheat,'String',num2str((handles.vars.params.range.heatmaps.occ.chase)));
    set(handles.edit_posheat,'String',num2str((handles.vars.params.range.heatmaps.occ.pos)));
    set(handles.edit_velheat,'String',num2str((handles.vars.params.range.heatmaps.occ.vel)));
    set(handles.edit_wingheat,'String',num2str((handles.vars.params.range.heatmaps.occ.wing)));
    set(handles.edit_circlheat,'String',num2str((handles.vars.params.range.heatmaps.occ.circl)));
    set(handles.edit_lungeheat,'String',num2str((handles.vars.params.range.heatmaps.occ.lunge)));
    set(handles.edit_wingthreatheat,'String',num2str((handles.vars.params.range.heatmaps.occ.wingthreat)));
    set(handles.edit_tusslheat,'String',num2str((handles.vars.params.range.heatmaps.occ.tussl)));
    set(handles.edit_walkstat,'String',num2str((handles.vars.params.range.stat.occ.walk)));
    set(handles.edit_flydiststat,'String',num2str((handles.vars.params.range.stat.occ.flydist)));
    set(handles.edit_velstat,'String',num2str((handles.vars.params.range.stat.occ.vel)));
    set(handles.edit_wingstat,'String',num2str((handles.vars.params.range.stat.occ.wing)));
    set(handles.edit_chasestat,'String',num2str((handles.vars.params.range.stat.occ.chase)));
    set(handles.edit_circlstat,'String',num2str((handles.vars.params.range.stat.occ.circl)));
    set(handles.edit_lungestat,'String',num2str((handles.vars.params.range.stat.occ.lunge)));
    set(handles.edit_wingthreatstat,'String',num2str((handles.vars.params.range.stat.occ.wingthreat)));
    set(handles.edit_tusslstat,'String',num2str((handles.vars.params.range.stat.occ.tussl)));
end


function load_Callback(hObject, eventdata, handles)
% SELECT AND LOAD CONFIG FILE
[FileName,PathName] = uigetfile({'*.cfg','Configuration Files (*.cfg)';'*.*','All Files (*.*)'},...
    'Load Configuration');
if FileName,
    load('-mat',[PathName FileName]);
    handles.vars = savevars;
    if handles.vars.file_index,
        handles = setall(hObject, handles);
        enable_all('on', hObject, handles);
        handles = refresh_genlist(handles,hObject,1);
    end
end
guidata(hObject, handles);

% --- Executes on selection change in listbox_files.
function listbox_files_Callback(hObject, eventdata, handles)
% SELECT A FILE IN FILE LISTBOX AND UPDATE 
% GENOTYPE LISTBOX & EDITBOX, N & CHAMBER & LISTBOX
if ~isempty(handles.vars.file),
    handles.vars.file_index = get(handles.listbox_files,'Value');
    handles.vars.editgen = 0;
    set(handles.listbox_chambers, 'String', handles.vars.file(handles.vars.file_index).Chambers, ...
        'Value', handles.vars.chamber_index(handles.vars.file_index));
    set(handles.edit_genotype, 'String', handles.vars.file(handles.vars.file_index).Genotypes{handles.vars.chamber_index(handles.vars.file_index)});
    in = find(strcmp(handles.vars.gens,handles.vars.file(handles.vars.file_index).Genotypes{handles.vars.chamber_index(handles.vars.file_index)}));
    set(handles.listbox_genotypes, 'String', handles.vars.gens, 'Value', in);
    set(handles.listbox_n, 'String', handles.vars.n, 'Value', in);
    guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function listbox_chambers_Callback(hObject, eventdata, handles)
% SELECT A CHAMBER (ARENA) AND UPDATE 
% GENOTYPE LISTBOX & EDITBOX, N LISTBOX
handles.vars.chamber_index(handles.vars.file_index) = get(handles.listbox_chambers,'Value');
handles.vars.editgen = 0; handles.vars.editchamber = 1;
set(handles.edit_genotype, 'String', handles.vars.file(handles.vars.file_index).Genotypes{handles.vars.chamber_index(handles.vars.file_index)});
in = find(strcmp(handles.vars.gens,handles.vars.file(handles.vars.file_index).Genotypes{handles.vars.chamber_index(handles.vars.file_index)}));
set(handles.listbox_genotypes, 'String', handles.vars.gens, 'Value', in);
set(handles.listbox_n, 'String', handles.vars.n, 'Value', in);
guidata(hObject, handles);

function edit_genotype_Callback(hObject, eventdata, handles)
% EDIT GENOTYPE NAMING IN EDITBOX
% Any genotypes present?
if ~handles.vars.editgen,
    handles.vars.file(handles.vars.file_index).Genotypes{handles.vars.chamber_index(handles.vars.file_index)} = get(handles.edit_genotype,'String');
    handles = refresh_genlist(handles,hObject,1);
else
    % Replace former by new name
    newgen = get(handles.edit_genotype,'String'); icnt = 0; in = 1;
    for i=1:length(handles.vars.gens),
        if strcmp(handles.vars.gens{i},handles.vars.oldgen),
            in = i;
            handles.vars.gens{i} = newgen;
        end
    end
    % Update for all files under this genotype
    for i=1:length(handles.vars.file),
        for ii=1:length(handles.vars.file(i).Genotypes)
            if strcmp(handles.vars.file(i).Genotypes{ii},handles.vars.oldgen),
                handles.vars.file(i).Genotypes{ii} = newgen;
                if ~icnt,
                    handles.vars.file_index = i;
                    icnt = 1;
                end
            end
        end
    end
    % Update genotype listbox
    set(handles.listbox_genotypes, 'String', handles.vars.gens, 'Value', in);
    handles = refresh_genlist(handles,hObject,1);
end
guidata(hObject, handles);

function listbox_n_Callback(hObject, eventdata, handles)
% SELECT A GENOTYPE IN N LISTBOX AND UPDATE SLIDER, 
% GENOTYPE LISTBOX & EDITBOX, CHAMBER & FILE LISTBOX
in = get(handles.listbox_n,'Value');
handles.vars.editgen = 1;
set(handles.listbox_genotypes, 'Value', in);
set(handles.listbox_n, 'Value', in);
set(handles.edit_genotype, 'String', handles.vars.gens{in});
handles.vars.oldgen = handles.vars.gens{in};
handles = get_fileindex(handles,hObject);
set(handles.listbox_files,'Value',handles.vars.file_index);
set(handles.listbox_chambers, 'String', handles.vars.file(handles.vars.file_index).Chambers, ...
    'Value', handles.vars.chamber_index(handles.vars.file_index));ind = find(strcmp(handles.vars.file(handles.vars.file_index).Genotypes,handles.vars.gens{in}));
set(handles.listbox_chambers,'Value',ind(1));
ngens = numel(handles.vars.gens);
if ngens > 1, set(handles.slider, 'SliderStep', [1/(ngens-1) 1/(ngens-1)], 'Value', (ngens-in)/(ngens-1)); end
guidata(hObject, handles);

function listbox_genotypes_Callback(hObject, eventdata, handles)
% SELECT A GENOTYPE IN GENOTYPE LISTBOX AND UPDATE 
% GENOTYPE EDITBOX, SLIDER, N LISTBOX, CHAMBER & FILE LISTBOX
in = get(handles.listbox_genotypes,'Value');
handles.vars.editgen = 1;
set(handles.listbox_genotypes, 'Value', in);
set(handles.listbox_n, 'Value', in);
set(handles.edit_genotype, 'String', handles.vars.gens{in});
handles.vars.oldgen = handles.vars.gens{in};
handles = get_fileindex(handles,hObject);
set(handles.listbox_files,'Value',handles.vars.file_index);
set(handles.listbox_chambers, 'String', handles.vars.file(handles.vars.file_index).Chambers, ...
    'Value', handles.vars.chamber_index(handles.vars.file_index));
ind = find(strcmp(handles.vars.file(handles.vars.file_index).Genotypes,handles.vars.gens{in}));
set(handles.listbox_chambers,'Value',ind(1));
ngens = numel(handles.vars.gens);
if ngens > 1, set(handles.slider, 'SliderStep', [1/(ngens-1) 1/(ngens-1)], 'Value', (ngens-in)/(ngens-1)); end
% handles.vars.file(handles.vars.file_index).Genotypes{handles.vars.chamber_index(handles.vars.file_index)} = handles.vars.gens{in};
% handles = refresh_genlist(handles,hObject,0);
guidata(hObject, handles);

function slider_Callback(hObject, eventdata, handles)
% SELECT GENOTYPE WITH THE SLIDER
handles.vars.slider_val = numel(handles.vars.gens)-round(get(handles.slider,'Value')*(numel(handles.vars.gens)-1));
handles = refresh_genlist(handles,hObject,0);
% set(handles.checkbox_analyzenew,'Value',1);
guidata(hObject, handles);

function addfolder_Callback(hObject, eventdata, handles)
% SELECT A FOLDER AND ADD ALL CONTAINED FEATURE FILES 
% TO THE FILELIST
PathName = uigetdir;
if PathName,
    handles.vars.PathName = [PathName handles.vars.params.slash];
    FileName = dir([handles.vars.PathName '*.feat']);
    handles = addfiles(FileName,handles,hObject);
    enable_all('on', hObject, handles);
end
% Set checkboxes to suggest reanalyzing all data
handles.vars.params.analyze_new = 1; 
set(handles.checkbox_analyzenew,'Value',handles.vars.params.analyze_new);
guidata(hObject, handles);

function addfiles_Callback(hObject, eventdata, handles)
% SELECT AT LEAST ONE FEATURE FILE AND ADD TO FILELIST
[FileName,PathName] = uigetfile('*.feat','Select Feature Files','MultiSelect','On');
fn = 0; if iscell(FileName), fn = 1; else fn = 0; end; if fn == 0, if FileName, fn = 1; end; end
if fn,
    handles.vars.PathName = PathName;
    handles = addfiles(FileName,handles,hObject);
    enable_all('on', hObject, handles);
end
handles.vars.params.analyze_new = 1; set(handles.checkbox_analyzenew,'Value',handles.vars.params.analyze_new);
guidata(hObject, handles);

function handles = removefile_Callback(hObject, eventdata, handles)
% REMOVE SELECTED FILE FROM FILE LIST AND 
% UPDATE THE OTHER LISTS ACCORDINGLY
handles.vars.file_index = get(handles.listbox_files,'Value');
ind_keep = [1:handles.vars.file_index-1 handles.vars.file_index+1:numel(handles.vars.file)];
handles.vars.file = handles.vars.file(ind_keep); 
handles.vars.file_list = handles.vars.file_list(ind_keep);
handles.vars.file_index = max(handles.vars.file_index-1,1);
handles.vars.chamber_index = handles.vars.chamber_index(ind_keep);
handles = refresh_genlist(handles,hObject,1);
if handles.vars.file_index && numel(handles.vars.gens),
    set(handles.listbox_files, 'String', handles.vars.file_list);
    set(handles.listbox_chambers, 'String', handles.vars.file(handles.vars.file_index).Chambers);
    set(handles.edit_genotype, 'String', handles.vars.file(handles.vars.file_index).Genotypes{handles.vars.chamber_index(handles.vars.file_index)});
else
    set(handles.listbox_files, 'String', '');
    set(handles.listbox_chambers, 'String', '');
    set(handles.edit_genotype, 'String', '');
end
handles.vars.params.analyze_new = 1; set(handles.checkbox_analyzenew,'Value',handles.vars.params.analyze_new);
handles = savedat(hObject, handles);

function handles = remove_gen_Callback(hObject, eventdata, handles)
% REMOVE SELECTED GENOTYPE OR CHAMBER (ARENA) FROM 
% GENOTYPE/CHAMBER LIST AND UPDATE 
% THE OTHER LISTS ACCORDINGLY
% Genotype or chamber selected?
if handles.vars.editgen,
    % Remove one genotype
    gen = handles.vars.gens{get(handles.listbox_genotypes,'Value')};
    current_file_index = handles.vars.file_index;
    handles.vars.file_index = 0;
    tmp_file = handles.vars.file; tmp_chamb = handles.vars.chamber_index;
    tmp_filel = handles.vars.file_list;
    handles.vars.file = []; handles.vars.chamber_index = []; handles.vars.file_list = [];
    for j=1:numel(tmp_file),
        in = [];
        for jj=1:numel(tmp_file(j).Genotypes),
            if ~strcmp(tmp_file(j).Genotypes{jj},gen),
                in = [in jj];
            end
        end
        if ~isempty(in),
            handles.vars.file_index = handles.vars.file_index + 1;
            handles.vars.file(handles.vars.file_index).FileName = tmp_file(j).FileName; % use curly braces to index into cell array
            handles.vars.file(handles.vars.file_index).PathName = tmp_file(j).PathName;
            handles.vars.file(handles.vars.file_index).Chambers = tmp_file(j).Chambers(in); % use curly braces to index into cell array
            for jj=1:numel(in),
                handles.vars.file(handles.vars.file_index).Genotypes{jj} = tmp_file(j).Genotypes{in(jj)};
            end
            handles.vars.chamber_index(handles.vars.file_index) = tmp_chamb(j);
            handles.vars.file_list{handles.vars.file_index} = tmp_filel{j};
        end
    end
    handles.vars.file_index = min(current_file_index,handles.vars.file_index);
    if handles.vars.file_index,
        set(handles.listbox_files, 'String', handles.vars.file_list, ...
            'Value',handles.vars.file_index);
        set(handles.listbox_chambers, 'String', handles.vars.file(handles.vars.file_index).Chambers);
        set(handles.edit_genotype, 'String', handles.vars.file(handles.vars.file_index).Genotypes{handles.vars.chamber_index(handles.vars.file_index)});
        guidata(hObject, handles);  % Update handles structure
    end
elseif handles.vars.editchamber,
    % Remove one chamber (or complete genotype in case it was
    % only one chamber remaining
    [chamb,chambi] = setdiff(handles.vars.file(handles.vars.file_index).Chambers,...
        handles.vars.file(handles.vars.file_index).Chambers(handles.vars.chamber_index(handles.vars.file_index)));
    handles.vars.file(handles.vars.file_index).Chambers = chamb;
    handles.vars.chamber_index(handles.vars.file_index) = min(handles.vars.chamber_index(handles.vars.file_index),numel(chamb));
    if handles.vars.chamber_index(handles.vars.file_index),
        tmp_gens = handles.vars.file(handles.vars.file_index).Genotypes;
        handles.vars.file(handles.vars.file_index).Genotypes = cell(1,numel(chamb));
        for jj=1:numel(chamb),
            handles.vars.file(handles.vars.file_index).Genotypes{jj} = tmp_gens{chambi(jj)};
        end
        set(handles.listbox_chambers, 'String', handles.vars.file(handles.vars.file_index).Chambers, ...
            'Value', handles.vars.chamber_index(handles.vars.file_index));
        guidata(hObject, handles);  % Update handles structure
    else
        handles = removefile_Callback(hObject, eventdata, handles);
    end
end
% Update all GUI lists
handles = refresh_genlist(handles,hObject,1);
in = get(handles.listbox_genotypes,'Value');
set(handles.listbox_genotypes, 'String', handles.vars.gens, 'Value', in);
% set(handles.listbox_genotypes, 'Value', 1);
set(handles.listbox_n, 'Value', in);
if ~handles.vars.file_index,
    set(handles.listbox_files, 'String', '');
    set(handles.listbox_chambers, 'String', '');
    set(handles.edit_genotype, 'String', '');
    set(handles.listbox_genotypes, 'String', '');
    set(handles.listbox_n, 'String', '');
end
handles.vars.params.analyze_new = 1; set(handles.checkbox_analyzenew,'Value',handles.vars.params.analyze_new);
handles = savedat(hObject, handles);

function pushbutton_saveas_Callback(hObject, eventdata, handles)
% ASSIGN PATH AND NAME OF CONFIG FILE BY
% SELECTION USING A GUI
[saveas1,SavePathName1] = uiputfile({'*.cfg','Configuration Files (*.cfg)';'*.*','All Files (*.*)'},...
    'Save Configuration',[handles.vars.SavePathName handles.vars.saveas]);
if saveas1,
    handles.vars.saveas = saveas1; handles.vars.SavePathName = SavePathName1;
end
set(handles.edit_save,'String',[handles.vars.SavePathName handles.vars.saveas]);
handles = savedat(hObject, handles);

function edit_save_Callback(hObject, eventdata, handles)
% ASSIGN PATH AND NAME OF CONFIG FILE BY
% EDITING THE SAVE EDITBOX
pathfile = get(handles.edit_save,'String');
in = strfind(pathfile,handles.vars.params.slash); in = in(end);
handles.vars.SavePathName = pathfile(1:in);
handles.vars.saveas = pathfile(in+1:end);
handles = savedat(hObject, handles);

function analyze_Callback(hObject, eventdata, handles)
% ANALYZE DATA
% Prepare needed variables for call of the main routine 
% 'main_analysis'
fnames = []; sorted_names = []; cnt = 0;
for j=1:numel(handles.vars.file),
    chamb = handles.vars.file(j).Chambers;
    for jj=1:numel(chamb),
        cnt = cnt + 1;
        fnames{cnt} = handles.vars.file(j).Genotypes{jj};
        sorted_names{cnt} = [handles.vars.file(j).PathName ...
            handles.vars.file(j).FileName '_' num2str(chamb(jj)) '.feat'];
%         paths{cnt} = handles.vars.file(j).PathName;
    end
end
mov_count = handles.vars.n; gen_ident = handles.vars.gens;
genotypes = zeros(1,numel(mov_count));
for j=1:numel(mov_count),
    genotypes(find(strcmp(fnames,gen_ident{j}))) = j;
end
for i=1:length(gen_ident), 
    gen_ident{i} = gen_ident{i}(setdiff(1:length(gen_ident{i}),strfind(gen_ident{i},'_')));
end
gen_ident1 = gen_ident;
% activate the following line in case the current plot list should be 
% reseted for each 'analyze' call
% handles.vars.params.fid = 0; 
% profile on; % performance profiler
% THIS IS THE MAIN CALL
handles.vars.params.fid = main_analysis(handles.vars.SavePathName,sorted_names,mov_count,gen_ident,gen_ident1,genotypes,handles.vars.chamber,handles.vars.params);
% p = profile('info');

% Deactivate some checkboxes after analysis
handles.vars.params.analyze_new = 0; handles.vars.params.feat_read_new = 0;
% set(handles.checkbox_readnew,'Value',handles.vars.params.feat_read_new);
set(handles.checkbox_analyzenew,'Value',handles.vars.params.analyze_new);
handles = savedat(hObject, handles);


function auto_Callback(hObject, eventdata, handles)
% TRY TO ENCODE GENOTYPE NAMES FROM FILE NAMES OR
% PATH NAMES AUTOMATICALLY 
for j=1:numel(handles.vars.file),
    chamb = handles.vars.file(j).Chambers;
    for jj=1:numel(chamb),
        fname = handles.vars.file(j).FileName;
        if handles.vars.auto || get(handles.checkbox_genfolder,'Value'),
            if ~get(handles.checkbox_genfolder,'Value'),
                % Genotypes from file names
                gen = genname(fname,chamb(jj));
            else
                % Genotypes from path names
                pth = handles.vars.file(j).PathName;
                in = strfind(pth,handles.vars.params.slash); in = in(end-1:end);
                gen = pth(in(1)+1:in(2)-1);
            end
        else
            gen = fname;
        end
        handles.vars.file(j).Genotypes{jj} = gen;
    end
end
if handles.vars.auto, handles.vars.auto = 0; else handles.vars.auto = 1; end
% Update all GUI lists
if get(handles.checkbox_genfolder,'Value'), set(handles.checkbox_genfolder,'Value',0); end
handles = refresh_genlist(handles,hObject,1);
set(handles.listbox_files, 'String', handles.vars.file_list);
set(handles.listbox_chambers, 'String', handles.vars.file(handles.vars.file_index).Chambers);
set(handles.edit_genotype, 'String', handles.vars.file(handles.vars.file_index).Genotypes{handles.vars.chamber_index(handles.vars.file_index)});
guidata(hObject, handles);

function about_Callback(hObject, eventdata, handles)
% DISPLAY GNU MESSAGE
msg = GNU_message_analysis;
h = msgbox(msg);
uiwait(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECKBOX CALLS

function checkbox_genfolder_Callback(hObject, eventdata, handles)

function checkbox_winlos_Callback(hObject, eventdata, handles)
handles.vars.params.winlos = get(handles.checkbox_winlos,'Value');
handles = savedat(hObject, handles);

function checkbox_courtship_Callback(hObject, eventdata, handles)
handles.vars.params.courtship = get(handles.checkbox_courtship,'Value');
handles.vars.params.analyze_new = 1;
set(handles.checkbox_analyzenew,'Value',handles.vars.params.analyze_new);
handles = savedat(hObject, handles);

function checkbox_oneobj_Callback(hObject, eventdata, handles)
handles.vars.params.oneobj = get(handles.checkbox_oneobj,'Value');
handles.vars.params.analyze_new = 1;
set(handles.checkbox_analyzenew,'Value',handles.vars.params.analyze_new);
handles = savedat(hObject, handles);

function checkbox_autoscale_Callback(hObject, eventdata, handles)
handles.vars.params.autoscale = get(handles.checkbox_autoscale,'Value');
handles = savedat(hObject, handles);

function edit_rows_Callback(hObject, eventdata, handles)
handles.vars.params.k(1) = abs(str2double(get(handles.edit_rows,'String')));
handles.vars.params.k(1) = max(handles.vars.params.k(1),1);
handles.vars.params.k(1) = min(handles.vars.params.k(1),numel(handles.vars.gens));
handles = savedat(hObject, handles);

function edit_cols_Callback(hObject, eventdata, handles)
handles.vars.params.k(2) = abs(str2double(get(handles.edit_cols,'String')));
handles.vars.params.k(2) = max(handles.vars.params.k(2),1);
handles.vars.params.k(2) = min(handles.vars.params.k(2),numel(handles.vars.gens));
handles = savedat(hObject, handles);

function checkbox_pdf_Callback(hObject, eventdata, handles)
handles.vars.params.pdf = get(handles.checkbox_pdf,'Value');
handles = savedat(hObject, handles);

function checkbox_text_Callback(hObject, eventdata, handles)
handles.vars.params.bool_xls = get(handles.checkbox_text,'Value');
handles = savedat(hObject, handles);

function checkbox_boxplot_Callback(hObject, eventdata, handles)
handles.vars.params.boxplot = get(handles.checkbox_boxplot,'Value');
handles = savedat(hObject, handles);

function checkbox_bonf_Callback(hObject, eventdata, handles)
handles.vars.params.bonf = get(handles.checkbox_bonf,'Value');
handles = savedat(hObject, handles);

function checkbox_meansem_Callback(hObject, eventdata, handles)
handles.vars.params.boxplot = 1-get(handles.checkbox_boxplot,'Value');
handles = savedat(hObject, handles);

% function checkbox_readnew_Callback(hObject, eventdata, handles)
% handles.vars.params.feat_read_new = get(handles.checkbox_readnew,'Value');
% if handles.vars.params.feat_read_new,
%     set(handles.checkbox_analyzenew,'Value',1);
%     handles.vars.params.analyze_new = 1;
% end
% handles = savedat(hObject, handles);

function checkbox_analyzenew_Callback(hObject, eventdata, handles)
handles.vars.params.analyze_new = get(handles.checkbox_analyzenew,'Value');
% Deactivate in case you want to separate 'read_feat' and 'analyze'
if handles.vars.params.analyze_new,
%     set(handles.checkbox_readnew,'Value',1);
    handles.vars.params.feat_read_new = 1;
end
%
handles = savedat(hObject, handles);

function checkbox_border_Callback(hObject, eventdata, handles)
handles.vars.params.border = get(handles.checkbox_border,'Value');
handles.vars.params.analyze_new = 1;
set(handles.checkbox_analyzenew,'Value',handles.vars.params.analyze_new);
if handles.vars.params.border,
    set(handles.checkbox_circular,'Value',0);
    handles.vars.params.radius = -99;
end
handles = savedat(hObject, handles);

function edit_axisfontsize_Callback(hObject, eventdata, handles)
handles.vars.params.axisfontsize = abs(str2double(get(handles.edit_axisfontsize,'String')));
handles.vars.params.axisfontsize = max(handles.vars.params.axisfontsize,4);
handles = savedat(hObject, handles);

function edit_chamberheight_Callback(hObject, eventdata, handles)
handles.vars.chamber.height = abs(str2double(get(handles.edit_chamberheight,'String')));
handles.vars.chamber.height = max(handles.vars.chamber.height,5);
handles = savedat(hObject, handles);

function edit_chamberwidth_Callback(hObject, eventdata, handles)
handles.vars.chamber.width = abs(str2double(get(handles.edit_chamberwidth,'String')));
handles.vars.chamber.width = max(handles.vars.chamber.width,5);
handles = savedat(hObject, handles);

function edit_timelimit_Callback(hObject, eventdata, handles)
handles.vars.params.max_frames = abs(str2double(get(handles.edit_timelimit,'String'))) * 60; %time limit [s]
handles.vars.params.k(2) = max(handles.vars.params.max_frames,60);
handles = savedat(hObject, handles);

function clean_Callback(hObject, eventdata, handles)
if get(handles.clean,'Value'),
    if handles.vars.params.fid > 0, 
        for i=1:handles.vars.params.fid,
            try
                close(i);
            catch err
            end
        end
        handles.vars.params.fid = 0;
    end
    set(handles.clean,'Value',0);
    guidata(hObject, handles);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION CALLS

% --- Executes during object creation, after setting all properties.
function listbox_genotypes_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function listbox_files_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function edit_genotype_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function listbox_n_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
% --- Executes during object creation, after setting all properties.
function edit_save_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function listbox_chambers_CreateFcn(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function edit_timelimit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function edit_rows_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function edit_cols_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function edit_axisfontsize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function edit_chamberwidth_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function edit_chamberheight_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function heatmap_position_Callback(hObject, eventdata, handles)
% CHECKBOX TO PLOT FLY POSITION HEATMAPS
handles.vars.params.plots.heatmaps.position = get(handles.heatmap_position,'Value');
handles.vars.params.selected = 'heatmap_position';
% Enable/Disable corresponding listbox for further selections & range editbox
if handles.vars.params.plots.heatmaps.position, onoff = 'on'; else onoff = 'off'; end
handles = setplots(onoff, hObject, handles);
set(handles.edit_posheat,'Enable',onoff);
handles = savedat(hObject, handles);

function heatmap_velocity_Callback(hObject, eventdata, handles)
% CHECKBOX TO PLOT FLY VELOCITY HEATMAPS
handles.vars.params.plots.heatmaps.velocity = get(handles.heatmap_velocity,'Value');
handles.vars.params.selected = 'heatmap_velocity';
% Enable/Disable corresponding listbox for further selections & range editbox
if handles.vars.params.plots.heatmaps.velocity, onoff = 'on'; else onoff = 'off'; end
handles = setplots(onoff, hObject, handles);
set(handles.edit_velheat,'Enable',onoff);
handles = savedat(hObject, handles);

function heatmap_ext_leftwing_Callback(hObject, eventdata, handles)
% CHECKBOX TO PLOT FLY LEFT WING EXTENSION HEATMAPS
handles.vars.params.plots.heatmaps.ext_leftwing = get(handles.heatmap_ext_leftwing,'Value');
handles.vars.params.selected = 'heatmap_wing';
% Enable/Disable corresponding listbox for further selections & range editbox
if handles.vars.params.plots.heatmaps.ext_leftwing || handles.vars.params.plots.heatmaps.ext_rightwing || ...
   handles.vars.params.plots.heatmaps.ext_onewing || handles.vars.params.plots.heatmaps.ext_diffleftrightwing, 
    onoff = 'on'; 
else
    onoff = 'off';
end
handles = setplots(onoff, hObject, handles);
set(handles.edit_wingheat,'Enable',onoff);
handles = savedat(hObject, handles);

function heatmap_ext_rightwing_Callback(hObject, eventdata, handles)
% CHECKBOX TO PLOT FLY RIGHT WING EXTENSION HEATMAPS
handles.vars.params.plots.heatmaps.ext_rightwing = get(handles.heatmap_ext_rightwing,'Value');
handles.vars.params.selected = 'heatmap_wing';
% Enable/Disable corresponding listbox for further selections & range editbox
if handles.vars.params.plots.heatmaps.ext_leftwing || handles.vars.params.plots.heatmaps.ext_rightwing || ...
   handles.vars.params.plots.heatmaps.ext_onewing || handles.vars.params.plots.heatmaps.ext_diffleftrightwing, 
    onoff = 'on'; 
else
    onoff = 'off';
end
handles = setplots(onoff, hObject, handles);
set(handles.edit_wingheat,'Enable',onoff);
handles = savedat(hObject, handles);

function heatmap_ext_diffleftrightwing_Callback(hObject, eventdata, handles)
% CHECKBOX TO PLOT FLY WING EXTENSION DIFFERENCE 
% (LEFT-RIGHT) HEATMAPS
handles.vars.params.plots.heatmaps.ext_diffleftrightwing = get(handles.heatmap_ext_diffleftrightwing,'Value');
handles.vars.params.selected = 'heatmap_wing';
% Enable/Disable corresponding listbox for further selections & range editbox
if handles.vars.params.plots.heatmaps.ext_leftwing || handles.vars.params.plots.heatmaps.ext_rightwing || ...
   handles.vars.params.plots.heatmaps.ext_onewing || handles.vars.params.plots.heatmaps.ext_diffleftrightwing, 
    onoff = 'on'; 
else
    onoff = 'off';
end
handles = setplots(onoff, hObject, handles);
set(handles.edit_wingheat,'Enable',onoff);
handles = savedat(hObject, handles);

function heatmap_ext_onewing_Callback(hObject, eventdata, handles)
% CHECKBOX TO PLOT FLY ONE WING EXTENSION HEATMAPS
handles.vars.params.plots.heatmaps.ext_onewing = get(handles.heatmap_ext_onewing,'Value');
handles.vars.params.selected = 'heatmap_wing';
% Enable/Disable corresponding listbox for further selections & range editbox
if handles.vars.params.plots.heatmaps.ext_leftwing || handles.vars.params.plots.heatmaps.ext_rightwing || ...
   handles.vars.params.plots.heatmaps.ext_onewing || handles.vars.params.plots.heatmaps.ext_diffleftrightwing, 
    onoff = 'on'; 
else
    onoff = 'off';
end
handles = setplots(onoff, hObject, handles);
set(handles.edit_wingheat,'Enable',onoff);
handles = savedat(hObject, handles);

function heatmap_wingrel_Callback(hObject, eventdata, handles)
% CHECKBOX TO PLOT RELATIVE FLY POSITION HEATMAPS 
% WHILE THE OPPONENT EXTENDS ONE WING
handles.vars.params.plots.heatmaps.wingrel = get(handles.heatmap_wingrel,'Value');
handles.vars.params.selected = 'heatmap_wingrel';
% Enable/Disable corresponding listbox for further selections & range editbox
if handles.vars.params.plots.heatmaps.wingrel, 
    onoff = 'on'; 
else
    onoff = 'off';
end
handles = setplots(onoff, hObject, handles);
set(handles.edit_wingheat,'Enable',onoff);
handles = savedat(hObject, handles);

function heatmap_chasing_Callback(hObject, eventdata, handles)
% CHECKBOX TO PLOT CHASING HEATMAPS
handles.vars.params.plots.heatmaps.chasing = get(handles.heatmap_chasing,'Value');
handles.vars.params.selected = 'heatmap_chase';
% Enable/Disable corresponding listbox for further selections & range editbox
if handles.vars.params.plots.heatmaps.chasing, onoff = 'on'; else onoff = 'off'; end
handles = setplots(onoff, hObject, handles);
set(handles.edit_chaseheat,'Enable',onoff);
handles = savedat(hObject, handles);

function heatmap_circling_Callback(hObject, eventdata, handles)
% CHECKBOX TO PLOT CHASING HEATMAPS
handles.vars.params.plots.heatmaps.circling = get(handles.heatmap_circling,'Value');
handles.vars.params.selected = 'heatmap_circl';
% Enable/Disable corresponding listbox for further selections & range editbox
if handles.vars.params.plots.heatmaps.circling, onoff = 'on'; else onoff = 'off'; end
handles = setplots(onoff, hObject, handles);
set(handles.edit_circlheat,'Enable',onoff);
handles = savedat(hObject, handles);

function heatmap_lunging_Callback(hObject, eventdata, handles)
% CHECKBOX TO PLOT LUNGING HEATMAPS
handles.vars.params.plots.heatmaps.lunging = get(handles.heatmap_lunging,'Value');
handles.vars.params.selected = 'heatmap_lunge';
% Enable/Disable corresponding listbox for further selections & range editbox
if handles.vars.params.plots.heatmaps.lunging, onoff = 'on'; else onoff = 'off'; end
handles = setplots(onoff, hObject, handles);
set(handles.edit_lungeheat,'Enable',onoff);
handles = savedat(hObject, handles);

function heatmap_wingthreat_Callback(hObject, eventdata, handles)
% CHECKBOX TO PLOT WING THREAT HEATMAPS
handles.vars.params.plots.heatmaps.wingthreat = get(handles.heatmap_wingthreat,'Value');
handles.vars.params.selected = 'heatmap_wingthreat';
% Enable/Disable corresponding listbox for further selections & range editbox
if handles.vars.params.plots.heatmaps.wingthreat, onoff = 'on'; else onoff = 'off'; end
handles = setplots(onoff, hObject, handles);
set(handles.edit_wingthreatheat,'Enable',onoff);
handles = savedat(hObject, handles);

function heatmap_tussling_Callback(hObject, eventdata, handles)
% CHECKBOX TO PLOT TUSSLING HEATMAPS
handles.vars.params.plots.heatmaps.tussling = get(handles.heatmap_tussling,'Value');
handles.vars.params.selected = 'heatmap_tussl';
% Enable/Disable corresponding listbox for further selections & range editbox
if handles.vars.params.plots.heatmaps.tussling, onoff = 'on'; else onoff = 'off'; end
handles = setplots(onoff, hObject, handles);
set(handles.edit_tusslheat,'Enable',onoff);
handles = savedat(hObject, handles);

function heatmap_selAll_Callback(hObject, eventdata, handles)
% CHECKBOX TO SELECT/UNSELECT ALL HEATMAP PLOTS
if get(handles.heatmap_selAll,'Value'), val = 1; onoff = 'on'; else val = 0; onoff = 'off'; end
handles.vars.params.plots.heatmaps.position = val;
handles.vars.params.plots.heatmaps.velocity = val;
handles.vars.params.plots.heatmaps.ext_leftwing = val;
handles.vars.params.plots.heatmaps.ext_rightwing = val;
handles.vars.params.plots.heatmaps.ext_diffleftrightwing = val;
handles.vars.params.plots.heatmaps.ext_onewing = val;
handles.vars.params.plots.heatmaps.chasing = val;
handles.vars.params.plots.heatmaps.circling = val;
handles.vars.params.plots.heatmaps.tussling = val;
handles.vars.params.plots.heatmaps.lunging = val;
handles.vars.params.plots.heatmaps.wingthreat = val;
handles.vars.params.plots.heatmaps.wingrel = val;

if val, dummyvec = 1:20; else dummyvec = 0; end
fieldn = fieldnames(handles.vars.params.plots.heatmaps);
for i=1:length(fieldn),
    if numel(strfind(fieldn{i},'vec')),
        handles.vars.params.plots.heatmaps.(fieldn{i}) = dummyvec;
    end
end
set(handles.edit_posheat,'Enable',onoff); set(handles.edit_velheat,'Enable',onoff);
set(handles.edit_wingheat,'Enable',onoff); set(handles.edit_wingheat,'Enable',onoff);
set(handles.edit_wingheat,'Enable',onoff); set(handles.edit_wingheat,'Enable',onoff);
set(handles.edit_wingheat,'Enable',onoff); set(handles.edit_chaseheat,'Enable',onoff);
set(handles.edit_circlheat,'Enable',onoff); set(handles.edit_lungeheat,'Enable',onoff);
set(handles.edit_wingthreatheat,'Enable',onoff); set(handles.edit_tusslheat,'Enable',onoff);

% set(handles.heatmap_position,'Value',onoff); set(handles.heatmap_velocity,'Value',onoff);
% set(handles.heatmap_ext_leftwing,'Value',onoff); set(handles.heatmap_ext_rightwing,'Value',onoff);
% set(handles.heatmap_ext_diffleftrightwing,'Value',onoff); set(handles.heatmap_ext_onewing,'Value',onoff);
% set(handles.heatmap_chasing,'Value',onoff); set(handles.heatmap_circling,'Value',onoff);
% set(handles.heatmap_tussling,'Value',onoff); set(handles.heatmap_lunging,'Value',onoff);
% set(handles.heatmap_wingrel,'Value',onoff);
handles = savedat(hObject, handles);


function cumul_time_Callback(hObject, eventdata, handles)
% CHECKBOX TO SELECT PLOTTING DATA AS 
% CUMULATIVE TIME OR OCCURRENCES
handles.vars.params.tim = get(handles.cumul_time,'Value');
handles = savedat(hObject, handles);


function stat_chasing_Callback(hObject, eventdata, handles)
% CHECKBOX TO PLOT CHASING STATISTICS
handles.vars.params.plots.stat.chasing = get(handles.stat_chasing,'Value');
if handles.vars.params.plots.stat.chasing, onoff = 'on'; else onoff = 'off'; end
handles = set_selplots_off(hObject, handles);
if handles.vars.params.plots.movieclips,
    % In case the 'movieclips' checkbox is activated allow only this
    % action to be selected
    if handles.vars.params.plots.stat.chasing, val = 1; else val = 0; end 
    handles.vars.params.plots.stat.chasing = val;
    handles.vars.params.plots.stat.lunging = 0;
    handles.vars.params.plots.stat.wingthreat = 0;
    handles.vars.params.plots.stat.onewing = 0; 
    handles.vars.params.plots.stat.circling = 0;
    handles.vars.params.plots.stat.copulation = 0;
    handles.vars.params.plots.stat.tussl = 0; 
    set(handles.edit_chasestat,'Enable','off');
else
    set(handles.edit_chasestat,'Enable',onoff);
end
handles = savedat(hObject, handles);

function stat_walkingdist_Callback(hObject, eventdata, handles)
% CHECKBOX TO PLOT WALKING DISTANCE STATISTICS
handles.vars.params.plots.stat.traveldist = get(handles.stat_walkingdist,'Value');
if handles.vars.params.plots.stat.traveldist, onoff = 'on'; else onoff = 'off'; end
set(handles.edit_walkstat,'Enable',onoff);
handles = set_selplots_off(hObject, handles);
handles = savedat(hObject, handles);

function stat_dist_Callback(hObject, eventdata, handles)
% CHECKBOX TO PLOT FLY DISTANCE STATISTICS
handles.vars.params.plots.stat.dist = get(handles.stat_dist,'Value');
handles.vars.params.selected = 'stat_flydist';
% Enable/Disable corresponding listbox for further selections & range editbox
if handles.vars.params.plots.stat.dist, onoff = 'on'; else onoff = 'off'; end
handles = setplots(onoff, hObject, handles);
set(handles.edit_flydiststat,'Enable',onoff);
handles = savedat(hObject, handles);

function stat_movestop_Callback(hObject, eventdata, handles)
% CHECKBOX TO PLOT MOVE/STOP STATISTICS
handles.vars.params.plots.stat.movestop = get(handles.stat_movestop,'Value');
handles = set_selplots_off(hObject, handles);
handles = savedat(hObject, handles);

function stat_bodysize_Callback(hObject, eventdata, handles)
% CHECKBOX TO PLOT FLY BODY SIZE STATISTICS
handles.vars.params.plots.stat.body = get(handles.stat_bodysize,'Value');
handles.vars.params.selected = 'stat_body';
% Enable/Disable corresponding listbox for further selections & range editbox
if handles.vars.params.plots.stat.body, onoff = 'on'; else onoff = 'off'; end
handles = setplots(onoff, hObject, handles);
handles = savedat(hObject, handles);

function stat_velocity_Callback(hObject, eventdata, handles)
% CHECKBOX TO PLOT FLY VELOCITY STATISTICS
handles.vars.params.plots.stat.velocity = get(handles.stat_velocity,'Value');
handles.vars.params.selected = 'stat_vel';
% Enable/Disable corresponding listbox for further selections & range editbox
if handles.vars.params.plots.stat.velocity, onoff = 'on'; else onoff = 'off'; end
handles = setplots(onoff, hObject, handles);
set(handles.edit_velstat,'Enable',onoff);
handles = savedat(hObject, handles);

function stat_dendro_Callback(hObject, eventdata, handles)
% CHECKBOX TO PLOT A DENDROGRAM
handles.vars.params.plots.stat.dendro = get(handles.stat_dendro,'Value');
handles = set_selplots_off(hObject, handles);
handles = savedat(hObject, handles);

function timeseries_Callback(hObject, eventdata, handles)
% CHECKBOX TO PLOT TIME SERIES
handles.vars.params.plots.timeseries = get(handles.timeseries,'Value');
handles.vars.params.selected = 'timeseries';
% Enable/Disable corresponding listbox for further selections & range editbox
if handles.vars.params.plots.timeseries, onoff = 'on'; else onoff = 'off'; end
handles = setplots(onoff, hObject, handles);
handles = savedat(hObject, handles);

function stat_actionrel_Callback(hObject, eventdata, handles)
% CHECKBOX TO PLOT ACTION RELATION STATISTICS
handles.vars.params.plots.stat.actionrel = get(handles.stat_actionrel,'Value');
handles.vars.params.selected = 'stat_actionrel';
% Enable/Disable corresponding listbox for further selections & range editbox
if handles.vars.params.plots.stat.actionrel, onoff = 'on'; else onoff = 'off'; end
handles = setplots(onoff, hObject, handles);
handles = savedat(hObject, handles);

function stat_tussl_Callback(hObject, eventdata, handles)
% CHECKBOX TO PLOT TUSSLING STATISTICS
handles.vars.params.plots.stat.tussl = get(handles.stat_tussl,'Value');
if handles.vars.params.plots.stat.tussl, onoff = 'on'; else onoff = 'off'; end
if handles.vars.params.plots.movieclips,
    % In case the 'movieclips' checkbox is activated allow only this
    % action to be selected
    if handles.vars.params.plots.stat.tussl, val = 1; else val = 0; end 
    handles.vars.params.plots.stat.chasing = 0;
    handles.vars.params.plots.stat.lunging = 0;
    handles.vars.params.plots.stat.wingthreat = 0;
    handles.vars.params.plots.stat.onewing = 0; 
    handles.vars.params.plots.stat.circling = 0;
    handles.vars.params.plots.stat.copulation = 0;
    handles.vars.params.plots.stat.tussl = val; 
    handles.vars.params.selected = '';
    set(handles.edit_tusslstat,'Enable','off');
else
    handles.vars.params.selected = 'stat_tussl';
    set(handles.edit_tusslstat,'Enable',onoff);
end
% Enable/Disable corresponding listbox for further selections & range editbox
handles = setplots(onoff, hObject, handles);
handles = savedat(hObject, handles);

function stat_lunging_Callback(hObject, eventdata, handles)
% CHECKBOX TO PLOT LUNGING STATISTICS
handles.vars.params.plots.stat.lunging = get(handles.stat_lunging,'Value');
if handles.vars.params.plots.stat.lunging, onoff = 'on'; else onoff = 'off'; end
if handles.vars.params.plots.movieclips,
    % In case the 'movieclips' checkbox is activated allow only this
    % action to be selected
    if handles.vars.params.plots.stat.lunging, val = 1; else val = 0; end 
    handles.vars.params.plots.stat.chasing = 0;
    handles.vars.params.plots.stat.lunging = val;
    handles.vars.params.plots.stat.wingthreat = 0;
    handles.vars.params.plots.stat.onewing = 0; 
    handles.vars.params.plots.stat.circling = 0;
    handles.vars.params.plots.stat.copulation = 0;
    handles.vars.params.plots.stat.tussl = 0; 
    handles.vars.params.selected = '';
    set(handles.edit_lungestat,'Enable','off');
else
    handles.vars.params.selected = 'stat_lunge';
    set(handles.edit_lungestat,'Enable',onoff);
end    
% Enable/Disable corresponding listbox for further selections & range editbox
handles = setplots(onoff, hObject, handles);
handles = savedat(hObject, handles);

function stat_wingthreat_Callback(hObject, eventdata, handles)
% CHECKBOX TO PLOT WING THREAT STATISTICS
handles.vars.params.plots.stat.wingthreat = get(handles.stat_wingthreat,'Value');
handles.vars.params.selected = 'stat_wingthreat';
if handles.vars.params.plots.stat.wingthreat, onoff = 'on'; else onoff = 'off'; end
% handles = setplots(onoff, hObject, handles);
handles = set_selplots_off(hObject, handles);
if handles.vars.params.plots.movieclips,
    % In case the 'movieclips' checkbox is activated allow only this
    % action to be selected
    if handles.vars.params.plots.stat.wingthreat, val = 1; else val = 0; end 
    handles.vars.params.plots.stat.chasing = 0;
    handles.vars.params.plots.stat.wingthreat = val;
    handles.vars.params.plots.stat.lunging = 0;
    handles.vars.params.plots.stat.onewing = 0; 
    handles.vars.params.plots.stat.circling = 0;
    handles.vars.params.plots.stat.copulation = 0;
    handles.vars.params.plots.stat.tussl = 0; 
    set(handles.edit_wingthreatstat,'Enable','off');
else
    set(handles.edit_wingthreatstat,'Enable',onoff);
end
handles = savedat(hObject, handles);

function stat_scat_Callback(hObject, eventdata, handles)
% CHECKBOX FOR SCATTERPLOTS
handles.vars.params.plots.stat.scat = get(handles.stat_scat,'Value');
handles.vars.params.selected = 'stat_scat';
% Enable/Disable corresponding listbox for further selections & range editbox
if handles.vars.params.plots.stat.scat, onoff = 'on'; else onoff = 'off'; end
handles = setplots(onoff, hObject, handles);
handles = savedat(hObject, handles);

function stat_leftwing_Callback(hObject, eventdata, handles)
% CHECKBOX TO PLOT LEFT WING EXTENSION STATISTICS
handles.vars.params.plots.stat.leftwing = get(handles.stat_leftwing,'Value');
handles.vars.params.selected = 'stat_wing';
% Enable/Disable corresponding listbox for further selections & range editbox
if handles.vars.params.plots.stat.leftwing, onoff = 'on'; else onoff = 'off'; end
handles = setplots(onoff, hObject, handles);
set(handles.edit_wingstat,'Enable',onoff);
handles = savedat(hObject, handles);

function stat_rightwing_Callback(hObject, eventdata, handles)
% CHECKBOX TO PLOT RIGHT WING EXTENSION STATISTICS
handles.vars.params.plots.stat.rightwing = get(handles.stat_rightwing,'Value');
handles.vars.params.selected = 'stat_wing';
% Enable/Disable corresponding listbox for further selections & range editbox
if handles.vars.params.plots.stat.rightwing, onoff = 'on'; else onoff = 'off'; end
handles = setplots(onoff, hObject, handles);
set(handles.edit_wingstat,'Enable',onoff);
handles = savedat(hObject, handles);

function stat_onewing_Callback(hObject, eventdata, handles)
% CHECKBOX TO PLOT ONE WING EXTENSION STATISTICS
handles.vars.params.plots.stat.onewing = get(handles.stat_onewing,'Value');
if handles.vars.params.plots.stat.onewing, onoff = 'on'; else onoff = 'off'; end
if handles.vars.params.plots.movieclips,
    % In case the 'movieclips' checkbox is activated allow only this
    % action to be selected
    if handles.vars.params.plots.stat.onewing, val = 1; else val = 0; end 
    handles.vars.params.plots.stat.chasing = 0;
    handles.vars.params.plots.stat.lunging = 0;
    handles.vars.params.plots.stat.wingthreat = 0;
    handles.vars.params.plots.stat.onewing = val; 
    handles.vars.params.plots.stat.circling = 0;
    handles.vars.params.plots.stat.copulation = 0;
    handles.vars.params.plots.stat.tussl = 0; 
    handles.vars.params.selected = '';
    set(handles.edit_wingstat,'Enable','off');
else
    handles.vars.params.selected = 'stat_wing';
    set(handles.edit_wingstat,'Enable',onoff);
end
% Enable/Disable corresponding listbox for further selections & range editbox
handles = setplots(onoff, hObject, handles);
handles = savedat(hObject, handles);

function stat_circling_Callback(hObject, eventdata, handles)
% CHECKBOX TO PLOT CIRCLING STATISTICS
handles.vars.params.plots.stat.circling = get(handles.stat_circling,'Value');
if handles.vars.params.plots.stat.circling, onoff = 'on'; else onoff = 'off'; end
handles = set_selplots_off(hObject, handles);
if handles.vars.params.plots.movieclips,
    % In case the 'movieclips' checkbox is activated allow only this
    % action to be selected
    if handles.vars.params.plots.stat.circling, val = 1; else val = 0; end 
    handles.vars.params.plots.stat.chasing = 0;
    handles.vars.params.plots.stat.lunging = 0;
    handles.vars.params.plots.stat.wingthreat = 0;
    handles.vars.params.plots.stat.onewing = 0; 
    handles.vars.params.plots.stat.circling = val;
    handles.vars.params.plots.stat.copulation = 0;
    handles.vars.params.plots.stat.tussl = 0; 
    set(handles.edit_circlstat,'Enable','off');
else
    set(handles.edit_circlstat,'Enable',onoff);
end
handles = savedat(hObject, handles);

function stat_wingextflydist_Callback(hObject, eventdata, handles)
% CHECKBOX TO PLOT FLY DISTANCE WHILE WING EXTENSION STATISTICS
handles.vars.params.plots.stat.wingextflydist = get(handles.stat_wingextflydist,'Value');
handles = set_selplots_off(hObject, handles);
set(handles.edit_wingstat,'Enable','off');
handles = savedat(hObject, handles);

function stat_copulation_Callback(hObject, eventdata, handles)
% CHECKBOX TO PLOT COPULATION STATISTICS
handles.vars.params.plots.stat.copulation = get(handles.stat_copulation,'Value');
if handles.vars.params.plots.stat.copulation, onoff = 'on'; else onoff = 'off'; end
if handles.vars.params.plots.movieclips,
    % In case the 'movieclips' checkbox is activated allow only this
    % action to be selected
    if handles.vars.params.plots.stat.copulation, val = 1; else val = 0; end 
    handles.vars.params.plots.stat.chasing = 0;
    handles.vars.params.plots.stat.lunging = 0;
    handles.vars.params.plots.stat.wingthreat = 0;
    handles.vars.params.plots.stat.onewing = 0; 
    handles.vars.params.plots.stat.circling = 0;
    handles.vars.params.plots.stat.copulation = val;
    handles.vars.params.plots.stat.tussl = 0; 
    handles.vars.params.selected = '';
else
    handles.vars.params.selected = 'stat_copul';
end
% Enable/Disable corresponding listbox for further selections & range editbox
handles = setplots(onoff, hObject, handles);
handles = savedat(hObject, handles);

function stat_selAll_Callback(hObject, eventdata, handles)
% CHECKBOX TO SELECT/UNSELECT ALL STATISTIC PLOTS
if get(handles.stat_selAll,'Value'), val = 1; onoff = 'on'; else val = 0; onoff = 'off'; end
handles.vars.params.plots.stat.velocity = val;
handles.vars.params.plots.stat.chasing = val;
handles.vars.params.plots.stat.traveldist = val;
handles.vars.params.plots.stat.dist = val;
handles.vars.params.plots.stat.movestop = val;
handles.vars.params.plots.stat.body = val;
handles.vars.params.plots.stat.tussl = val;
handles.vars.params.plots.stat.lunging = val;
handles.vars.params.plots.stat.wingthreat = val;
handles.vars.params.plots.stat.scat = val;
handles.vars.params.plots.stat.leftwing = val;
handles.vars.params.plots.stat.rightwing = val;
handles.vars.params.plots.stat.onewing = val;
handles.vars.params.plots.stat.circling = val;
handles.vars.params.plots.stat.wingextflydist = val;
handles.vars.params.plots.stat.copulation = val;
handles.vars.params.plots.stat.actionrel = val;

if val, dummyvec = 1:20; else dummyvec = 0; end
fieldn = fieldnames(handles.vars.params.plots.stat);
for i=1:length(fieldn),
    if numel(strfind(fieldn{i},'vec')),
        handles.vars.params.plots.stat.(fieldn{i}) = dummyvec;
    end
end
set(handles.edit_chasestat,'Enable',onoff); set(handles.edit_walkstat,'Enable',onoff);
set(handles.edit_flydiststat,'Enable',onoff); set(handles.edit_velstat,'Enable',onoff);
set(handles.edit_tusslstat,'Enable',onoff); set(handles.edit_lungestat,'Enable',onoff);
set(handles.edit_wingthreatstat,'Enable',onoff); set(handles.edit_wingstat,'Enable',onoff);
set(handles.edit_wingstat,'Enable',onoff); set(handles.edit_wingstat,'Enable',onoff);
set(handles.edit_circlstat,'Enable',onoff); set(handles.edit_wingstat,'Enable','off');

% set(handles.stat_velocity,'Value',onoff); set(handles.stat_chasing,'Value',onoff);
% set(handles.stat_walkingdist,'Value',onoff); set(handles.stat_dist,'Value',onoff);
% set(handles.stat_movestop,'Value',onoff); set(handles.stat_bodysize,'Value',onoff);
% set(handles.stat_dendro,'Value',onoff); set(handles.stat_tussl,'Value',onoff);
% set(handles.stat_lunging,'Value',onoff); set(handles.stat_scat,'Value',onoff);
% set(handles.stat_leftwing,'Value',onoff); set(handles.stat_rightwing,'Value',onoff);
% set(handles.stat_onewing,'Value',onoff); set(handles.stat_circling,'Value',onoff);
% set(handles.stat_wingextflydist,'Value',onoff); set(handles.stat_copulation,'Value',onoff);
handles = savedat(hObject, handles);

function handles = enable_heatm_stat(onoff, hObject, handles)
% CHECK/UNCHECK ALL HEATMAP PLOT CHECKBOXES
set(handles.heatmap_position,'Enable',onoff); set(handles.heatmap_velocity,'Enable',onoff);
set(handles.heatmap_ext_leftwing,'Enable',onoff); set(handles.heatmap_ext_rightwing,'Enable',onoff);
set(handles.heatmap_ext_diffleftrightwing,'Enable',onoff); set(handles.heatmap_ext_onewing,'Enable',onoff);
set(handles.heatmap_chasing,'Enable',onoff); set(handles.heatmap_circling,'Enable',onoff);
set(handles.heatmap_tussling,'Enable',onoff); set(handles.heatmap_lunging,'Enable',onoff);
set(handles.heatmap_wingthreat,'Enable',onoff); set(handles.heatmap_wingrel,'Enable',onoff);
set(handles.heatmap_selAll,'Enable',onoff);

onoff2 = 'off';
set(handles.edit_posheat,'Enable',onoff2); set(handles.edit_velheat,'Enable',onoff2);
set(handles.edit_wingheat,'Enable',onoff2); set(handles.edit_wingheat,'Enable',onoff2);
set(handles.edit_wingheat,'Enable',onoff2); set(handles.edit_wingheat,'Enable',onoff2);
set(handles.edit_wingheat,'Enable',onoff2); set(handles.edit_chaseheat,'Enable',onoff2);
set(handles.edit_circlheat,'Enable',onoff2); set(handles.edit_lungeheat,'Enable',onoff2);
set(handles.edit_wingthreatheat,'Enable',onoff2); set(handles.edit_tusslheat,'Enable',onoff2);

set(handles.edit_chasestat,'Enable',onoff2); set(handles.edit_walkstat,'Enable',onoff2);
set(handles.edit_flydiststat,'Enable',onoff2); set(handles.edit_velstat,'Enable',onoff2);
set(handles.edit_tusslstat,'Enable',onoff2); set(handles.edit_lungestat,'Enable',onoff2);
set(handles.edit_wingthreatstat,'Enable',onoff2); set(handles.edit_wingstat,'Enable',onoff2);
set(handles.edit_wingstat,'Enable',onoff2); set(handles.edit_wingstat,'Enable',onoff2);
set(handles.edit_circlstat,'Enable',onoff2); set(handles.edit_wingstat,'Enable','off');

set(handles.stat_velocity,'Enable',onoff); set(handles.stat_chasing,'Enable',onoff);
set(handles.stat_walkingdist,'Enable',onoff); set(handles.stat_dist,'Enable',onoff);
set(handles.stat_movestop,'Enable',onoff); set(handles.stat_bodysize,'Enable',onoff);
set(handles.stat_dendro,'Enable',onoff); set(handles.stat_tussl,'Enable',onoff);
set(handles.stat_actionrel,'Enable',onoff);
set(handles.timeseries,'Enable',onoff);
set(handles.stat_lunging,'Enable',onoff); set(handles.stat_wingthreat,'Enable',onoff); set(handles.stat_scat,'Enable',onoff);
set(handles.stat_leftwing,'Enable',onoff); set(handles.stat_rightwing,'Enable',onoff);
set(handles.stat_onewing,'Enable',onoff); set(handles.stat_circling,'Enable',onoff);
set(handles.stat_wingextflydist,'Enable',onoff); set(handles.stat_copulation,'Enable',onoff);
set(handles.stat_selAll,'Enable',onoff);

if strcmp(onoff,'off'),
    val = 0;
    handles.vars.params.plots.heatmaps.position = val;
    handles.vars.params.plots.heatmaps.velocity = val;
    handles.vars.params.plots.heatmaps.ext_leftwing = val;
    handles.vars.params.plots.heatmaps.ext_rightwing = val;
    handles.vars.params.plots.heatmaps.ext_diffleftrightwing = val;
    handles.vars.params.plots.heatmaps.ext_onewing = val;
    handles.vars.params.plots.heatmaps.chasing = val;
    handles.vars.params.plots.heatmaps.circling = val;
    handles.vars.params.plots.heatmaps.tussling = val;
    handles.vars.params.plots.heatmaps.lunging = val;
    handles.vars.params.plots.heatmaps.wingrel = val;          
    
    handles.vars.params.plots.stat.velocity = val;
    handles.vars.params.plots.stat.chasing = val;
    handles.vars.params.plots.stat.traveldist = val;
    handles.vars.params.plots.stat.dist = val;
    handles.vars.params.plots.stat.movestop = val;
    handles.vars.params.plots.stat.body = val;
    handles.vars.params.plots.stat.dendro = val;
    handles.vars.params.plots.stat.actionrel = val;
    handles.vars.params.plots.stat.tussl = val;
    handles.vars.params.plots.stat.lunging = val;
    handles.vars.params.plots.stat.wingthreat = val;
    handles.vars.params.plots.stat.scat = val;
    handles.vars.params.plots.stat.leftwing = val;
    handles.vars.params.plots.stat.rightwing = val;
    handles.vars.params.plots.stat.onewing = val;
    handles.vars.params.plots.stat.circling = val;
    handles.vars.params.plots.stat.wingextflydist = val;
    handles.vars.params.plots.stat.copulation = val;

    handles.vars.params.plots.timeseries = val;
    handles.vars.params.plots.ethogram = val;    
end
handles = enable_ranges(onoff, hObject, handles);
handles = savedat(hObject, handles);


function enable_all(onoff, hObject, handles)
% ENABLE/DISABLE THE GUI FUNCTIONS
handles = enable_heatm_stat(onoff, hObject, handles);

set(handles.reset_button,'Enable',onoff);
set(handles.cumul_time,'Enable',onoff);

set(handles.edit_genotype,'Enable',onoff); set(handles.edit_save,'Enable',onoff);
set(handles.pushbutton_saveas,'Enable',onoff);

set(handles.edit_timelimit,'Enable',onoff); set(handles.edit_rows,'Enable',onoff);
set(handles.edit_cols,'Enable',onoff); set(handles.edit_chamberwidth,'Enable',onoff);
set(handles.edit_chamberheight,'Enable',onoff); set(handles.edit_axisfontsize,'Enable',onoff);
set(handles.edit_border,'Enable',onoff);
set(handles.checkbox_pdf,'Enable',onoff); 
set(handles.checkbox_text,'Enable',onoff); 
set(handles.checkbox_boxplot,'Enable',onoff);
set(handles.checkbox_meansem,'Enable',onoff);
set(handles.checkbox_bonf,'Enable',onoff);
% set(handles.checkbox_readnew,'Enable',onoff);
set(handles.checkbox_analyzenew,'Enable',onoff); 
set(handles.checkbox_winlos,'Enable',onoff); 
set(handles.checkbox_border,'Enable',onoff);
set(handles.checkbox_courtship,'Enable',onoff);
set(handles.checkbox_oneobj,'Enable',onoff);
set(handles.checkbox_autoscale,'Enable',onoff);
set(handles.pushbutton_tuning,'Enable',onoff);
set(handles.checkbox_circular,'Enable',onoff);
set(handles.clean,'Enable',onoff),
    
set(handles.listbox_genotypes,'Enable',onoff);
set(handles.listbox_files,'Enable',onoff);
set(handles.listbox_chambers,'Enable',onoff);
set(handles.listbox_n,'Enable',onoff);
set(handles.edit_genotype,'Enable',onoff);
set(handles.slider,'Enable',onoff);

set(handles.checkbox_genfolder,'Enable',onoff); set(handles.removefile,'Enable',onoff); 
set(handles.auto,'Enable',onoff); set(handles.remove_gen,'Enable',onoff); 

set(handles.analyze,'Enable',onoff); 
set(handles.ethogram,'Enable',onoff); 

%%% MOVIE CLIPS ==>  for enabling set 'off' to onoff
set(handles.movieclips,'Enable',onoff); 
set(handles.showhtml,'Enable',onoff);

set(handles.text6,'Enable',onoff); set(handles.text7,'Enable',onoff);
set(handles.text8,'Enable',onoff); set(handles.text9,'Enable',onoff);
set(handles.text_mm_arena,'Enable',onoff); set(handles.text11,'Enable',onoff);
set(handles.text_mm_border,'Enable',onoff); set(handles.text_border,'Enable',onoff);

if strcmp(onoff,'off'),
    set(handles.heatmap_selAll,'Value',0);
    set(handles.stat_selAll,'Value',0);
    set(handles.checkbox_circular,'Value',0);
end

handles = enable_ranges('off', hObject, handles);


function ethogram_Callback(hObject, eventdata, handles)
% CHECKBOX TO PLOT ETHOGRAMS
handles.vars.params.plots.ethogram = get(handles.ethogram,'Value');
handles.vars.params.selected = 'stat_ethogram';
% Enable/Disable corresponding listbox for further selections & range editbox
if handles.vars.params.plots.ethogram, onoff = 'on'; else onoff = 'off'; end
handles = setplots(onoff, hObject, handles);
handles = savedat(hObject, handles);

function movieclips_Callback(hObject, eventdata, handles)
% CHECKBOX FOR MOVIE CLIP EXTRACTION
handles.vars.params.plots.movieclips = get(handles.movieclips,'Value');
if handles.vars.params.plots.movieclips,
    handles.vars.params.plots.showhtml = 1;
    handles.vars.params.selected = '';
    handles = setplots('off', hObject, handles);
    handles = enable_heatm_stat('off', hObject, handles);
    onoff = 'on';
    set(handles.stat_chasing,'Enable',onoff);
    set(handles.stat_tussl,'Enable',onoff); set(handles.stat_lunging,'Enable',onoff);
    set(handles.stat_onewing,'Enable',onoff); set(handles.stat_circling,'Enable',onoff);
    set(handles.stat_wingthreat,'Enable',onoff); set(handles.stat_copulation,'Enable',onoff);
    handles.vars.params.plots.stat.chasing = 0;
    handles.vars.params.plots.stat.lunging = 0;
    handles.vars.params.plots.stat.wingthreat = 0;
    handles.vars.params.plots.stat.onewing = 0; 
    handles.vars.params.plots.stat.circling = 0;
    handles.vars.params.plots.stat.copulation = 0;
    handles.vars.params.plots.stat.tussl = 0; 
else
    handles.vars.params.plots.showhtml = 0;
    handles = enable_heatm_stat('on', hObject, handles);
end
handles = savedat(hObject, handles);

function showhtml_Callback(hObject, eventdata, handles)
% CHECKBOX FOR WEBSITE GENERATION FROM ALL AVAILABLE
% EXTRACTED MOVIE CLIPS
handles.vars.params.plots.showhtml = get(handles.showhtml,'Value');
if handles.vars.params.plots.showhtml,
    if ~handles.vars.params.plots.movieclips,
        handles.vars.params.plots.showhtml = 0;
        write_html(handles.vars.SavePathName);
    end
end
handles = savedat(hObject, handles);

function reset_button_Callback(hObject, eventdata, handles)
% RESET_BUTTON GUI
handles.reset = 1; guidata(hObject, handles);
for i=1:length(handles.vars.file),
    set(handles.listbox_genotypes,'Value',1);
    handles = removefile_Callback(hObject, eventdata, handles);
end
handles.vars.chamber.width = handles.vars.chamber.width0;
handles.vars.chamber.height = handles.vars.chamber.height0;
set(handles.edit_chamberwidth,'String',num2str(handles.vars.chamber.width,'%3.0f'));
set(handles.edit_chamberheight,'String',num2str(handles.vars.chamber.height,'%3.0f'));
handles.vars.params.border_width = handles.vars.params.border_width0;
set(handles.edit_border,'String',num2str(handles.vars.params.border_width,'%4.1f'));
handles.vars.SavePathName = []; set(handles.edit_save,'String',[]);
heatmap_selAll_Callback(hObject, eventdata, handles);
stat_selAll_Callback(hObject, eventdata, handles);
enable_all('off', hObject, handles);
set_selplots_off(hObject, handles);
handles.reset = 0; guidata(hObject, handles);

function edit_chaseheat_Callback(hObject, eventdata, handles)
% EDITBOX FOR UPPER BOUNDARY OF CHASING HEATMAP
tmp = abs(str2double(get(handles.edit_chaseheat,'String')));
if handles.vars.params.tim,
    handles.vars.params.range.heatmaps.tim.chase = tmp;
else
    handles.vars.params.range.heatmaps.occ.chase = tmp;
end
handles = savedat(hObject, handles);

function edit_chaseheat_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_posheat_Callback(hObject, eventdata, handles)
% EDITBOX FOR UPPER BOUNDARY OF POSITION HEATMAP
tmp = abs(str2double(get(handles.edit_posheat,'String')));
% handles.vars.params.selected = 'heatmap_position';
% handles = setplots('on', hObject, handles);
if handles.vars.params.tim,
    handles.vars.params.range.heatmaps.tim.pos = tmp;
else
    handles.vars.params.range.heatmaps.occ.pos = tmp;
end
handles = savedat(hObject, handles);

function edit_posheat_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_velheat_Callback(hObject, eventdata, handles)
% EDITBOX FOR UPPER BOUNDARY OF VELOCITY HEATMAP
tmp = (abs(str2double(get(handles.edit_velheat,'String'))));
if handles.vars.params.tim,
    handles.vars.params.range.heatmaps.tim.vel = tmp;
else
    handles.vars.params.range.heatmaps.occ.vel = tmp;
end
handles = savedat(hObject, handles);

function edit_velheat_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_wingheat_Callback(hObject, eventdata, handles)
% EDITBOX FOR UPPER BOUNDARY OF WING EXTENSION HEATMAP
tmp = (abs(str2double(get(handles.edit_wingheat,'String'))));
if handles.vars.params.tim,
    handles.vars.params.range.heatmaps.tim.wing = tmp;
else
    handles.vars.params.range.heatmaps.occ.wing = tmp;
end
handles = savedat(hObject, handles);
function edit_wingheat_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_circlheat_Callback(hObject, eventdata, handles)
% EDITBOX FOR UPPER BOUNDARY OF CIRCLING HEATMAP
tmp = (abs(str2double(get(handles.edit_circlheat,'String'))));
if handles.vars.params.tim,
    handles.vars.params.range.heatmaps.tim.circl = tmp;
else
    handles.vars.params.range.heatmaps.occ.circl = tmp;
end
handles = savedat(hObject, handles);
function edit_circlheat_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_lungeheat_Callback(hObject, eventdata, handles)
% EDITBOX FOR UPPER BOUNDARY OF LUNGING HEATMAP
tmp = (abs(str2double(get(handles.edit_lungeheat,'String'))));
if handles.vars.params.tim,
    handles.vars.params.range.heatmaps.tim.lunge = tmp;
else
    handles.vars.params.range.heatmaps.occ.lunge = tmp;
end
handles = savedat(hObject, handles);
function edit_lungeheat_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_wingthreatheat_Callback(hObject, eventdata, handles)
% EDITBOX FOR UPPER BOUNDARY OF WING THREAT HEATMAP
tmp = (abs(str2double(get(handles.edit_wingthreatheat,'String'))));
if handles.vars.params.tim,
    handles.vars.params.range.heatmaps.tim.wingthreat = tmp;
else
    handles.vars.params.range.heatmaps.occ.wingthreat = tmp;
end
handles = savedat(hObject, handles);
function edit_wingthreatheat_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_tusslheat_Callback(hObject, eventdata, handles)
% EDITBOX FOR UPPER BOUNDARY OF TUSSLING HEATMAP
tmp = (abs(str2double(get(handles.edit_tusslheat,'String'))));
if handles.vars.params.tim,
    handles.vars.params.range.heatmaps.tim.tussl = tmp;
else
    handles.vars.params.range.heatmaps.occ.tussl = tmp;
end
handles = savedat(hObject, handles);
function edit_tusslheat_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_chasestat_Callback(hObject, eventdata, handles)
% EDITBOX FOR UPPER BOUNDARY OF CHASING STATISTICS
tmp = (abs(str2double(get(handles.edit_chasestat,'String'))));
if handles.vars.params.tim,
    handles.vars.params.range.stat.tim.chase = tmp;
else
    handles.vars.params.range.stat.occ.chase = tmp;
end
handles = savedat(hObject, handles);
function edit_chasestat_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_circlstat_Callback(hObject, eventdata, handles)
% EDITBOX FOR UPPER BOUNDARY OF CIRCLING STATISTICS
tmp = (abs(str2double(get(handles.edit_circlstat,'String'))));
if handles.vars.params.tim,
    handles.vars.params.range.stat.tim.circl = tmp;
else
    handles.vars.params.range.stat.occ.circl = tmp;
end
handles = savedat(hObject, handles);
function edit_circlstat_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_lungestat_Callback(hObject, eventdata, handles)
% EDITBOX FOR UPPER BOUNDARY OF LUNGING STATISTICS
tmp = (abs(str2double(get(handles.edit_lungestat,'String'))));
if handles.vars.params.tim,
    handles.vars.params.range.stat.tim.lunge = tmp;
else
    handles.vars.params.range.stat.occ.lunge = tmp;
end
handles = savedat(hObject, handles);
function edit_lungestat_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_wingthreatstat_Callback(hObject, eventdata, handles)
% EDITBOX FOR UPPER BOUNDARY OF WING THREAT STATISTICS
tmp = (abs(str2double(get(handles.edit_wingthreatstat,'String'))));
if handles.vars.params.tim,
    handles.vars.params.range.stat.tim.wingthreat = tmp;
else
    handles.vars.params.range.stat.occ.wingthreat = tmp;
end
handles = savedat(hObject, handles);
function edit_wingthreatstat_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_tusslstat_Callback(hObject, eventdata, handles)
% EDITBOX FOR UPPER BOUNDARY OF TUSSLING STATISTICS
tmp = (abs(str2double(get(handles.edit_tusslstat,'String'))));
if handles.vars.params.tim,
    handles.vars.params.range.stat.tim.tussl = tmp;
else
    handles.vars.params.range.stat.occ.tussl = tmp;
end
handles = savedat(hObject, handles);
function edit_tusslstat_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_walkstat_Callback(hObject, eventdata, handles)
% EDITBOX FOR UPPER BOUNDARY OF WALKING STATISTICS
tmp = (abs(str2double(get(handles.edit_walkstat,'String'))));
% if handles.vars.params.tim,
    handles.vars.params.range.stat.tim.walk = tmp;
% else
    handles.vars.params.range.stat.occ.walk = tmp;
% end
handles = savedat(hObject, handles);
function edit_walkstat_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_flydiststat_Callback(hObject, eventdata, handles)
% EDITBOX FOR UPPER BOUNDARY OF FLY DISTANCE STATISTICS
tmp = (abs(str2double(get(handles.edit_flydiststat,'String'))));
handles.vars.params.range.stat.tim.flydist = tmp;
handles.vars.params.range.stat.occ.flydist = tmp;
handles = savedat(hObject, handles);
function edit_flydiststat_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_velstat_Callback(hObject, eventdata, handles)
% EDITBOX FOR UPPER BOUNDARY OF VELOCITY STATISTICS
tmp = (abs(str2double(get(handles.edit_velstat,'String'))));
if handles.vars.params.tim,
    handles.vars.params.range.stat.tim.vel = tmp;
else
    handles.vars.params.range.stat.occ.vel = tmp;
end
handles = savedat(hObject, handles);
function edit_velstat_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_wingstat_Callback(hObject, eventdata, handles)
% EDITBOX FOR UPPER BOUNDARY OF WING EXTENSION STATISTICS
tmp = (abs(str2double(get(handles.edit_wingstat,'String'))));
if handles.vars.params.tim,
    handles.vars.params.range.stat.tim.wing = tmp;
else
    handles.vars.params.range.stat.occ.wing = tmp;
end
handles = savedat(hObject, handles);
function edit_wingstat_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function handles = enable_ranges(onoff, hObject, handles)
% ENABLE/DISABLE RANGE EDITBOXES
if handles.vars.params.plots.heatmaps.chasing,
    set(handles.edit_chaseheat,'Enable',onoff);
end
if handles.vars.params.plots.heatmaps.position,
    set(handles.edit_posheat,'Enable',onoff);
end
if handles.vars.params.plots.heatmaps.velocity,
    set(handles.edit_velheat,'Enable',onoff);
end
if handles.vars.params.plots.heatmaps.ext_leftwing || ...
        handles.vars.params.plots.heatmaps.ext_rightwing || ...
        handles.vars.params.plots.heatmaps.ext_diffleftrightwing || ...
        handles.vars.params.plots.heatmaps.ext_onewing || ...
        handles.vars.params.plots.heatmaps.wingrel,
    set(handles.edit_wingheat,'Enable',onoff);
end
if handles.vars.params.plots.heatmaps.circling,
    set(handles.edit_circlheat,'Enable',onoff);
end
if handles.vars.params.plots.heatmaps.lunging,
    set(handles.edit_lungeheat,'Enable',onoff);
end
if handles.vars.params.plots.heatmaps.wingthreat,
    set(handles.edit_wingthreatheat,'Enable',onoff);
end
if handles.vars.params.plots.heatmaps.tussling,
    set(handles.edit_tusslheat,'Enable',onoff);
end

if handles.vars.params.plots.stat.traveldist,
    set(handles.edit_walkstat,'Enable',onoff);
end
if handles.vars.params.plots.stat.dist,
    set(handles.edit_flydiststat,'Enable',onoff);
end
if handles.vars.params.plots.stat.velocity,
    set(handles.edit_velstat,'Enable',onoff);
end
if handles.vars.params.plots.stat.leftwing || ...
        handles.vars.params.plots.stat.rightwing || ...
        handles.vars.params.plots.stat.onewing,
    set(handles.edit_wingstat,'Enable',onoff);
end
if handles.vars.params.plots.stat.chasing,
    set(handles.edit_chasestat,'Enable',onoff);
end
if handles.vars.params.plots.stat.circling,
    set(handles.edit_circlstat,'Enable',onoff);
end
if handles.vars.params.plots.stat.lunging,
    set(handles.edit_lungestat,'Enable',onoff);
end
if handles.vars.params.plots.stat.wingthreat,
    set(handles.edit_wingthreatstat,'Enable',onoff);
end
if handles.vars.params.plots.stat.tussl,
    set(handles.edit_tusslstat,'Enable',onoff);
end

% function edit_rangemax_Callback(hObject, eventdata, handles)
% tmp = (abs(str2double(get(handles.edit_rangemax,'String'))));
% if strcmp(handles.vars.params.selected,'heatmap_position'),
%     if handles.vars.params.tim,
%         handles.vars.params.range.heatmaps.tim.pos = tmp;
%     else
%         handles.vars.params.range.heatmaps.occ.pos = tmp;
%     end
% end
% handles = savedat(hObject, handles);
% 
% function edit_rangemax_CreateFcn(hObject, eventdata, handles)
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end

function listbox_selplots_Callback(hObject, eventdata, handles)
% WHICH PLOTS ARE SELECTED IN 'SELPLOTS' LISTBOX?
if strcmp(handles.vars.params.selected,'heatmap_position'),
    handles.vars.params.plots.heatmaps.position_vec = get(handles.listbox_selplots,'Value');
end
if strcmp(handles.vars.params.selected,'heatmap_velocity'),
    handles.vars.params.plots.heatmaps.velocity_vec = get(handles.listbox_selplots,'Value');
end
if strcmp(handles.vars.params.selected,'heatmap_wing'),
    handles.vars.params.plots.heatmaps.wing_vec = get(handles.listbox_selplots,'Value');
end
if strcmp(handles.vars.params.selected,'heatmap_wingrel'),
    handles.vars.params.plots.heatmaps.wingrel_vec = get(handles.listbox_selplots,'Value');
end
if strcmp(handles.vars.params.selected,'heatmap_chase'),
    handles.vars.params.plots.heatmaps.chase_vec = get(handles.listbox_selplots,'Value');
end
if strcmp(handles.vars.params.selected,'heatmap_circl'),
    handles.vars.params.plots.heatmaps.circl_vec = get(handles.listbox_selplots,'Value');
end
if strcmp(handles.vars.params.selected,'heatmap_lunge'),
    handles.vars.params.plots.heatmaps.lunge_vec = get(handles.listbox_selplots,'Value');
end
if strcmp(handles.vars.params.selected,'heatmap_wingthreat'),
    handles.vars.params.plots.heatmaps.wingthreat_vec = get(handles.listbox_selplots,'Value');
end
if strcmp(handles.vars.params.selected,'heatmap_tussl'),
    handles.vars.params.plots.heatmaps.tussl_vec = get(handles.listbox_selplots,'Value');
end
if strcmp(handles.vars.params.selected,'stat_flydist'),
    handles.vars.params.plots.stat.flydist_vec = get(handles.listbox_selplots,'Value');
end
if strcmp(handles.vars.params.selected,'stat_body'),
    handles.vars.params.plots.stat.body_vec = get(handles.listbox_selplots,'Value');
end
if strcmp(handles.vars.params.selected,'stat_vel'),
    handles.vars.params.plots.stat.vel_vec = get(handles.listbox_selplots,'Value');
end
if strcmp(handles.vars.params.selected,'stat_wing'),
    handles.vars.params.plots.stat.wing_vec = get(handles.listbox_selplots,'Value');
end
if strcmp(handles.vars.params.selected,'stat_lunge'),
    handles.vars.params.plots.stat.lunge_vec = get(handles.listbox_selplots,'Value');
end
if strcmp(handles.vars.params.selected,'stat_wingthreat'),
    handles.vars.params.plots.stat.wingthreat_vec = get(handles.listbox_selplots,'Value');
end
if strcmp(handles.vars.params.selected,'stat_tussl'),
    handles.vars.params.plots.stat.tussl_vec = get(handles.listbox_selplots,'Value');
end
if strcmp(handles.vars.params.selected,'stat_scat'),
    handles.vars.params.plots.stat.scat_vec = get(handles.listbox_selplots,'Value');
end
if strcmp(handles.vars.params.selected,'stat_copul'),
    handles.vars.params.plots.stat.copul_vec = get(handles.listbox_selplots,'Value');
end
if strcmp(handles.vars.params.selected,'stat_actionrel'),
    handles.vars.params.plots.stat.actionrel_vec = get(handles.listbox_selplots,'Value');
end
if strcmp(handles.vars.params.selected,'timeseries'),
    handles.vars.params.plots.timeseries_vec = get(handles.listbox_selplots,'Value');
end
if strcmp(handles.vars.params.selected,'stat_ethogram'),
    handles.vars.params.plots.ethogram_vec = get(handles.listbox_selplots,'Value');
end
handles = savedat(hObject, handles);

function listbox_selplots_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function handles = set_selplots_off(hObject, handles)
% DISABLE AND CLEAR 'SELPLOTS' LISTBOX
set(handles.listbox_selplots,'Enable','off'); 
set(handles.listbox_selplots,'String',' '); 
set(handles.listbox_selplots,'Value',1);
guidata(hObject, handles);


function handles = setplots(onoff, hObject, handles)
% SHOW AVAILABLE PLOTS IN 'SELPLOTS' LISTBOX
% if strcmp(onoff,'on'),
    if strcmp(handles.vars.params.selected,'heatmap_position'),
        str = {'Mean','SEM'};
        if ~handles.vars.params.plots.heatmaps.position_vec, 
            handles.vars.params.plots.heatmaps.position_vec = 1:length(str);
        end
        val = handles.vars.params.plots.heatmaps.position_vec;
    end
    if strcmp(handles.vars.params.selected,'heatmap_velocity'),
        str = {'Mean','Direction'};
        if ~handles.vars.params.plots.heatmaps.velocity_vec, 
            handles.vars.params.plots.heatmaps.velocity_vec = 1:length(str);
        end
        val = handles.vars.params.plots.heatmaps.velocity_vec;
    end
    if strcmp(handles.vars.params.selected,'heatmap_wing'),
        str = {'Frequency','Temporal Distr.'};
        if ~handles.vars.params.plots.heatmaps.wing_vec, 
            handles.vars.params.plots.heatmaps.wing_vec = 1:length(str);
        end
        val = handles.vars.params.plots.heatmaps.wing_vec;
    end
    if strcmp(handles.vars.params.selected,'heatmap_wingrel'),
        str = {'Left & Right Wing','Left-Right Wing Difference',...
            'Circling Direction',...
            'Head Positions','Head Position Time Series',...
            'Wing Threats','Wing Threat Time Series'};
        if ~handles.vars.params.plots.heatmaps.wingrel_vec, 
            handles.vars.params.plots.heatmaps.wingrel_vec = 1:length(str);
        end
        val = handles.vars.params.plots.heatmaps.wingrel_vec;
    end
    if strcmp(handles.vars.params.selected,'heatmap_chase'),
        str = {'Frequency','Temporal Distr.'};
        if ~handles.vars.params.plots.heatmaps.chase_vec, 
            handles.vars.params.plots.heatmaps.chase_vec = 1:length(str);
        end
        val = handles.vars.params.plots.heatmaps.chase_vec;
    end
    if strcmp(handles.vars.params.selected,'heatmap_circl'),
        str = {'Frequency','Temporal Distr.'};
        if ~handles.vars.params.plots.heatmaps.circl_vec, 
            handles.vars.params.plots.heatmaps.circl_vec = 1:length(str);
        end
        val = handles.vars.params.plots.heatmaps.circl_vec;
    end
    if strcmp(handles.vars.params.selected,'heatmap_lunge'),
        str = {'Frequency','Temporal Distr.'};
        if ~handles.vars.params.plots.heatmaps.lunge_vec, 
            handles.vars.params.plots.heatmaps.lunge_vec = 1:length(str);
        end
        val = handles.vars.params.plots.heatmaps.lunge_vec;
    end
    if strcmp(handles.vars.params.selected,'heatmap_wingthreat'),
        str = {'Frequency','Temporal Distr.'};
        if ~handles.vars.params.plots.heatmaps.wingthreat_vec, 
            handles.vars.params.plots.heatmaps.wingthreat_vec = 1:length(str);
        end
        val = handles.vars.params.plots.heatmaps.wingthreat_vec;
    end
    if strcmp(handles.vars.params.selected,'heatmap_tussl'),
        str = {'Frequency','Temporal Distr.'};
        if ~handles.vars.params.plots.heatmaps.tussl_vec, 
            handles.vars.params.plots.heatmaps.tussl_vec = 1:length(str);
        end
        val = handles.vars.params.plots.heatmaps.tussl_vec;
    end
    if strcmp(handles.vars.params.selected,'stat_flydist'),
        str = {'Distance from Center','Relative Fly Distance',...
               'Fly Proximity Period','Time in Diff. Distances'};
        if ~handles.vars.params.plots.stat.flydist_vec, 
            handles.vars.params.plots.stat.flydist_vec = 1:length(str);
        end
        val = handles.vars.params.plots.stat.flydist_vec;
    end
    if strcmp(handles.vars.params.selected,'stat_body'),
        str = {'Fly Length','Fly Area','Fly Length Difference','Fly Area Difference'};
        if ~handles.vars.params.plots.stat.body_vec, 
            handles.vars.params.plots.stat.body_vec = 1:length(str);
        end
        val = handles.vars.params.plots.stat.body_vec;
    end
    if strcmp(handles.vars.params.selected,'stat_vel'),
        str = {'Vel. under Diff. Fly Proximities','Total Duration of Velocities','Approach Velocities','Escape Velocities'};
        if ~handles.vars.params.plots.stat.vel_vec, 
            handles.vars.params.plots.stat.vel_vec = 1:length(str);
        end
        val = handles.vars.params.plots.stat.vel_vec;
    end
    if strcmp(handles.vars.params.selected,'stat_wing'),
        str = {'Frequency','Duration'};
        if ~handles.vars.params.plots.stat.wing_vec, 
            handles.vars.params.plots.stat.wing_vec = 1:length(str);
        end
        val = handles.vars.params.plots.stat.wing_vec;
    end
    if strcmp(handles.vars.params.selected,'stat_lunge'),
        str = {'Min. Lat. 1st Lunge','Cumulative Lunge Number','Latency 1st Lunge','Fly Pairs > 1 Lunge','Fly Pairs over Lunge Number','Lunge Frequency over Time','Frequency','Lunges per Meter','Lunges per Minute'};
        if ~handles.vars.params.plots.stat.lunge_vec, 
            handles.vars.params.plots.stat.lunge_vec = 1:length(str);
        end
        val = handles.vars.params.plots.stat.lunge_vec;
    end
    if strcmp(handles.vars.params.selected,'stat_tussl'),
        str = {'Min. Lat. 1st Tussling','Latency 1st Tussling','Mean Cumulative Tussling','Fly Pais > 1 Tussling','Tussling Bout Period over Time','Tussling Frequency over Time','Frequency','Tussling per Meter','Tussling per Minute'};
        if ~handles.vars.params.plots.stat.tussl_vec, 
            handles.vars.params.plots.stat.tussl_vec = 1:length(str);
        end
        val = handles.vars.params.plots.stat.tussl_vec;
    end
    if strcmp(handles.vars.params.selected,'stat_scat'),
        str = {'Lunges <-> Tussling','Wing Thr. <-> Lunging','Wing Thr. <-> Tussling','Chasing <-> Lunging','Chasing <-> Tussling','Wing Ext. <-> Lunging'};
        if ~handles.vars.params.plots.stat.scat_vec, 
            handles.vars.params.plots.stat.scat_vec = 1:length(str);
        end
        val = handles.vars.params.plots.stat.scat_vec;
    end
    if strcmp(handles.vars.params.selected,'stat_copul'),
        str = {'Begin Copulation','Copulation Period','End Copulation','Latency to 1st Wing','Copulation Index'};
        if ~handles.vars.params.plots.stat.copul_vec, 
            handles.vars.params.plots.stat.copul_vec = 1:length(str);
        end
        val = handles.vars.params.plots.stat.copul_vec;
    end
    if strcmp(handles.vars.params.selected,'stat_actionrel'),
        str = {'Incl. Chasing','Percentage','Total','Percentage Time Series'};
        if ~handles.vars.params.plots.stat.actionrel_vec, 
            handles.vars.params.plots.stat.actionrel_vec = 1:length(str);
        end
        val = handles.vars.params.plots.stat.actionrel_vec;
    end
    if strcmp(handles.vars.params.selected,'timeseries'),
        str = {'Velocities'};
        if ~handles.vars.params.plots.timeseries_vec, 
            handles.vars.params.plots.timeseries_vec = 1:length(str);
        end
        val = handles.vars.params.plots.timeseries_vec;
    end
    if strcmp(handles.vars.params.selected,'stat_ethogram'),
        str = {'Transitions','Reactions','Trans. Time Series','React. Time Series',...
               'Trans. Matrices','React. Matrices'};
        if ~handles.vars.params.plots.ethogram_vec, 
            handles.vars.params.plots.ethogram_vec = 1:length(str);
        end
        val = handles.vars.params.plots.ethogram_vec;
    end
    if isempty(handles.vars.params.selected) || (numel(str) < max(val)),
        str = ' '; val = 1; onoff = 'off';
    end
    set(handles.listbox_selplots,'String',str);
    set(handles.listbox_selplots,'Value',val);
% else
%     set(handles.listbox_selplots,'String',' ');
%     set(handles.listbox_selplots,'Value',1);
% end
set(handles.listbox_selplots,'Enable',onoff);
guidata(hObject, handles);

function pushbutton13_Callback(hObject, eventdata, handles)

function pushbutton14_Callback(hObject, eventdata, handles)

function pushbutton_tuning_Callback(hObject, eventdata, handles)
% CALL GUI FOR PARAMETER TUNING
old_tune = handles.vars.params.tune;

handles.vars.params.tune = tuning(handles.vars.params.tune);

if handles.vars.params.tune.reanalyze,
    handles.vars.params.analyze_new = 1;
    handles.vars.params.feat_read_new = 1;
end
set(handles.checkbox_analyzenew,'Value',handles.vars.params.analyze_new);
handles = savedat(hObject, handles);


function checkbox_circular_Callback(hObject, eventdata, handles)
% CALL GUI FOR CIRCULAR AREA OF INTEREST SELECTION
handles.vars.params.circular = get(handles.checkbox_circular,'Value');
old_radius = handles.vars.params.radius;
if handles.vars.params.circular,
    handles.vars.params.border = 0;
    
    [handles.vars.params.radius handles.vars.params.radius0] = ...
        border(handles.vars);
    
    handles.vars.params.border_width = handles.vars.params.radius0 - ...
        handles.vars.params.radius;
else
    handles.vars.params.border = 1;
    handles.vars.params.radius = -99;
end
if old_radius ~= handles.vars.params.radius && handles.vars.params.radius > 0,
    handles.vars.params.analyze_new = 1;
    handles.vars.params.feat_read_new = 1;
end
set(handles.checkbox_analyzenew,'Value',handles.vars.params.analyze_new);
set(handles.checkbox_border,'Value',handles.vars.params.border);
set(handles.edit_border,'String',num2str(handles.vars.params.border_width));
handles.vars.chamber.height = handles.vars.params.radius0*2;
handles.vars.chamber.width = handles.vars.params.radius0*2;
handles = savedat(hObject, handles);


function edit_border_Callback(hObject, eventdata, handles)
handles.vars.params.border_width = abs(str2double(get(handles.edit_border,'String')));
handles.vars.params.border_width = ...
    min(handles.vars.params.border_width,...
    min(handles.vars.chamber.width,handles.vars.chamber.height)/2-2.5);
if handles.vars.params.radius > 0,
    handles.vars.params.radius = handles.vars.params.radius0 - ...
        handles.vars.params.border_width;
end
handles = savedat(hObject, handles);

function edit_border_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
