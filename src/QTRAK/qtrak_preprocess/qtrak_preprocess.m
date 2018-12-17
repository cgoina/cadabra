% qtrak_preprocess : preprocess movie files. 
%
% Usage:  qtrak_prepocess
%
% qtrak_preprocess is derived from qtrak.m. In order to be adapted for the 
% cluster computing, qtrak.m was modifed to two parts: qtrak_preprocess and 
% qtrak_cluster. qtrak_preprocess inludes all the GUI for user to input 
% experimental parameters. It generates the ROI mat file, errlog file, empty 
% feature files, and a configure file, preprocessed.cofnig. These files are 
% necessary for the execustion of the second part, qtrak_cluster.
%
% * Original implementation by Heiko Dankert, modified by Jinyang Liu
% * Optimization and documentation by Edwin Soedarmadji

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
% * Original implementation by Heiko Dankert, modified by Jinyang Liu
% * Optimization and documentation by Edwin Soedarmadji
%
%%
% This is the main entry point of the entire tracking application. Most 
% actions performed in this module is self-explanatory, and are explained 
% in inline comments. The function makes heavy use of global variables, a 
% feature that can be revised in the next iteration of this code. 

function qtrak_preprocess

%#function ProcessFrameCapWin
%#function PlotImage_Captioning
%#function histfit 
%#function fhist gmmactiv gmmpost
%#function FFGrab seg_fullfly_inflexbox
%#function fhistc ellipsefit write_sequence
%#function gauss write_sequence_crop callback
%#function seg_fullfly_inconstbox segfly fotsu
%#function intersect_flies fly_bodymeas switch_flies
%#function fly_headtailwings mask2ellipses gmm gmminit consist
%#function kmeans dist2 gmmem gmmprob histfit
%#function auto_cal auto_cal_circ CircularHough_Grd fhistc openavi
%#function calibrate calibrate fhistfit GNU_message

global params chamber scale;
global slashstr Files NFiles strInVideoFNameArray;
global FigureHandle Panels PanelSize;
global DimX DimY DimZ Spacer;
global intStartFrm intNFrms dt;
global nframes_mean mean_image std_image frmindx;

    intStartFrm = 30;                           % The first frame processed
    nframes_mean = 4000;                        % The number of frames used to calculate the mean image

    params.cs = 8;                              % Border area (pixel) excluded from our analysis
    params.dist_mm = 10.5;                      % Calibrated distance between two known points  for Heisenberg Chamber
    params.distcirc_mm = 20;                    % Calibrated distance between two known points  for Heisenberg Chamber
    params.twopi = 2 * pi;                      % Define the useful constant 2\pi
    params.hpi = pi / 2;                        % and \frac{1}{2}\pi
    params.nchambers = 2;                       % The number of chambers in the frame 
    params.firstfile = 0;
    params.findex = 0;                          % FilterIndex return value of Matlab's uigetfile 
    params.filechg = 0;                         % Flag to indicate that user has selected a file
    params.QuitPrg = 0;                         % If set, user desires to quit through user interface
    params.ClosingEvent = 0;
    params.bool_plottrack = 0;                  % If set, shows tracks on live video 
    params.bool_plotcount = 1;                  % If set, updates frame count on the screen
    params.bool_boundbox = 1;                   % If set, shows the calibration bounds
    params.bool_asked = 0;                      % Ask questions? 0 = yes
    params.bool_cal = 0;                        % Fixed a bug. When choose recalibrate, the movie should be calibrated 
                                                % no matter roi.mat exists or not
    params.bool_recal = 1;                      % Recalibrate? 1 = yes
    params.chamberNum = 1;                      % User input the number of the chamber in the movie
    params.unprocessedMvNum = 0;                % When the number of chamber deteced by auto_ca;_circ is different from 
                                                % the number of chamber users input, the movie files failed the calibration 
                                                % and leave it for the user to calibrate later. params.unprocessedMvNum is 
                                                % a variable to count the number of unprocessed movie number.
    params.unprocessedMvName = cell(0);         % params.unprocessedMvName records the names of the unprocessed movies
    params.stopProcess = 0;                     % params.stopProcess records the state of the processing movie. if it is one,
                                                % the movie failed the calibraton.
    params.batchprocess = 1;
    
%     params.bool_dot = 0;                        % If set, one fly is marked with a dot
%     params.bool_court = 0;                      % Courtship events expected? 1 = yes
%     params.bool_nfly = 1;                       % More than 1 fly in a chamber? 1 = yes

    if (ispc)                                   % Is this program running on a PC?
        slashstr = '\';                         % If yes, then the directory separator is a '\'
    else
        slashstr = '/';                         % Otherwise, it is a forward slash '/' 
    end
    
    warning off all;                            % Surpress all warnings
    
    msg = GNU_message(1);                       % Print short GNU license message into command line window
    for i=1:numel(msg), fprintf([msg{i} '\n']); end
    
    if ~params.bool_asked

        %%
        % Build the user interface by creating a new figure.
        % Next, calculate the various dimensions for our user interface, 
        % and adjust the position of our figure by modifying its 
        % rectangular coordinate to the dimensions we previously calculated. 

        FigureHandle = figure;                      
        OldDimY = 576;                              
        OldDimX = 720;                              
        OldDimZ = 3;                                
        DimZ = OldDimZ;
        DimY = 480;
        DimX = round(OldDimX * DimY / OldDimY);
        Spacer =  0.02;
        Panels.PanelDimX = 160;

        rect = get(FigureHandle,'Position');        
        if (ispc)
            rect(1) = 100;
            rect(2) = 100;
        end
        rect(4) = DimY + (4 * DimY * Spacer);
        rect(3) = DimX + Panels.PanelDimX + (6 * DimY * Spacer);
        set(FigureHandle,'Position',rect);

        %%
        % Having positioned and sized our figure, calculate the dimensions 
        % for various panels to be placed within the figure, and store 
        % these dimensions in the PanelSize variable. Four panels are used:
        % Image, Info, Menu, and Events, each with its own flag that indicates
        % whether or not the panel needs to be refreshed. 

        NormVec = rect([3,4,3,4]);    
        PanelSize.Image.Left = (2 * DimY * Spacer);
        PanelSize.Image.Low = (2 * DimY * Spacer);
        PanelSize.Image.Height = DimY;
        PanelSize.Image.Width = DimX;
        PanelSize.Image.PVec = [PanelSize.Image.Left, PanelSize.Image.Low, ...
                                PanelSize.Image.Width, PanelSize.Image.Height]./NormVec;
        PanelSize.Info.Left = (4 * DimY * Spacer) + DimX;
        PanelSize.Info.Low = (6 * DimY * Spacer) + (DimY - 4 * DimY * Spacer) * 14 / 15;
        PanelSize.Info.Height = (DimY - 2 * DimY * Spacer) * 1 / 15;
        PanelSize.Info.Width = Panels.PanelDimX;
        PanelSize.Info.PVec = [PanelSize.Info.Left, PanelSize.Info.Low, ... 
                               PanelSize.Info.Width, PanelSize.Info.Height]./NormVec;                       
        PanelSize.Menu.Left = (4 * DimY * Spacer) + DimX;
        PanelSize.Menu.Low = (4 * DimY * Spacer) + (DimY - 4 * DimY * Spacer) * 5 / 15;
        PanelSize.Menu.Height = (DimY - 4 * DimY * Spacer) * 9 / 15;
        PanelSize.Menu.Width = Panels.PanelDimX;
        PanelSize.Menu.PVec = [PanelSize.Menu.Left, PanelSize.Menu.Low, ... 
                               PanelSize.Menu.Width, PanelSize.Menu.Height]./NormVec;
        PanelSize.Events.Left = (4 * DimY * Spacer) + DimX;
        PanelSize.Events.Low = (2 * DimY * Spacer);
        PanelSize.Events.Height = (DimY - 4 * DimY * Spacer) * 5 / 15;
        PanelSize.Events.Width = Panels.PanelDimX;
        PanelSize.Events.PVec = [PanelSize.Events.Left, PanelSize.Events.Low, ... 
                                 PanelSize.Events.Width, PanelSize.Events.Height]./NormVec;

        %% 
        % Next, we actually create arguments to Matlab functions that will
        % actually create the GUI controls, as well as creating some of the 
        % controls themselves (with the subplot function) and saving the handles.  

        Panels.EventsPanelDimY = PanelSize.Events.Height;
        Panels.MenuPanelDimY = PanelSize.Menu.Height;
        Panels.InfoPanelDimY = PanelSize.Info.Height;
        Panels.BlankEventsPanel = ones(round(Panels.EventsPanelDimY),Panels.PanelDimX,DimZ)/4;
        Panels.BlankMenuPanel = ones(round(Panels.MenuPanelDimY),Panels.PanelDimX,DimZ)/4;
        Panels.BlankInfoPanel = ones(round(Panels.InfoPanelDimY),Panels.PanelDimX,DimZ)/4;
        Panels.HImage = subplot('position',PanelSize.Image.PVec);
        Panels.HEvents = subplot('position',PanelSize.Events.PVec);
        Panels.HMenu = subplot('position',PanelSize.Menu.PVec);
        Panels.HInfo = subplot('position',PanelSize.Info.PVec);

        
        %%
        % Place the controls on FigureHandle, specifically asking
        % Matlab not to draw the tick marks. Also set other properties. 

        set(Panels.HImage,'Parent', FigureHandle,'XTick',[],'YTick',[]);
        axes(Panels.HMenu); image(Panels.BlankMenuPanel);
        set(Panels.HMenu, 'Parent', FigureHandle,'XTick',[],'YTick',[]);
        axes(Panels.HInfo); image(Panels.BlankInfoPanel);
        set(Panels.HInfo, 'Parent', FigureHandle,'XTick',[],'YTick',[]);
        axes(Panels.HEvents); image(Panels.BlankEventsPanel);
        set(Panels.HEvents, 'Parent', FigureHandle,'XTick',[],'YTick',[]);
        set(FigureHandle,'DoubleBuffer','on');
        set(FigureHandle,'Resize','off');
        set(FigureHandle,'Color','k');
        set(FigureHandle,'ToolBar','none');
        set(FigureHandle,'MenuBar','none');
        set(FigureHandle,'NumberTitle','off');

        %% 
        % Create the menu and submenu items along with specifying
        % a shared, generic callback function to all the menu items.

        figure(FigureHandle);
        f = uimenu('Label','File');
        g = uimenu('Label','View');
        uimenu('Label','About','enable','on','Callback','callback');
        uimenu(f,'Label','Open','enable','on','Callback','callback');
        uimenu(f,'Label','Quit [Esc]','enable','on','Callback','callback');
        uimenu(g,'Label','Show Tracking','enable','on','Checked', 'off','Callback','callback');
        uimenu(g,'Label','Show Counting','enable','on','Checked', 'on','Callback','callback');
        uimenu(g,'Label','Show Bounding Box','enable','on','Checked', 'on','Callback','callback');    

        % Print short GNU license message into command line window
        msg = GNU_message(1);
        sep1 = strfind(msg{1},'Heiko'); sep2 = strfind(msg{1},'Cali');
        text(PanelSize.Events.Width/2,PanelSize.Events.Height/2-15,msg{1}(1:sep1-1),'FontSize',7,'color',[.7 .7 .7],'HorizontalAlignment','Center');
        text(PanelSize.Events.Width/2,PanelSize.Events.Height/2,msg{1}(sep1:sep2-1),'FontSize',7,'color',[.7 .7 .7],'HorizontalAlignment','Center');
        text(PanelSize.Events.Width/2,PanelSize.Events.Height/2+15,msg{1}(sep2:end),'FontSize',7,'color',[.7 .7 .7],'HorizontalAlignment','Center');
    end

    %%
    % Open a dialog box for choosing the video file(s) to be processed. 
    % During benchmarking, openavi is executed in a non-interactive mode,
    % which is done by setting 0 in its input argument. 
    
    params.findex = 1;
    openavi(1);
    params.filechg = 0;        

    
    FileCount = 1;
    
    if NFiles <= 1
        params.batchprocess = 0;
    else
        params.batchprocess = 1;
    end

    %%
    % Loop over all the files selected for processing:
    
    while FileCount <= NFiles
        tic; % telapsed
        %%
        % Process file and path names (basename, extension, background 

        if (NFiles > 1)
            Files.strInVideoFName = cell2mat(strInVideoFNameArray(FileCount)); 
        else
            Files.strInVideoFName = strInVideoFNameArray; 
        end
        Files.strInVideoFName = Files.strInVideoFName(1:end-4);
        Files.strInFName = [Files.strInVideoPath slashstr Files.strInVideoFName '.' Files.strVideoFExt];
        BackgroundFileName = [Files.strInVideoPath slashstr '' Files.strInVideoFName '_meanimg.dat'];
        fprintf('\n================================\n%s',Files.strInFName);
        if ~nargin, fprintf('\nfile: %g / %g',FileCount,NFiles); end
        fprintf('\n--------------------------------\n');
        if ~params.bool_asked
            set(FigureHandle,'Name',['Drosophila Behavior Recognition - ' Files.strInFName]);
        end

        %%
        % Total number of frames, frame rate, dt
        a = mmread(Files.strInFName);
        intNFrms = abs(a.nrFramesTotal);            % The number of frames to process
        nframes_mean = min(nframes_mean,intNFrms);
        intFps = intNFrms / a.totalDuration;        % Frame rate 
        dt = 1 / intFps;                            % Inter-frame separation
         
        %%
        % Set the appropriate flag if the background file already exists. 

        ResumeBackground = 0;
        BackgroundExist = exist([BackgroundFileName '_win'],'file');
        if BackgroundExist
            disp('>> background image file exists');
            ResumeBackground = 1;
        end

        if ~ResumeBackground

            %%
            % Begin frame capture for background image calculation
            % Only nframes_mean frames starting from intStartFrm are used 
            
            % Randomly select 'nframes_mean' frames
            frmindx = 0;
            randframes = randperm(max(intNFrms/2,nframes_mean));
            randframes = sort(randframes(1:nframes_mean) + intStartFrm);

            try
                mmread(Files.strInFName,randframes, ...
                        [],false,true,'ProcessFrameMeanWin', false);
            catch err
                if ~strcmp(err.message(end-14:end),'STOP PROCESSING')
                    rethrow(err); 
                end
            end

            %%
            % Calculate the actual mean and std images, and adjust them
            nframes_mean = frmindx;
            mean_image = mean_image / nframes_mean;
            std_image = std_image / nframes_mean - mean_image.^2;
            mean_image = 0.11 * mean_image(:,:,1) + ... 
                         0.59 * mean_image(:,:,2) + 0.3 * mean_image(:,:,3);
            std_image = 0.11 * std_image(:,:,1) + ... 
                        0.59 * std_image(:,:,2) + 0.3 * std_image(:,:,3);

            %%
            % Save the result to a file

            [BackgroundFID] = fopen([BackgroundFileName '_win'],'w');
            fwrite(BackgroundFID,size(mean_image,1),'double');
            fwrite(BackgroundFID,size(mean_image,2),'double');
            fwrite(BackgroundFID,mean_image,'double');
            fwrite(BackgroundFID,std_image,'double');
            fclose(BackgroundFID);

        else

            %%
            % Read the background image from the file 
            
            [BackgroundFID] = fopen([BackgroundFileName '_win'],'r');
            nrows = fread(BackgroundFID,1,'double');
            ncols = fread(BackgroundFID,1,'double');
            mean_image = fread(BackgroundFID,[nrows ncols],'double');
            std_image = fread(BackgroundFID,[nrows ncols], 'double');
            fclose(BackgroundFID);
            
        end

        %%
        % Begin frame capture for feature identification

        if ispc
            clear('mexDDGrab.mexw32'); % Clear previous instance of video grabber
        end

        try            
            mmread( Files.strInFName, intStartFrm-1:intStartFrm-1, ...
                [], false, true, 'ProcessFrameCapWin', false );
            
            if params.stopProcess == 1
                break;
            end
            
        catch err
            if ~strcmp(err.message(end-14:end),'STOP PROCESSING')
                rethrow(err); 
            end
        end


        %%
        % Close and save the file containing all the features identified 
        % and errors detected from all chambers in the configuration. 

        for i=1:params.nchambers
            fclose( Files.ErrorFID{i} );
        end

        FileCount = FileCount + 1; 
                
        configureFileName = [Files.strInVideoPath slashstr '' Files.strInVideoFName '.config'];
        save(configureFileName,'params', 'Files', 'NFiles',...
         'DimX', 'DimY', 'DimZ', 'Spacer', 'chamber',...
         'intStartFrm', 'intNFrms', 'dt', 'nframes_mean',...
         'frmindx', 'scale');
  
        %%
        % Display the final time measurement

        telapsed = toc;
        fprintf('time elapsed: %02.0f:%02.0f mm:ss\n', floor(telapsed/60), mod(telapsed,60));

    end %for
        
    %%
    % Close the figure
    close all force;
    
    if params.unprocessedMvNum > 0
        dlgname = 'List of unprocessed movie files';
        warningstring{1} = 'The following movies failed to pass the calibration:';
        for i=1:params.unprocessedMvNum
            warningstring{i+1} = params.unprocessedMvName{i};                
        end
        warningstring{params.unprocessedMvNum+2} = 'Please run qtrak_preprocess for each movie one by one.';
        warndlg(warningstring,dlgname);        
    end
    
