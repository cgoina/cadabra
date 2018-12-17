%% Calibrate
% Copyright (C) 2008 Heiko Dankert, California Institute of Technology

% This file is part of QTRAK.
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
% * Original implementation by Heiko Dankert
% * Optimization and documentation by Edwin Soedarmadji
% * Changed By Heiko Dankert 11/11/2008
%
%%
% This function performs three major tasks: letting the user (1) calibrate 
% the scaling factor between the pixel unit on the screen, and the actual, 
% physical unit in mm (2) specify the regions of interest that define the 
% experimental chambers to be tracked, OR use an existing calibration file 
% that contains an equivalent information if one with a matching name is 
% found, and (3) specify other running time parameters such as whether or 
% not one of the flies in the chamber has an identifying dot. 

function [ROI,scale,params] = ...
calibrate( path, fname, Image1, mean_image, Panels, FigureHandle, params )
global Files slashstr;
global diam scale;
global intNFrms;    
    %%
    % Initialize file name parameter if not provided by user, and 
    % initialize the calibration file name. 

    s_arr = size(Image1);
    params.bool_cal = 0; % When ROI.mat exists and use one calibraton, this should be 0 in order to avoid recalcuate scale 
    params.loadfile = 'No';
    params.stopProcess = 0;
    if ~params.bool_recal, %if use one calibration
        calfname = params.firstfname; 
    else
        calfname = [path slashstr fname '_1_roi.mat']; 
    end
    ROI_FID = fopen(calfname, 'r');
    
    ROI = struct;

    %%
    % Initialize the calibration prompts and the display panel.

    figure(FigureHandle);
    axes(Panels.HInfo); cla(Panels.HInfo);
    image(Panels.BlankInfoPanel);
    set(Panels.HInfo,'XTick',[],'YTick',[]);
    axes(Panels.HImage); cla(Panels.HImage);
    image(Image1);
    axis image;
    set(Panels.HImage,'XTick',[],'YTick',[]);
    drawnow;
    
    %%
    % Open calibration file, and ask if a new recalibration is needed 
    % The alternative is for the user to reuse another calibration file.
    % Ask if the user wants to use another calibration file, or 
    % or if the user wants to specify the ROI's and calibration manually.

    if (~params.bool_asked),
        params.firstfname = [path slashstr fname '_1_roi.mat']; 
        ROI_FID = fopen(calfname, 'r');
        if (ROI_FID > 0),
            dot = questdlg('Re-Calibrate?','Calibration File Found','Yes','No','Quit','No');
            if strcmp('Yes',dot), 
                params.bool_cal = 1; 
            elseif strcmp('Quit', dot)
                params.QuitPrg = 1;
                params.stopProcess = 1;
                return;
            end    
            
        end
        params.loadfile = questdlg('Load Different Calibration File?',...
            'Load File','Yes','No','Quit','No');
        if strcmp('Quit', params.loadfile)
            params.QuitPrg = 1;
            params.stopProcess = 1;
            return;
        end
    end

    %%
    % If the default calibration file cannot be opened or if recalibration is 
    % requested and user wants to manually perform calibration, then 

    if (ROI_FID < 0 || params.bool_cal) && strcmp(params.loadfile,'No'),
                    
        %%
        % Confirm if it is a Heisenberg configuration, and if it is,
        % determine automatically if it is a one- or two-chamber setup. 
        % The coordinate obtained by getroi is measured in pixels.
        % The scaling factors in mm/pixel are then calculated.

        params.nchambers = questdlg('Heisenberg Chamber?', ...
            'Mode','Yes','No','Quit','Yes');
                    
        if strcmp('Quit', params.nchambers)
                params.QuitPrg = 1;
                params.stopProcess = 1;
                return;
        end
        
        prompt = {'Number of Chamber'};
        dlg_title = 'Calibration';
        def = {num2str(params.chamberNum)};
        params.chamberNum = str2double(cell2mat(inputdlg(prompt,dlg_title,1,def)));
       
        if strcmp(params.nchambers,'Yes'),

            %% Try Auto-Calibration            
            coord = auto_cal(mean_image);
            
            % Did the auto-calibration find the food patch?
            if isstruct(coord),
                if numel(coord.x) == 4, params.nchambers = 2; else params.nchambers = 1; end;
                
                for i=1:params.nchambers,
                    rectangle('Position',[coord.x(2*i-1) coord.y(1) ...
                        coord.x(2*i)-coord.x(2*i-1) coord.y(2)-coord.y(1)],...
                        'EdgeColor','b');
                end
                button = questdlg('Correct?','ROI','Yes','No','Yes');
            else
                button = 'No';
            end

            % Auto-calibration inaccurate?
            if strcmp(button,'No') || ~numel(button),
                bool_auto = 0;
                axes(Panels.HMenu); cla(Panels.HMenu);
                image(Panels.BlankMenuPanel);
                set(Panels.HMenu,'XTick',[],'YTick',[]);
                h = text(round(Panels.PanelDimX/2),105,'SELECT FOOD AREA');
                set(h,'HorizontalAlignment','center','FontSize',8,'color','w');
                h = text(round(Panels.PanelDimX/2),120,'LEFT Vial');
                set(h,'HorizontalAlignment','center','FontSize',8,'color','w');
                figure(FigureHandle); axes(Panels.HImage);
                coord = [];
                coord = getroi(Image1,Panels,FigureHandle,params,coord);
                if coord.x(1) < s_arr(2)/3,
                    params.nchambers = 2;
                else
                    params.nchambers = 1;
                end
            else
                bool_auto = 1;
            end
            dist_pix = coord.x(2) - coord.x(1);
            scale.x = params.dist_mm / dist_pix; 
            scale.y = params.dist_mm / dist_pix;
            scale.xy = scale.x * scale.y; 
            scale.r = sqrt(scale.x^2+scale.y^2);
        else
            params.cs = 4; % this adds 4 pixels onto each side of the circular ROIs. Try changing to 0.
        end
        
        %%
        % If there is only one chamber, ROI can be calculated automatically
        % from the region provided by user, known configuration dimensions, 
        % and scaling factors. The variable ROI is measured in pixels.

        if (params.nchambers == 1),

            bound = cell(1,1); 
            bound{1} = floor( [ 
                coord.x(1) - (47-params.dist_mm)/2/scale.x - params.cs, ...
                coord.y(1) - (37-params.dist_mm)/2/scale.x - params.cs, ...
                47/scale.x + 2*params.cs-1, ... 
                37/scale.x+2*params.cs-1] );
            ROI.rows = bound{1}(2):bound{1}(2)+bound{1}(4);
            ROI.cols = bound{1}(1):bound{1}(1)+bound{1}(3);

        %%
        % If there are two chambers, user needs to specify the distance 
        % between the vials. First, show the image for calibration.
        % Next, let user specify the region (and allow for revision).
        % Finally, calculate the ROI in pixels.

        elseif (params.nchambers == 2),

            if ~bool_auto,
                axes(Panels.HMenu);
                cla(Panels.HMenu);
                image(Panels.BlankMenuPanel);
                set(Panels.HMenu,'XTick',[],'YTick',[]);
                h = text(round(Panels.PanelDimX/2),95,'CLICK on LEFT');
                set(h,'HorizontalAlignment','center','FontSize',8,'color','w');
                h = text(round(Panels.PanelDimX/2),110,'FOOT-BORDER');
                set(h,'HorizontalAlignment','center','FontSize',8,'color','w');
                h = text(round(Panels.PanelDimX/2),125,'of RIGHT Vial');
                set(h,'HorizontalAlignment','center','FontSize',8,'color','w');
                figure(FigureHandle);
                axes(Panels.HImage);
                
                buttonid = 0;
                while ~buttonid,
                    image(Image1);
                    set(Panels.HImage,'XTick',[],'YTick',[]);
                    drawnow;
                    hold on;
                    [x,y] = ginput(1);
                    if (x > 1) && (y > 1),
                        plot([x-5 x+5],[y y],'b','MarkerSize',10);
                        plot([x x],[y-5 y+5],'b','MarkerSize',10);
                    end
                    button = questdlg('Correct?','ROI','Yes');
                    if strcmp(button,'Yes'),
                        buttonid = 1;
                    end
                end
                shift.x = [0 x-coord.x(1)];
            else
                shift.x = [0 coord.x(3)-coord.x(1)];
            end
            
            area_width.x = 37; 
            area_width.y = 47; % [mm]
            shift.y = [0 0] ./ scale.y;


            bound = cell(2,1); 
            for i=1:params.nchambers,
                bound{i} = floor( [ 
                    coord.x(1) + shift.x(i) - (area_width.x - params.dist_mm)/2/scale.x - params.cs, ...
                    coord.y(1) + shift.y(i) - (area_width.y - params.dist_mm)/2/scale.x - params.cs, ...
                    area_width.x/scale.x + 2*params.cs - 1, ... 
                    area_width.y/scale.y + 2*params.cs - 1 ] );
                ROI(i).rows = bound{i}(2):bound{i}(2)+bound{i}(4);
                ROI(i).cols = bound{i}(1):bound{i}(1)+bound{i}(3);
            end

        %%
        % If the configuration has more than two chambers, the user has to 
        % specify the ROI manually by drawing the regions on the screen.

        else

            prompt = {'Diameter [mm] (0 = non-circular/squared)'}; 
            dlg_title = 'Arena size';
            def = {num2str(16)};
            diam = str2double(cell2mat(inputdlg(prompt,dlg_title,1,def)));

            axes(Panels.HMenu); 
            cla(Panels.HMenu);
            image(Panels.BlankMenuPanel);
            set(Panels.HMenu,'XTick',[],'YTick',[]);
            h = text(round(Panels.PanelDimX/2),95,'CALIBRATION');
            set(h,'HorizontalAlignment','center','FontSize',8,'color','w');
            h = text(round(Panels.PanelDimX/2),110,'Select');
            set(h,'HorizontalAlignment','center','FontSize',8,'color','w');
            h = text(round(Panels.PanelDimX/2),125,'2 Points');
            set(h,'HorizontalAlignment','center','FontSize',8,'color','w');
            figure(FigureHandle); 
            axes(Panels.HImage); 

            buttonid = 0;
            xc = [0 0]; 
            yc = [0 0];
            while ~buttonid,
                image(Image1); 
                axis image;
                set(Panels.HImage,'XTick',[],'YTick',[]); 
                drawnow; 
                hold on;
                for i=1:2,
                    [x,y] = ginput(1);
                    if (x > 1) && (y > 1), 
                        plot([x-5 x+5],[y y],'r','MarkerSize',30,'LineWidth',2); 
                        plot([x x],[y-5 y+5],'r','MarkerSize',30,'LineWidth',2); 
                    end
                    xc(i) = x; yc(i) = y;
                end
                button = questdlg('Correct?','ROI','Yes');
                if strcmp(button,'Yes'), 
                    buttonid = 1; 
                end
            end

            prompt = {'Distance [mm]'}; 
            dlg_title = 'Calibration';
            def = {num2str(params.distcirc_mm)};
            params.dist_mm = str2double(cell2mat(inputdlg(prompt,dlg_title,1,def)));

            dist_pix = sqrt((xc(2)-xc(1)).^2+(yc(2)-yc(1)).^2);
            scale.x = params.dist_mm / dist_pix; 
            scale.y = params.dist_mm / dist_pix;
            scale.xy = scale.x * scale.y; 
            scale.r = sqrt(scale.x^2+scale.y^2);

            if diam > 0,
                % DETECT CIRCULAR ARENAS (ALSO COVERS SQUARED ARENAS)
                diam_pix = diam/scale.x;
                coord = auto_cal_circ(mean_image,Image1,diam_pix,Panels,FigureHandle);
            else
                % RECTANGULAR ARENAS (draw ROIs manually)
                axes(Panels.HMenu);
                cla(Panels.HMenu);
                image(Panels.BlankMenuPanel);
                set(Panels.HMenu,'XTick',[],'YTick',[]);
                h = text(round(Panels.PanelDimX/2),105,'DRAW ROIs');
                set(h,'HorizontalAlignment','center','FontSize',8,'color','w');
                h = text(round(Panels.PanelDimX/2),120,'around Arenas');
                set(h,'HorizontalAlignment','center','FontSize',8,'color','w');
                figure(FigureHandle); axes(Panels.HImage); cla(Panels.HImage);
                button = [];
                coord = [];
                while ~strcmp(button,'Finish'),
                    [coord,button] = getroi(Image1,Panels,FigureHandle,params,coord);
                end
            end
            
            params.nchambers = size(coord,2);
            
            if params.chamberNum ~= params.nchambers
                params.unprocessedMvNum = params.unprocessedMvNum + 1;
                params.unprocessedMvName {params.unprocessedMvNum} = Files.strInFName;
                params.stopProcess = 1;
                return;
            end
            
            %%
            % Parameters and bounds are calculated here:

            bound = cell(length(coord),1);
            for i=1:length(coord),
                bound{i} = floor([
                    coord(i).x(1) - params.cs ... 
                    coord(i).y(1) - params.cs ...
                    coord(i).x(2) - coord(i).x(1) + 2*params.cs ...
                    coord(i).y(2) - coord(i).y(1) + 2*params.cs ]);
            end

        end
        
        % Clear menu and image
        axes(Panels.HMenu); cla(Panels.HMenu); image(Panels.BlankMenuPanel);        
        set(Panels.HMenu,'XTick',[],'YTick',[]);
        figure(FigureHandle); axes(Panels.HImage); image(Image1);
        set(Panels.HImage,'XTick',[],'YTick',[]);
        drawnow;
    %%
    % If the default calibration file is found, and manual calibration is 
    % not requested, or if user specifically wants to use a calibration file
    
    else
        if ~params.bool_asked,
            if strcmp(params.loadfile,'Yes'),
                [fname, path] = ...
                    uigetfile( {'*_roi.mat','Calibration file'}, 'Open calibration file');
                calfname = [path fname];
                params.firstfname = calfname;
            end
        end
        
        %params.nchambers = length(dir([calfname(1:end-9) '*_roi.mat']));
        
        %JL10012009 According to Eric's requirement, using one calibration
        %means all movie share the scale but not ROI. So for the second and
        %later movies, ROI should be calculated seperatedly.
        % The .mat file contains the previous definition and data for ROI.
        diam_pix = diam/scale.x;
        coord = auto_cal_circ(mean_image,Image1,diam_pix,Panels,FigureHandle);
        
        params.nchambers = size(coord,2);
        
        if params.chamberNum ~= params.nchambers
            params.unprocessedMvNum = params.unprocessedMvNum + 1;
            params.unprocessedMvName {params.unprocessedMvNum} = Files.strInFName;
            params.stopProcess = 1;
            return;
        end
        %%
        % Parameters and bounds for circular arena ROIs are calculated here:
        
        bound = cell(length(coord),1);
        for i=1:length(coord),
            bound{i} = floor([
                coord(i).x(1) - params.cs ...
                coord(i).y(1) - params.cs ...
                coord(i).x(2) - coord(i).x(1) + 2*params.cs ...
                coord(i).y(2) - coord(i).y(1) + 2*params.cs ]);
        end
        
        
    end

    %%
    % Adjust the bounds for display. Overlay the ROIs on the image, and 
    % finally return the ROIs. 
    
    if params.nchambers == 2,
        diff = bound{1}(1:2) - 2;
        if diff(1) < 0,
            for i=1:params.nchambers,
                bound{i}(1) = bound{i}(1) - diff(1);
                bound{i}(3) = bound{i}(3) + 2*diff(1);
            end
        end
        if diff(2) < 0,
            for i=1:params.nchambers,
                bound{i}(2) = bound{i}(2) - diff(2);
                bound{i}(4) = bound{i}(4) + 2*diff(2);
            end
        end

        diff = s_arr(2) - (bound{end}(1) + bound{end}(3));
        if diff < 0,
            for i=1:params.nchambers,
                bound{i}(1) = bound{i}(1) - diff;
                bound{i}(3) = bound{i}(3) + 2*diff;
            end
        end
        diff = s_arr(1) - (bound{end}(2) + bound{end}(4));
        if diff < 0,
            for i=1:params.nchambers,
                bound{i}(2) = bound{i}(2) - diff;
                bound{i}(4) = bound{i}(4) + 2*diff;
            end
        end
    else
        % added by H Dankert, 11/11/2008
        for i=1:params.nchambers,
            if bound{i}(1) < 1, bound{i}(1) = 1; end
            if bound{i}(2) < 1, bound{i}(2) = 1; end
            diff = s_arr(2) - (bound{i}(1) + bound{i}(3));
            if diff < 0, bound{i}(3) = bound{i}(3) + diff; end
            diff = s_arr(1) - (bound{i}(2) + bound{i}(4));
            if diff < 0, bound{i}(4) = bound{i}(4) + diff; end
        end
    end

    r = struct('rows',{},'cols',{});
    for i=1:params.nchambers,
        if params.bool_boundbox,
            text(bound{i}(1)+10,bound{i}(2)+10,num2str(i),'Color','b');
            rectangle('Position',bound{i},'EdgeColor','b');
        end
        ROI = [];
        ROI.rows = bound{i}(2):bound{i}(2)+bound{i}(4);
        ROI.cols = bound{i}(1):bound{i}(1)+bound{i}(3);
        r(i) = ROI;
    end
    ROI = r;    

% save roi.mat files
%    save ('params.mat', params)
    %%
    % Save each region's region of interest into a separate file.
    % In addition, create separate feature and error files. 
    
%      FeatureFileName = cell(params.nchambers,1); 
     ErrorFileName = cell(params.nchambers,1);
%      Files.FeatureFID = cell(params.nchambers,1); 
     Files.ErrorFID = cell(params.nchambers,1);

    for i=1:params.nchambers,
        
        tmp = ROI;
        ROI = tmp(i); %#ok<NASGU>
 %        FeatureFileName{i} = [Files.strInVideoPath slashstr Files.strInVideoFName '_' num2str(i) '.feat'];
         ErrorFileName{i} = [Files.strInVideoPath slashstr Files.strInVideoFName '_' num2str(i) '.err'];
         save( [ErrorFileName{i}(1:end-4) '_roi.mat'], 'ROI', 'scale');
        ROI = tmp;

%         Files.FeatureFID{i} = fopen(FeatureFileName{i},'w');
%         if Files.FeatureFID{i} == -1,
%             disp('ERROR: The Following Feature file could not be Created');
%             disp(['      s   NAME: ', FeatureFileName{i}]);
%             disp(['         PATH: ', Files.strInVideoPath slashstr]);
%             return
%         end

        Files.ErrorFID{i} = fopen(ErrorFileName{i},'w');
        if Files.ErrorFID{i} == -1,
            disp('ERROR: The Following Error file could not be Created');
            disp(['         NAME: ', ErrorFileName{i}]);
            disp(['         PATH: ', Files.strInVideoPath slashstr]);
            return
        end
        
    end
    
    if ~params.bool_asked,
        
        %select analysis area
        [params.analysis.radius params.analysis.radius0] = ...
            border( );
        
        params.analysis.border_width = params.analysis.radius0 - ...
            params.analysis.radius;
        
        answer = questdlg('Use one Calibration?','Calibration','Yes','No','Quit','Yes');
        if strcmp('Yes',answer),
            params.bool_recal = 0;
        elseif strcmp('Quit', answer)
            params.QuitPrg = 1;
            params.stopProcess = 1;
            return;
        else
            params.bool_recal = 1;
        end
        params.bool_asked = 1;
        answer = questdlg('Does one Fly have a Dot?','Dotted Fly?','Yes','No','Quit','Yes');
        if strcmp('Yes',answer),
            params.bool_dot = 1;
            params.bool_nfly = 1;
            params.bool_court = 0;
        elseif strcmp('Quit', answer)
            params.QuitPrg = 1;
            params.stopProcess = 1;
            return;
        else
            params.bool_dot = 0;
            answer = questdlg('Female-Male?','Female&Male?','Yes','No','Quit','Yes');
            if strcmp('Yes',answer),
                params.bool_court = 1;
                params.bool_nfly = 1;
            elseif strcmp('Quit', answer)
                params.QuitPrg = 1;
                params.stopProcess = 1;
                return;
            else
                params.bool_court = 0;
                answer = questdlg('2 flies/chamber?','One Fly?','Yes','No','Quit','Yes');
                if strcmp('Yes',answer),
                    params.bool_nfly = 1;
                elseif strcmp('Quit', answer)
                    params.QuitPrg = 1;
                    params.stopProcess = 1;
                    return;
                else
                    params.bool_nfly = 0;
                end
            end
        end
        
        prompt={'Enter the tuning threshold:'};
        name='Tuning threshold';
        numlines=1;
        defaultanswer={'0.70'};
        answer=inputdlg(prompt,name,numlines,defaultanswer);
        params.analysis.tuningthreshold = str2num(answer{1});
        
        answer = questdlg('Correct fly orientation?','Dection parameters','Yes','No','Quit','Yes');
        if strcmp('Yes',answer),
            params.analysis.correct_orient = 1;
        elseif strcmp('Quit', answer)
            params.QuitPrg = 1;
            params.stopProcess = 1;
            return;
        else
            params.analysis.correct_orient = 0;
        end
        
        answer = questdlg('Correct fly position?','Dection parameters','Yes','No','Quit','Yes');
        if strcmp('Yes',answer),
            params.analysis.correct_positions = 1;
        elseif strcmp('Quit', answer)
            params.QuitPrg = 1;
            params.stopProcess = 1;
            return;
        else
            params.analysis.correct_positions = 0;
        end
        
        params.analysis.max_frames = intNFrms;
        
    end

end

%% 
% This subroutine simply allows the user to draw regions of interest on 
% the screen, confirm that the regions are indeed what he/she wants, and 
% return the result to _calibrate_.

function [coord,button] = getroi( Image1, Panels, FigureHandle, params, coord ) %#ok<INUSL>

%%
% If the coordinate structure is empty, initialize it with 
% empty arrays, and set the length to zero. Otherwise, 
% update the length to reflect the number of coordinates stored

if isempty(coord),
    coord.x = [];
    coord.y = [];
    ncoord = 0;
else
    ncoord = length(coord);
end

%%
% Execute the following loop until user is satisfied with
% the rectangle created on the screen, which is captured 
% by Matlab's function getrect() and stored in the variable bnd.

buttonid = 0;
while ~buttonid,
    image(Image1); 
    for i=1:ncoord,
        bnd = [ coord(i).x(1) 
                coord(i).y(1) 
                coord(i).x(end) - coord(i).x(1) 
                coord(i).y(end) - coord(i).y(1) ];
        rectangle('Position',bnd,'EdgeColor','b');
        text(bnd(1)+10,bnd(2)+10,num2str(i),'Color','b');
    end
    set(Panels.HImage,'XTick',[],'YTick',[]); 
    drawnow;
    bnd = getrect(FigureHandle);
    if (bnd(3) > 1) && (bnd(4) > 1), 
        rectangle('Position',bnd,'EdgeColor','b'); 
        text(bnd(1)+10,bnd(2)+10,num2str(ncoord+1),'Color','b');
    end
    button = questdlg('Correct?','ROI','Finish','Yes','No','Finish');
    if strcmp(button,'Yes') || strcmp(button,'Finish'), 
        buttonid = 1; 
    end
end

%%
% The rectangle corners are then added into the 
% coordinate array as a new region of interest 

coord(ncoord+1).x = [ bnd(1), bnd(1)+bnd(3) ];
coord(ncoord+1).y = [ bnd(2), bnd(2)+bnd(4) ];

end