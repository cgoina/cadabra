%% updateFrame
% Copyright (C) 2009 Janelia Farm Reserach Campus 

% This file is part of QTRAK
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
% * Original implementation by Heiko Dankert 2006/2007
% * Optimization and documentation by Edwin Soedarmadji
% * Changes and further optimization by Heiko Dankert 11/2008
% 
% This routine is a callback function that performs the main object 
% tracking and analysis functions. The function is called from the 
% frame grabber for every frame that are selected through the 
% _mmread_ function. When finished, this function returns control
% back to the frame grabber. 
%
%% 
% <html>
% <h1> A. Outer Loop </h1>
% </html>
%% A.1 Function Arguments
% The variable _data_ contains the full image. The variable _images_ 
% contains the processed chamber images, while _width_ and _image_ 
% measures the dimension of _data_. To allow this function to update 
% the display with the image and annotations from the current frame, 
% the variable _FigureHandle_ is passed from the _ProcessFrameCapWin_ 
% function, along with the variable _FrameNumber_ that keeps track of 
% the current frame number. 

function updateFrame( data, images, width, height, FigureHandle, FrameNumber )
global intStartFrm;
global Files Panels;
global scale tcnt;
global ROI;
global object_1 object_2;
global ic;
global params;
global currentFrame_O1 currentFrame_O2;
global slashstr;
global currentFrame aviobj a;

    %% A.2 Calibration, ROI Definition, and Initialization
    % This code section is only executed on the starting frame. 
    % The purpose is to calibrate the dimensions and define the regions of 
    % interest, or load them from an existing file. 
    if (FrameNumber == intStartFrm-1)
        calfname = [Files.strInVideoPath slashstr '' Files.strInVideoFName '*_roi.mat'];
        
        params.nchambers = length(dir(calfname));
        
        for i=1:params.nchambers,         
            load([calfname(1:end-9) '_' num2str(i) '_roi.mat']);
            temp(i) = ROI;
            ROI = temp;
            tcnt{i} = 1;
        end
        

        
        %% A.3 Tracking and Analysis
        % The main analysis and tracking functions are performed in the
        % following section of the code. If the user specifies live tracking,
        % then the code begins with refreshing the image, followed with drawing
        % the boundaries of each chamber's ROI.
    else
        mov = ...
            struct('cdata', zeros(a.height, a.width, 3, 'uint8'),...
            'colormap', []);
        
        
        figure(FigureHandle); axes(Panels.HImage); cla(Panels.HImage);
        Image1 = data;
        image(Image1);
        axis image;
        if (params.bool_boundbox),
            for i=1:params.nchambers,
                pos = [ ROI(i).cols(1), ROI(i).rows(1), ...
                    ROI(i).cols(end) - ROI(i).cols(1), ...
                    ROI(i).rows(end) - ROI(i).rows(1)];
                rectangle('Position',pos,'EdgeColor','b');
                text(pos(1)+10,pos(2)+10,num2str(i),'Color','b');
            end
        end
        set(Panels.HImage,'XTick',[],'YTick',[]); hold on;

        
        %% A.4 The Six-Frame Object Buffer
        % update ObjBuf
        if (params.bool_plottrack),
            for ic = 1 : params.nchambers,             
                if (FrameNumber > intStartFrm+5)
                    ObjBuf{ic}.i1111.xc = currentFrame_O1{ic}.xc(FrameNumber);
                    ObjBuf{ic}.i1111.yc = currentFrame_O1{ic}.yc(FrameNumber);
                    ObjBuf{ic}.i1111.head = currentFrame_O1{ic}.head(FrameNumber);
                    ObjBuf{ic}.i1111.xh = currentFrame_O1{ic}.xh(FrameNumber);
                    ObjBuf{ic}.i1111.xt = currentFrame_O1{ic}.xt(FrameNumber);
                    ObjBuf{ic}.i1111.yh = currentFrame_O1{ic}.yh(FrameNumber);
                    ObjBuf{ic}.i1111.yt = currentFrame_O1{ic}.yt(FrameNumber);
                    ObjBuf{ic}.i1111.phir = currentFrame_O1{ic}.phir(FrameNumber);
                    ObjBuf{ic}.i1111.phil = currentFrame_O1{ic}.phil(FrameNumber);
                    ObjBuf{ic}.i1111.r = currentFrame_O1{ic}.r(FrameNumber);
                    ObjBuf{ic}.i1111.l = currentFrame_O1{ic}.l(FrameNumber);
                    
                    ObjBuf{ic}.i1112.xc = currentFrame_O2{ic}.xc(FrameNumber);
                    ObjBuf{ic}.i1112.yc = currentFrame_O2{ic}.yc(FrameNumber);
                    ObjBuf{ic}.i1112.head = currentFrame_O2{ic}.head(FrameNumber);
                    ObjBuf{ic}.i1112.xh = currentFrame_O2{ic}.xh(FrameNumber);
                    ObjBuf{ic}.i1112.xt = currentFrame_O2{ic}.xt(FrameNumber);
                    ObjBuf{ic}.i1112.yh = currentFrame_O2{ic}.yh(FrameNumber);
                    ObjBuf{ic}.i1112.yt = currentFrame_O2{ic}.yt(FrameNumber);
                    ObjBuf{ic}.i1112.phir = currentFrame_O2{ic}.phir(FrameNumber);
                    ObjBuf{ic}.i1112.phil = currentFrame_O2{ic}.phil(FrameNumber);
                    ObjBuf{ic}.i1112.r = currentFrame_O2{ic}.r(FrameNumber);
                    ObjBuf{ic}.i1112.l = currentFrame_O2{ic}.l(FrameNumber);
                    
                    %update tracking 
                    
                    [tcnt{ic}] = drawTracking(ObjBuf{ic}, object_1{ic}, object_2{ic}, ...
                        ROI(ic), tcnt{ic}, params, scale, FrameNumber, FigureHandle);
                end
            end
        end
        
        mov = capturescreen(gcf);
        aviobj = addframe(aviobj,mov);
        currentFrame = currentFrame + 1;
    end
    return
end


