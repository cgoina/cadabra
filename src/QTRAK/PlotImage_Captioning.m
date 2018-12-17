%% PlotImage_Captioning.m
% Copyright (C) 2008 Heiko Dankert, California Institute of Technology

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

function PlotImage_Captioning( data, images, width, height, FigureHandle, FrameNumber ) %#ok<INUSL>
global dt intStartFrm intNFrms;
global Files Panels ObjBuf;
global mean_image std_image mea std;
global v1 v2 scale tcnt;
global ROI;
global object object_1 object_2;
global chamber ic;
global params;
global currentFrame_O1 currentFrame_O2;

    %% A.2 Calibration, ROI Definition, and Initialization
    % This code section is only executed on the starting frame. 
    % The purpose is to calibrate the dimensions and define the regions of 
    % interest, or load them from an existing file. 
    

    if (FrameNumber == intStartFrm-1)
       
        ic = 0;
        Image1 = data;

        [ROI, scale, params] = calibrate( ...
            Files.strInVideoPath, Files.strInVideoFName, ... 
            Image1, mean_image, Panels, FigureHandle, params );

        % The _arr_lim_ variables provide what considered
        % to be the standard fly size in mm (with or without dot, or courtship).
        % This maybe made more flexible by sampling the fly sizes 
        % from frames (or by measuring the flies in advance)
        
        if ~params.bool_dot,
            params.arr_lim = 3.2;
        else
            params.arr_lim = 2.2;
        end
        if params.bool_court,
            params.arr_lim = 3.7;
        end
        
        %% 
        % Initialize various data structures for each chamber. The corners 
        % of these ROIs are stored in _roiCorners_ , which are later passed 
        % back into the frame grabber for frame #2 (where tracking really 
        % begins). Tracked objects are stored in _ObjBuf_.
        
        chamber = struct;
        roiCorners = zeros( params.nchambers * 4 , 1, 'double' );
        mea = cell(params.nchambers,1); 
        std = cell(params.nchambers,1);
        ObjBuf = cell(params.nchambers,1);
        
        for i=1:params.nchambers,
            
            roiCorners( (i-1)*4+1 : (i-1)*4+4 ) = ... 
                [ ROI(i).rows(1), ROI(i).rows(end), ... 
                  ROI(i).cols(1), ROI(i).cols(end) ];
            s_arr = [length(ROI(i).rows) length(ROI(i).cols)];
            cs_ind.r = params.cs : s_arr(1)-params.cs;
            cs_ind.c = params.cs : s_arr(2)-params.cs;
            [x_arr, y_arr] = meshgrid( 1 : s_arr(2) , 1 : s_arr(1) );
            ind_border = find( ...
                ((x_arr >= params.cs) & (x_arr <= s_arr(2)-params.cs)) & ...
                ((y_arr >= params.cs) & (y_arr <= s_arr(1)-params.cs)) );
            x_arr = x_arr * scale.x;
            y_arr = y_arr * scale.y;

            chamber(i).nrows = s_arr(1);
            chamber(i).ncols = s_arr(2);
            chamber(i).s_arr = s_arr;
            chamber(i).cs_ind = cs_ind;
            chamber(i).ind_border = ind_border;
            chamber(i).x_arr = x_arr;
            chamber(i).y_arr = y_arr;

            %% 
            % Calculate the mean and std image for each chamber, as well 
            % as the region of interest by cropping the larger matrices. 
            % Note that these 'images' are shaped as vectors. 
            % The actual mean image is then calculated in (1) below
            % 3*std (standard deviation) helps in cases of mostly 
            % none-moving flies; factor .95 increases the contrast

            mea{i} = mean_image(ROI(i).rows,ROI(i).cols);
            std{i} = std_image(ROI(i).rows,ROI(i).cols);
            mea{i} = 1./(mea{i} + 3 * std{i}); %...............(1)
            
            %%
            % Initialize various tracking variables. The most important 
            % variables are the structures: _object_, _object_1_, and 
            % _object_2_. The structures _object1_ and _object2_ track the 
            % (absolute) motion parameters of the first and the second 
            % flies, respectively, and _object_ tracks the relative motion 
            % parameters between them. Each field is a vector of length 
            % _infNFrms_ that is indexed by the frame number.
            
            tcnt{i} = 1;
            ObjBuf{i} = struct;
            object{i}.dirdiff = 0;
            object{i}.mvdirdiff = 0;
            object{i}.distc = 0;
            object{i}.disth = 0;
            object{i}.distt = 0;
            object{i}.der_distc = 0;
            object{i}.der_disth = 0;
            object{i}.der_distt = 0;
            object_1{i}.disthto2 = 0;
            object_1{i}.pos_x = zeros(1,intNFrms);
            object_1{i}.pos_y = zeros(1,intNFrms);
            object_1{i}.headdir = zeros(1,intNFrms);
            object_1{i}.movedir = 0;
            object_1{i}.to2mvdirdiff = 0;
            object_1{i}.orient = 0;
            object_1{i}.mea = 0;
            object_1{i}.vel = 0;
            object_1{i}.acc = 0;
            object_1{i}.area = 0;
            object_1{i}.Area = 0;
            object_1{i}.length = 0;
            object_1{i}.r = 0;
            object_2{i} = object_1{i};
            
            
            currentFrame_O1{i}.xc = zeros(1,intNFrms);
            currentFrame_O1{i}.yc = zeros(1,intNFrms);
            currentFrame_O1{i}.head = zeros(1,intNFrms);
            currentFrame_O1{i}.xh = zeros(1,intNFrms);
            currentFrame_O1{i}.xt = zeros(1,intNFrms);
            currentFrame_O1{i}.yh = zeros(1,intNFrms);
            currentFrame_O1{i}.yt = zeros(1,intNFrms);
            currentFrame_O1{i}.phir = zeros(1,intNFrms);
            currentFrame_O1{i}.phil = zeros(1,intNFrms);
            currentFrame_O1{i}.r = zeros(1,intNFrms);
            currentFrame_O1{i}.l = zeros(1,intNFrms);
            
            currentFrame_O2{i}.xc = zeros(1,intNFrms);
            currentFrame_O2{i}.yc = zeros(1,intNFrms);
            currentFrame_O2{i}.head = zeros(1,intNFrms);
            currentFrame_O2{i}.xh = zeros(1,intNFrms);
            currentFrame_O2{i}.xt = zeros(1,intNFrms);
            currentFrame_O2{i}.yh = zeros(1,intNFrms);
            currentFrame_O2{i}.yt = zeros(1,intNFrms);
            currentFrame_O2{i}.phir = zeros(1,intNFrms);
            currentFrame_O2{i}.phil = zeros(1,intNFrms);
            currentFrame_O2{i}.r = zeros(1,intNFrms);
            currentFrame_O2{i}.l = zeros(1,intNFrms);            
            
            
            %% 
            % Print the header line into the feature file. Note that there 
            % is a separate feature file for each chamber. 
            
            fprintf( Files.FeatureFID{i}, '%14s', ... 
                'frame'       ,'time [s]'    ,'fly1_dir'    ,'fly2_dir'    , ... 
                'fly1_mvdir'  ,'fly2_mvdir'  ,'fly1_ori'    ,'fly2_ori'    , ... 
                'dir_diff'    ,'mvdir_diff'  ,'1to2_mvdird' ,'2to1_mvdird' , ... 
                'fly1_mean'   ,'fly2_mean'   ,'e_len_1 [mm]','e_len_2 [mm]', ... 
                'e_ar1 [mm2]' ,'e_ar2 [mm2]' ,'vel1 [mm/s]' ,'vel2 [mm/s]' , ...
                'ac1 [mm/s^2]','ac2 [mm/s^2]','c1dist [mm]' ,'c2dist [mm]' , ...
                'dist [mm]'   ,'h1-t2 [mm]'  ,'h2-t1 [mm]'  ,'h-h [mm]'    , ...
                't-t [mm]'    ,'ddist [mm]'  ,'d h-h [mm]'  ,'d t-t [mm]'  , ...
                'Area1 [mm2]' ,'Area2 [mm2]' ,'Length1 [mm]','Length2 [mm]', ...
                'fly_closing' ,'fly_facing'  ,'fly_acting'  ,'fly_conn'    , ...
                'fly1_wingr'  ,'fly1_wingl'  ,'fly1_wing'   ,'fly2_wingr'  , ...
                'fly2_wingl'  ,'fly2_wing'   ,'chpos1 [mm]' ,'chpos2 [mm]' , ...
                'fly1 x [mm]' ,'fly1 y [mm]' ,'fly2 x [mm]' ,'fly2 y [mm]' , ...
                'phi1_r'      ,'phi1_l'      ,'fly1_r [mm]' ,'fly1_l [mm]' , ... 
                'phi2_r'      ,'phi2_l'      ,'fly2_r [mm]' ,'fly2_l [mm]' );
            fprintf(Files.FeatureFID{i},'\n');

        end
        
        %%
        % Send the mean image and ROI corner coordinates, as well as the 
        % indicator of whether or not one of the flies contain a dot, to 
        % the frame grabber, where the information will be used in 
        % preprocessing the pixels.
        
        mearoidat = [Files.strInVideoPath '/' Files.strInVideoFName '_mearoi.mat'];
        save(mearoidat,'mea', 'roiCorners');
        
    %% A.3 Tracking and Analysis 
    % The main analysis and tracking functions are performed in the 
    % following section of the code. If the user specifies live tracking, 
    % then the code begins with refreshing the image, followed with drawing
    % the boundaries of each chamber's ROI. 
    
    else
        
        if (params.bool_plottrack),
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
        end

        %% A.4 The Six-Frame Object Buffer
        %
        % For each chamber, perform analysis and tracking on the current 
        % image, taking the previous measurement as an argument. The
        % results are stored in _ObjBuf's_ fields that are named according 
        % to a convention defined by Heiko as follows. The fields
        % 
        %               000X 00X 0X 1X 11X 111X
        % 
        % correspond to the results measured from the frames grabbed at 
        % T-5, T-4, T-3, T-2, T-1, and T (where T is the current frame in
        % this loop). On frame T, results from T-2 are calculated and 
        % recorded into the feature files. The letter X could either be a 
        % 1 or 2 to indicate the fly number. In the following _for_ loop, 
        % the first five if statements handle the task of populating the 
        % fields 000X, 00X, 0X, 1X, 11X, and 111X. 
        
        for ic = 1 : params.nchambers,
            if (FrameNumber == intStartFrm)
                [ObjBuf{ic}.i0001, ObjBuf{ic}.i0002] = meas_dist( images(ic) );
            elseif (FrameNumber == intStartFrm+1)
                [ObjBuf{ic}.i001, ObjBuf{ic}.i002 ] = meas_dist( images(ic), ... 
                    ObjBuf{ic}.i0001, ObjBuf{ic}.i0002 );
            elseif (FrameNumber == intStartFrm+2)
                [ObjBuf{ic}.i01, ObjBuf{ic}.i02 ] = meas_dist( images(ic), ... 
                    ObjBuf{ic}.i001, ObjBuf{ic}.i002 );
            elseif (FrameNumber == intStartFrm+3)
                [ObjBuf{ic}.i11, ObjBuf{ic}.i12 ] = meas_dist( images(ic), ... 
                    ObjBuf{ic}.i01, ObjBuf{ic}.i02 );
                object_1{ic}.pos_x(FrameNumber) = ObjBuf{ic}.i11.xc;
                object_1{ic}.pos_y(FrameNumber) = ObjBuf{ic}.i11.yc;
                object_2{ic}.pos_x(FrameNumber) = ObjBuf{ic}.i12.xc;
                object_2{ic}.pos_y(FrameNumber) = ObjBuf{ic}.i12.yc;
            elseif (FrameNumber == intStartFrm+4)
                [ObjBuf{ic}.i111, ObjBuf{ic}.i112 ] = meas_dist( images(ic), ...
                    ObjBuf{ic}.i11, ObjBuf{ic}.i12 );
            elseif (FrameNumber == intStartFrm+5)
                [ObjBuf{ic}.i1111, ObjBuf{ic}.i1112 ] = meas_dist( images(ic), ... 
                    ObjBuf{ic}.i111, ObjBuf{ic}.i112 );
            else
                
                %% 
                % From the sixth frame on, all the fields are full. Thus, 
                % the field with the oldest information (i.e., 000X) is 
                % updated with new information; with the effect of shifting
                % all information back to create room for the newest result 
                % to be stored in 111X. 
                
                ObjBuf{ic}.i0001 = ObjBuf{ic}.i001;
                ObjBuf{ic}.i0002 = ObjBuf{ic}.i002;
                ObjBuf{ic}.i001  = ObjBuf{ic}.i01;
                ObjBuf{ic}.i002  = ObjBuf{ic}.i02;
                ObjBuf{ic}.i01   = ObjBuf{ic}.i11;
                ObjBuf{ic}.i02   = ObjBuf{ic}.i12;
                ObjBuf{ic}.i11   = ObjBuf{ic}.i111;
                ObjBuf{ic}.i12   = ObjBuf{ic}.i112;
                ObjBuf{ic}.i111  = ObjBuf{ic}.i1111;
                ObjBuf{ic}.i112  = ObjBuf{ic}.i1112;
                [ObjBuf{ic}.i1111, ObjBuf{ic}.i1112 ] = meas_dist( images(ic), ... 
                     ObjBuf{ic}.i111, ObjBuf{ic}.i112 );
                 
            currentFrame_O1{ic}.xc(FrameNumber) = ObjBuf{ic}.i1111.xc;
            currentFrame_O1{ic}.yc(FrameNumber) = ObjBuf{ic}.i1111.yc;
            currentFrame_O1{ic}.head(FrameNumber) = ObjBuf{ic}.i1111.head;
            currentFrame_O1{ic}.xh(FrameNumber) = ObjBuf{ic}.i1111.xh;
            currentFrame_O1{ic}.xt(FrameNumber) = ObjBuf{ic}.i1111.xt;
            currentFrame_O1{ic}.yh(FrameNumber) = ObjBuf{ic}.i1111.yh;
            currentFrame_O1{ic}.yt(FrameNumber) = ObjBuf{ic}.i1111.yt;
            currentFrame_O1{ic}.phir(FrameNumber) = ObjBuf{ic}.i1111.phir;
            currentFrame_O1{ic}.phil(FrameNumber) = ObjBuf{ic}.i1111.phil;
            currentFrame_O1{ic}.r(FrameNumber) = ObjBuf{ic}.i1111.r;
            currentFrame_O1{ic}.l(FrameNumber) = ObjBuf{ic}.i1111.l;
            
            currentFrame_O2{ic}.xc(FrameNumber) = ObjBuf{ic}.i1112.xc;
            currentFrame_O2{ic}.yc(FrameNumber) = ObjBuf{ic}.i1112.yc;
            currentFrame_O2{ic}.head(FrameNumber) = ObjBuf{ic}.i1112.head;
            currentFrame_O2{ic}.xh(FrameNumber) = ObjBuf{ic}.i1112.xh;
            currentFrame_O2{ic}.xt(FrameNumber) = ObjBuf{ic}.i1112.xt;
            currentFrame_O2{ic}.yh(FrameNumber) = ObjBuf{ic}.i1112.yh;
            currentFrame_O2{ic}.yt(FrameNumber) = ObjBuf{ic}.i1112.yt;
            currentFrame_O2{ic}.phir(FrameNumber) = ObjBuf{ic}.i1112.phir;
            currentFrame_O2{ic}.phil(FrameNumber) = ObjBuf{ic}.i1112.phil;
            currentFrame_O2{ic}.r(FrameNumber) = ObjBuf{ic}.i1112.r;
            currentFrame_O2{ic}.l(FrameNumber) = ObjBuf{ic}.i1112.l;

                %%
                % If counter is not shown on the display panel, provide 
                % progress feedback to the user through the console by 
                % displaying elapsed processing time and movie time every
                % 1800 frames (and only for the first chamber). 
                
                if (~mod(FrameNumber,1800) && ~params.bool_plotcount && (ic == 1)),
                    hou  = floor(( FrameNumber-2)*dt / 3600);
                    minu = floor(((FrameNumber-2)*dt - hou*3600)/60);
                    sec  = floor(((FrameNumber-2)*dt - hou*3600)- minu*60);
                    hun  = floor(((FrameNumber-2)*dt - hou*3600 - minu*60 - sec)*100);
                    t = toc;
                    hou1  = floor( t/3600);
                    minu1 = floor((t-hou1*3600)/60);
                    sec1  = floor( t-hou1*3600 - minu1*60);
                    fprintf('Timestamp: %02d:%02d:%02d:%02d  Systemtime: %02d:%02d:%02d\n', ... 
                        hou,minu,sec,hun,hou1,minu1,sec1);
                end

                %% A.5 Motion Parameter Calculation
                % Having finished the measurements of various object 
                % parameters such as position, size, and so on, the
                % actual motion parameters are calculated based on the 
                % historical information stored in the object fields.
                % First, the coordinates can be used straight away (2). 
                
                % ...................(2)
                object_1{ic}.pos_x(FrameNumber-2) = ObjBuf{ic}.i11.xc;
                object_1{ic}.pos_y(FrameNumber-2) = ObjBuf{ic}.i11.yc;
                object_2{ic}.pos_x(FrameNumber-2) = ObjBuf{ic}.i12.xc;
                object_2{ic}.pos_y(FrameNumber-2) = ObjBuf{ic}.i12.yc;

                %% 
                % *Velocities* are calculated as follows. First, from the six 
                % positions (calculated from six frames), we can compute
                % five position differences, which when divided by dt, will
                % give us five velocity measurements (3). A convolution
                % with a kernel of size 3 gives us three velocities that are 
                % smoother in time. The velocity is then based on the 
                % mid-point value (i.e., 0X and 1X, or T-3 and T-2). 
                
                % ...................(3)
                v1x = conv2([ObjBuf{ic}.i001.xc  - ObjBuf{ic}.i0001.xc, ... 
                             ObjBuf{ic}.i01.xc   - ObjBuf{ic}.i001.xc, ... 
                             ObjBuf{ic}.i11.xc   - ObjBuf{ic}.i01.xc, ... 
                             ObjBuf{ic}.i111.xc  - ObjBuf{ic}.i11.xc, ... 
                             ObjBuf{ic}.i1111.xc - ObjBuf{ic}.i111.xc]/dt, ...
                            [.25 .5 .25],'valid');
                v1y = conv2([ObjBuf{ic}.i001.yc  - ObjBuf{ic}.i0001.yc, ... 
                             ObjBuf{ic}.i01.yc   - ObjBuf{ic}.i001.yc, ... 
                             ObjBuf{ic}.i11.yc   - ObjBuf{ic}.i01.yc, ... 
                             ObjBuf{ic}.i111.yc  - ObjBuf{ic}.i11.yc, ... 
                             ObjBuf{ic}.i1111.yc - ObjBuf{ic}.i111.yc]/dt, ... 
                            [.25 .5 .25],'valid');
                v2x = conv2([ObjBuf{ic}.i002.xc  - ObjBuf{ic}.i0002.xc, ... 
                             ObjBuf{ic}.i02.xc   - ObjBuf{ic}.i002.xc, ... 
                             ObjBuf{ic}.i12.xc   - ObjBuf{ic}.i02.xc, ... 
                             ObjBuf{ic}.i112.xc  - ObjBuf{ic}.i12.xc, ... 
                             ObjBuf{ic}.i1112.xc - ObjBuf{ic}.i112.xc]/dt, ... 
                            [.25 .5 .25],'valid');
                v2y = conv2([ObjBuf{ic}.i002.yc  - ObjBuf{ic}.i0002.yc, ... 
                             ObjBuf{ic}.i02.yc   - ObjBuf{ic}.i002.yc, ... 
                             ObjBuf{ic}.i12.yc   - ObjBuf{ic}.i02.yc, ... 
                             ObjBuf{ic}.i112.yc  - ObjBuf{ic}.i12.yc, ... 
                             ObjBuf{ic}.i1112.yc - ObjBuf{ic}.i112.yc]/dt, ... 
                            [.25 .5 .25],'valid');
                v1 = sqrt(v1x.^2 + v1y.^2);
                v2 = sqrt(v2x.^2 + v2y.^2);
                object_1{ic}.vel = v1(2);
                object_2{ic}.vel = v2(2);

                %% 
                % *Accelerations* are calculated from velocities using the 
                % same principle. From the three velocity values stored 
                % in _v1_ and _v2_, we can obtain a single acceleration 
                % sample by performing a smoothing convolution (4).
                
                % ...................(4)
                ac1 = sqrt((conv2(v1x,[-1 0 1],'valid')/2/dt).^2 + ... 
                           (conv2(v1y,[-1 0 1],'valid')/2/dt).^2);
                ac2 = sqrt((conv2(v2x,[-1 0 1],'valid')/2/dt).^2 + ... 
                           (conv2(v2y,[-1 0 1],'valid')/2/dt).^2);
                object_1{ic}.acc = ac1;
                object_2{ic}.acc = ac2;

                %% 
                % *Distances* between the center, head, and tail of the 
                % two flies are then calculated for time 0X (5a) and 
                % time 1X (5b). The formula used is the standard L2 norm.
                % Next, the differences between the positions at 0X and 
                % 1X, and between one fly's head to the other fly's tail 
                % at 1X are calculated in (5c) and (5d), respectively. 
                
                % ...................(5a)                
                dist0c = dist([ObjBuf{ic}.i01.xc ObjBuf{ic}.i01.yc], ... 
                              [ObjBuf{ic}.i02.xc ObjBuf{ic}.i02.yc]');
                dist0h = dist([ObjBuf{ic}.i01.xh ObjBuf{ic}.i01.yh], ... 
                              [ObjBuf{ic}.i02.xh ObjBuf{ic}.i02.yh]');
                dist0t = dist([ObjBuf{ic}.i01.xt ObjBuf{ic}.i01.yt], ... 
                              [ObjBuf{ic}.i02.xt ObjBuf{ic}.i02.yt]');
                          
                % ...................(5b)                
                dist1c = dist([ObjBuf{ic}.i11.xc ObjBuf{ic}.i11.yc], ... 
                              [ObjBuf{ic}.i12.xc ObjBuf{ic}.i12.yc]');
                dist1h = dist([ObjBuf{ic}.i11.xh ObjBuf{ic}.i11.yh], ... 
                              [ObjBuf{ic}.i12.xh ObjBuf{ic}.i12.yh]');
                dist1t = dist([ObjBuf{ic}.i11.xt ObjBuf{ic}.i11.yt], ... 
                              [ObjBuf{ic}.i12.xt ObjBuf{ic}.i12.yt]');
                
                % ...................(5c)                
                ddistc = dist1c - dist0c;
                ddisth = dist1h - dist0h;
                ddistt = dist1t - dist0t;
                
                % ...................(5d)                
                object_1{ic}.disthto2 = dist( ... 
                              [ObjBuf{ic}.i11.xh ObjBuf{ic}.i11.yh], ... 
                              [ObjBuf{ic}.i12.xt ObjBuf{ic}.i12.yt]');
                object_2{ic}.disthto1 = dist( ... 
                              [ObjBuf{ic}.i12.xh ObjBuf{ic}.i12.yh], ... 
                              [ObjBuf{ic}.i11.xt ObjBuf{ic}.i11.yt]');

                %% 
                % *Directions* are calculated from the center positions 
                % of both flies at 1X. The angles between the flies' heads 
                % and their directions are similarly calculated in (6a). 
                % The results are then adjusted to ensure all angles fall
                % between -180 and 180 degrees (6b).
                
                % ...................(6a)                
                dirdiff = abs(ObjBuf{ic}.i11.head-ObjBuf{ic}.i12.head);
                mv_dirdiff = abs(ObjBuf{ic}.i11.mv_phi-ObjBuf{ic}.i12.mv_phi);
                obj12dir = atan2( ObjBuf{ic}.i12.yc - ObjBuf{ic}.i11.yc, ... 
                                  ObjBuf{ic}.i12.xc - ObjBuf{ic}.i11.xc );
                object_1{ic}.to2mvdirdiff = 180/pi*(ObjBuf{ic}.i11.head - obj12dir);
                object_2{ic}.to1mvdirdiff = 180/pi*(ObjBuf{ic}.i12.head - (obj12dir-pi));
                object_1{ic}.orient = 180/pi*(ObjBuf{ic}.i11.phi);
                object_2{ic}.orient = 180/pi*(ObjBuf{ic}.i12.phi);
                object{ic}.dirdiff = 180/pi*(dirdiff);
                object{ic}.mvdirdiff = 180/pi*(mv_dirdiff);
                
                % ...................(6b)                
                if object_1{ic}.to2mvdirdiff>180, 
                    object_1{ic}.to2mvdirdiff = object_1{ic}.to2mvdirdiff - 360; 
                end
                if object_1{ic}.to2mvdirdiff<(-180), 
                    object_1{ic}.to2mvdirdiff = object_1{ic}.to2mvdirdiff + 360; 
                end
                if object_2{ic}.to1mvdirdiff>180, 
                    object_2{ic}.to1mvdirdiff = object_2{ic}.to1mvdirdiff - 360; 
                end
                if object_2{ic}.to1mvdirdiff<(-180), 
                    object_2{ic}.to1mvdirdiff = object_2{ic}.to1mvdirdiff + 360; 
                end
                if object_1{ic}.orient<0, 
                    object_1{ic}.orient = object_1{ic}.orient + 180; 
                end
                if object_2{ic}.orient<0, 
                    object_2{ic}.orient = object_2{ic}.orient + 180; 
                end
                if object{ic}.dirdiff>180, 
                    object{ic}.dirdiff = object{ic}.dirdiff - 360; 
                end
                if object{ic}.mvdirdiff>180, 
                    object{ic}.mvdirdiff = object{ic}.mvdirdiff - 360; 
                end

                %% A.6 Measurement Output
                % After the motion parameters are calculated, prepare them 
                % for eventual transfer to the feature file. Most results 
                % are just copied directly to the object variables (which 
                % suggests a potential for code cleanup), except for a few 
                % translation from radian into degree. The feature file has
                % a fixed-width format that includes ten empty spaces.

                object_1{ic}.pos_change = ObjBuf{ic}.i11.mv_d;
                object_2{ic}.pos_change = ObjBuf{ic}.i12.mv_d;
                object_1{ic}.headdir(FrameNumber-2) = 180/pi*(ObjBuf{ic}.i11.head);
                object_2{ic}.headdir(FrameNumber-2) = 180/pi*(ObjBuf{ic}.i12.head);
                object_1{ic}.movedir = 180/pi*(ObjBuf{ic}.i11.mv_phi);
                object_2{ic}.movedir = 180/pi*(ObjBuf{ic}.i12.mv_phi);
                object_1{ic}.area = ObjBuf{ic}.i11.area;
                object_2{ic}.area = ObjBuf{ic}.i12.area;
                object_1{ic}.Area = ObjBuf{ic}.i11.Area;
                object_2{ic}.Area = ObjBuf{ic}.i12.Area;
                object_1{ic}.length = ObjBuf{ic}.i11.length;
                object_2{ic}.length = ObjBuf{ic}.i12.length;
                object_1{ic}.mea = ObjBuf{ic}.i11.mean;
                object_2{ic}.mea = ObjBuf{ic}.i12.mean;
                object{ic}.distc = dist1c;
                object{ic}.disth = dist1h;
                object{ic}.distt = dist1t;
                object{ic}.der_distc = ddistc;
                object{ic}.der_disth = ddisth;
                object{ic}.der_distt = ddistt;
                object_1{ic}.r = ObjBuf{ic}.i11.ra;
                object_2{ic}.r = ObjBuf{ic}.i12.ra;

                fprintf( Files.FeatureFID{ic}, ... 
                   ['%14.0f %13.4f %13.0f %13.0f %13.0f %13.0f %13.0f %13.0f %13.0f %13.0f ', ...
                    '%13.0f %13.0f %13.6f %13.6f %13.6f %13.6f %13.6f %13.6f %13.6f %13.6f ', ...
                    '%13.6f %13.6f %13.6f %13.6f %13.6f %13.6f %13.6f %13.6f %13.6f %13.6f ', ...
                    '%13.6f %13.6f %13.6f %13.6f %13.6f %13.6f %13.0f %13.0f %13.0f %13.0f ', ...
                    '%13.0f %13.0f %13.0f %13.0f %13.0f %13.0f %13.6f %13.6f %13.6f %13.6f ', ...
                    '%13.6f %13.6f %13.6f %13.6f %13.6f %13.6f %13.6f %13.6f %13.6f %13.6f'], ...
                    FrameNumber-2,                          (FrameNumber-2)*dt, ...
                    object_1{ic}.headdir(FrameNumber-2),    object_2{ic}.headdir(FrameNumber-2), ...
                    object_1{ic}.movedir,                   object_2{ic}.movedir, ...
                    object_1{ic}.orient,                    object_2{ic}.orient, ...
                    object{ic}.dirdiff,                     object{ic}.mvdirdiff, ...
                    object_1{ic}.to2mvdirdiff,              object_2{ic}.to1mvdirdiff, ...
                    object_1{ic}.mea,                       object_2{ic}.mea, ...
                    2*ObjBuf{ic}.i11.A,                     2*ObjBuf{ic}.i12.A, ...
                    object_1{ic}.area,                      object_2{ic}.area, ...
                    object_1{ic}.vel,                       object_2{ic}.vel, ...
                    object_1{ic}.acc,                       object_2{ic}.acc, ...
                    object_1{ic}.r,                         object_2{ic}.r, ...
                    object{ic}.distc,                       object_1{ic}.disthto2, ...
                    object_2{ic}.disthto1,                  object{ic}.disth, ...
                    object{ic}.distt,                       object{ic}.der_distc, ...
                    object{ic}.der_disth,                   object{ic}.der_distt,...
                    object_1{ic}.Area,                      object_2{ic}.Area, ...
                    object_1{ic}.length,                    object_2{ic}.length,...
                    zeros(1,10),... 
                    object_1{ic}.pos_change,                object_2{ic}.pos_change, ...
                    object_1{ic}.pos_x(FrameNumber-2),      object_1{ic}.pos_y(FrameNumber-2), ...
                    object_2{ic}.pos_x(FrameNumber-2),      object_2{ic}.pos_y(FrameNumber-2), ...
                    180/pi*(ObjBuf{ic}.i11.phir),           180/pi*(ObjBuf{ic}.i11.phil), ...
                    ObjBuf{ic}.i11.r,                       ObjBuf{ic}.i11.l, ...
                    180/pi*(ObjBuf{ic}.i12.phir),           180/pi*(ObjBuf{ic}.i12.phil), ...
                    ObjBuf{ic}.i12.r,                       ObjBuf{ic}.i12.l );
                fprintf( Files.FeatureFID{ic}, '\n' ); 

                %% A.7 Visualization
                % If requested, the following lines would draw markers on 
                % the frame image displayed on the panel. Visualization 
                % incurs heavy performance penalty, and should not be used 
                % in production tracking. Three types of markers are 
                % supported: center position, direction, and motion 
                % parameter history. 

                if (params.bool_plottrack),

                    [tcnt{ic}] = write_sequence(ObjBuf{ic}, object_1{ic}, object_2{ic}, ...
                        ROI(ic), tcnt{ic}, params, scale, FrameNumber, FigureHandle);
                end
            end
        end
    end
    
return
end

%% 
% <html>
% <h1> B. Inner Loop </h1>
% The outer loop activates the inner loop once for every chamber. 
% The function <i>meas_dist()</i> implements <br> the inner loop functionalities. 
% </html>
%
%% B.1 Function Arguments
% The variable _Image0_ is a structure that contains subimages from the 
% chamber currently processed. These subimages are produced by the 
% optimized video grabber, taking into account whether or not one of the 
% fly contains a dot. The variable _obj_1_ and _obj_2_ are previous 
% measurements of the two flies. Other parameters are passed implicitly
% through the global variables _params_, _scale_, _ic_, and _chamber_. 

function [obj1,obj2] = meas_dist( Image0, obj_1, obj_2 )
global scale chamber ic params;

%% B.2 Initialization
% The following two if statements are executed only in the first few
% frames where the input arguments are not supplied (because the objects
% have not been measured yet. In this case, the objects are initialized
% to contain zero values. 

if nargin<2,
    obj_1.xc = 0;         obj_1.yc = 0;           obj_1.area = 0;
    obj_1.ra = 0;         obj_1.mean = 0;         obj_1.mv_d = 0;
    obj_1.mv_phi = 0;     obj_1.head = -99;       obj_1.Area = 0;
    obj_1.length = 0;     obj_1.l = 0;            obj_1.phil = 0;
    obj_1.r = 0;          obj_1.phir = 0;         obj_1.xh = 0;
    obj_1.yh = 0;         obj_1.xt = 0;           obj_1.yt = 0;
    obj_1.xcc = 0;        obj_1.ycc = 0;          obj_1.rcc = 0;
    obj_1.mu = [0 3 5];   obj_1.cov = [1 1 1].^2; obj_1.exc = 0;
    obj_1.eyc = 0;        obj_1.A = 0;            obj_1.B = 0;
    obj_1.phi = 0;
end
if nargin<3,
    obj_2.xc = 0;         obj_2.yc = 0;           obj_2.area = 0;
    obj_2.ra = 0;         obj_2.mean = 0;         obj_2.mv_d = 0;
    obj_2.mv_phi = 0;     obj_2.head = -99;       obj_2.Area = 0;
    obj_2.length = 0;     obj_2.l = 0;            obj_2.phil = 0;
    obj_2.r = 0;          obj_2.phir = 0;         obj_2.xh = 0;
    obj_2.yh = 0;         obj_2.xt = 0;           obj_2.yt = 0;
    obj_2.xcc = 0;        obj_2.ycc = 0;          obj_2.rcc = 0;
    obj_2.mu = [0 3 5];   obj_2.cov = [1 1 1].^2; obj_2.exc = 0;
    obj_2.eyc = 0;        obj_2.A = 0;            obj_2.B = 0;
    obj_2.phi = 0;
end
params.dist2 = 99;
arr_lim = params.arr_lim;

%% B.3-B.5 Fly Body Pixel Segmentation
% In this  step we want to segment the fly body pixel only and find out
% how many flies are present (B.3), and whether flies appear merged or not 
% and may need to be separated (B.4,B.5)
%
% Using an optimized histograming function _histfit_, a binary
% image _bin0_ and a masked image _img_ are computed (7a). Note that
% _bin0_ is smaller than _img_ by _cs_ (all around the border).
%
% The next step is to label the connected components in the binary
% image (7b). For example, if there are two components, the first one
% would be labeled 1, the second one 2, etc. The number of components
% is _l0_ and the labeled image is _lab0_.
%
% Next, the size (in pixel) of these _l0_ components are computed
% and stored in _arr0_. The hope is that at least one (or more) of
% these components correspond to the fly (or flies). This is determined
% by first scaling the component area into mm^2 and finding those
% components that exceed 0.2mm^2 (7c).

% ...................(7a)
img2 = Image0.img2;
[ img, bin0, obj1.mu, obj1.cov ] = histfit( Image0.img, 0, obj_1.mu, obj_1.cov);
obj2.mu = obj1.mu; obj2.cov = obj1.cov;
bin0 = bin0(chamber(ic).cs_ind.r, chamber(ic).cs_ind.c);

% ...................(7b)
[lab0,l0] = bwlabel(bin0);
arr0 = zeros(1,l0);
for i=1:l0,
    arr0(i) = numel(find(lab0 == i));
end

% ...................(7c)
indl0 = find(arr0 * scale.xy > arr_lim/16);
l0 = numel(indl0);

%%
% If there are two connected components and none of the flies are
% tagged with a dot, then perform binary morphological operation on
% the image to clean up the blobs (removing tiny dots, or holes).
% Otherwise, if one of the flies have a dot, then first sort the
% fly components (_indl0_) by its physical area (not pixel area).
% If there are more than two fly components, OR if there are two
% fly components (_l0_ == 0), but there is an indication that the dot
% splits the fly blob into two pieces, perform another type of binary
% morphological operation to glue back the pieces. The number of "flies"
% are updated after either one of these morphological operations.

ind_morph = 0;
% ...................(8a)
if l0 > 2 && ~params.bool_dot,
    ind_morph = 1;
    bin0 = bwmorph(bwmorph(bin0,'dilate',2),'erode',2);
    [lab0,l0] = bwlabel(bin0);
    arr0 = zeros(1,l0);
    for i=1:l0,
        arr0(i) = numel(find(lab0 == i));
    end
    % ...................(8b)
elseif params.bool_dot,
    arr0s = sort(arr0(indl0) * scale.xy);
    if (l0 ~= 2 || (l0 == 2 && (arr0s(end-1)/arr0s(end))<.55)),
        ind_morph = 1;
        bin0 = bwmorph(bwmorph(bwmorph(bin0,'dilate',2),'close',1),'erode',2);
        [lab0,l0] = bwlabel(bin0);
        arr0 = zeros(1,l0);
        for i=1:l0,
            arr0(i) = numel(find(lab0 == i));
        end
    end
end
% ...................(8c)
l0 = numel(find(arr0 * scale.xy > .4));



%% B.4 Merged Flies - Try a simple separation via mophological operations
% The final number of components (_l0_) is the primary indicator of whether
% or not one (or more) fly is detected in the image. First, identify the 
% largest blob (in mm^2 scale). If we only detect one fly, and the blob 
% size is smaller than the size limit, and the user specifically informs 
% the program that there should be more than one fly in the chamber, then
% more processing is required. The hope is to obtain more than one blobs
% larger than 0.4mm^2 (9b), with the largest blob isolated and measured. 

if (l0),
    ind_head = 0;
    max01 = find(arr0 == max(arr0));
    max_arr = arr0(max01(1)) * scale.xy;
    % ...................(9a)
    if (l0 == 1) && (max_arr <= arr_lim) && params.bool_nfly,
        bin0 = bwmorph(bwmorph(bin0,'dilate',2),'erode',2);
        [lab0,l0] = bwlabel(bin0);
        arr0 = zeros(1,l0);
        for i=1:l0, 
            arr0(i) = numel(find(lab0 == i)); 
        end
        % ...................(9b)
        l0 = numel(find(arr0 * scale.xy > .4));
        max01 = find(arr0 == max(arr0));
        max_arr = arr0(max01(1)) * scale.xy;
        arr_lim = arr_lim * 1.05;
    end

    %% B.5 Merged Flies (two flies in One blob) - Separation via GMM
    % If from the previous result, the largest blob is bigger than the size
    % threshold, and the blob is the only one with area exceeding 1mm^2, 
    % and the user explicitly specifies there should be more than one fly, 
    % we have the case of two flies very close to each other. The proposed 
    % solution to this is to perform an 'open' (10a) operation on the blob, 
    % and identify the blobs with area exceeding 0.4 mm^2 (10b).
    
    if( max_arr > arr_lim ) && ( sum(arr0 * scale.xy > 1) == 1 ) && params.bool_nfly,

        % ...................(10a)
        bin01 = bwmorph(bin0,'open');
        [lab01,l01] = bwlabel(bin01);
        arr01 = zeros(1,l01); 
        for i=1:l01, 
            arr01(i) = numel(find(lab01 == i)); 
        end
        
        % ...................(10b)
        l01 = numel(find(arr01 * scale.xy > .4));
        if l01 > l0,
            l0 = l01; 
            bin0 = bin01; 
            lab0 = lab01; 
            arr0 = arr01;
        end

        %%
        % As before, find the largest object, and create the index _f01_ 
        % for this candidate fly (11a). If this blob is bigger than the 
        % threshold, and the blob is the only one with area exceeding 1mm^2, 
        % then indicate that we have found the "head" of the fly, and 
        % create a small rectangle indexed by _fly_ind_. 
        
        % ...................(11a)
        max01 = find(arr0 == max(arr0));
        max_arr = arr0(max01(1)) * scale.xy;
        f01 = chamber(ic).ind_border(lab0 == max01(1));
        if (max_arr > arr_lim) && (sum(arr0 * scale.xy > 1) == 1),
            ind_head = 1;
            fly_ind.r = uint16( ... 
                min(chamber(ic).y_arr(f01)/scale.y) - params.cs + 1 : ... 
                max(chamber(ic).y_arr(f01)/scale.y) + params.cs - 1 );
            fly_ind.c = uint16( ... 
                min(chamber(ic).x_arr(f01)/scale.x) - params.cs + 1 : ... 
                max(chamber(ic).x_arr(f01)/scale.x) + params.cs - 1 );
            
            %%
            % The image enclosed by this small rectangle is then analyzed by
            % performing another masking operation (11b) (or if there is a 
            % dot, just simply copy the whole image. Finally, after a series 
            % of morphological operations, the image is then re-labeled (11c). 

            % ...................(11b)
            if ~params.bool_dot,
                tmp1 = histfit(img(fly_ind.r,fly_ind.c),-99);
            else
                tmp1 = img * 0;
                tmp1(chamber(ic).cs_ind.r,chamber(ic).cs_ind.c) = bin0;
                tmp1 = tmp1(fly_ind.r,fly_ind.c);
            end

            % ...................(11c)
            if ind_morph,
                [lab0_,l0] = bwlabel(bwmorph(bwmorph(tmp1 > 0,'dilate',2),'erode',2));
            else
                [lab0_,l0] = bwlabel(bwmorph(tmp1 > 0,'close'));
            end

            %%
            % Next, the largest blobs that are larger than 0.4mm^2 is then 
            % isolated (12a). If there is only one such blob, we will try 
            % to fit two ellipses on the image (12b). The centers of these 
            % ellipses will then be the proposed centers of the flies. If 
            % this fitting fails, then return the previous results. If the 
            % fit returns a positive result, we now have two flies. Find 
            % the blob with largest area (12c). 

            % ...................(12a)
            arr0 = zeros(1,l0);
            for i=1:l0, arr0(i)=numel(find(lab0_ == i)); end
            l01 = numel(find(arr0 * scale.xy > .4));

            % ...................(12b)
            if (l01 == 1),
                [lab01,STATS01] = mask2ellipses(lab0_,img(fly_ind.r,fly_ind.c).^2,params.bool_dot);
                nume = numel(STATS01);
                l0 = 0;
                lab0_ = lab0_ * 0;
                for i=1:nume,
                    if (STATS01(i).MajorAxisLength > 0) && ... 
                       (STATS01(i).MajorAxisLength < 100),
                        l0 = l0 + 1;
                        % erase smaller split-up segments
                        [a,a0] = bwlabel(lab01 == i);
                        if a0 > 1,
                            ma = zeros(a0,1);
                            for ii=1:a0, ma(ii) = numel(find(a == ii)); end
                            inma = find(ma == max(ma)); inma = inma(1);
                            lab01(a > 0 & a < inma) = 0;
                        end
                        %
                        lab0_(lab01 == i) = l0;
                    end
                end
                if ~l0,
                    obj1 = obj_1; obj2 = obj_2;
                    return;
                end
            end
            
            % ...................(12c)
            lab0 = zeros(chamber(ic).s_arr(1),chamber(ic).s_arr(2));
            lab0(fly_ind.r,fly_ind.c) = lab0_;
            lab0 = lab0(chamber(ic).cs_ind.r,chamber(ic).cs_ind.c);
            arr0 = zeros(l0,1);
            arr0(1:l0) = 0;
            for i=1:l0, 
                arr0(i)=numel(find(lab0 == i)); 
            end
            max01 = find(arr0 == max(arr0));
        end
    end
    % The output of the previous step is _max01_, the pixel indices _f011_ 
    % of the largest blob (the body) of what is now considered to be fly 1.    
    f011 = chamber(ic).ind_border(lab0 == max01(1));

    
    %% B.6 Segmentation of Full Fly Pixel (fly 1)
    
    [f01,fly_ind1] = seg_fullfly_inconstbox(img2,f011,chamber(ic),params,scale,ind_head);
    
    %% 
    % Having identified the first fly, now we have to associate another 
    % blob with the other fly (if the user specifies that there are more 
    % than one fly in the chamber). If there are two blobs with the same 
    % maximum size, out of which one was identified as the first fly, then 
    % the other blob can be considered as the second fly (_max02_). 
    % Otherwise, find the next largest blob. If such a blob exist, then 
    % it is _max02_, from which its pixel indices _f02_ are obtained (14)
    
    if (l0 > 1) && params.bool_nfly,
        
        if (numel(max01) > 1),
            max02 = max01(2);
        else
            max02 = find(arr0 == max(arr0(arr0 < max(arr0))));
            if (numel(max02) > 1), 
                max02 = max02(1); 
            end
        end
        
        % ...................(14)
        f022 = chamber(ic).ind_border(lab0 == max02);
        if isempty(f022), 
            l0 = 1; 
        end
    end
    
    if (l0 > 1) && params.bool_nfly,    
        
        %% Segmentation Full Fly Pixel (fly 2)
        
        [f02,fly_ind2] = seg_fullfly_inconstbox(img2,f022,chamber(ic),params,scale,ind_head);

        %% Test for Intersection of Flies
        
        [f01,f02] = intersect_flies(f01,f02,f011,f022);

        % In the event that fly 2 turns out to be too small 
        % to be a fly (probably flew away or hidden) -
        % take the last known information of this object
        area2 = numel(f022) * scale.xy;
        if (area2 < .5) && ~params.bool_dot,
            if isempty(obj_2),
                obj_2.xc = 0; 
                obj_2.yc = 0; 
                obj_2.head = -99;
            end;
            obj2 = obj_2;
            ind_twoobj = 0;
        else
            ind_twoobj = 1; 
        end
        
    %% 
    % In the event that only one fly is detected, the _head_ variable of 
    % the second object is marked with a special value -99, and the flag 
    % indicating the presence of two flies are set to 0. In addition, the 
    % pixels _f02_ and _f022_ are set to be identical to _f01_ and _f011_. 
    
    else
        if isempty(obj_2),
            obj_2.xc = 0;
            obj_2.yc = 0;
            obj_2.head = -99;
        end;
        obj2 = obj_2;
        ind_twoobj = 0;
        f02 = f01;
        f022 = f011;
    end

    %% B.7 Measurements on Segmented Fly Body Pixel
    
    [obj1,obj_1] = fly_bodymeas(obj1,obj_1,img,f011,fly_ind1,chamber(ic),params,scale,arr_lim,ind_head);
    
    params.dist2 = 0;
    if ind_twoobj && params.bool_nfly,
        [obj2,obj_2] = fly_bodymeas(obj2,obj_2,img,f022,fly_ind2,chamber(ic),params,scale,arr_lim,ind_head);
        params.dist2 = sqrt((obj2.xc - obj1.xc).^2 + (obj2.yc - obj1.yc).^2);
    end
    
    %% B.8 Fly Switchover?

    [obj1,obj2,f01,f02,ind_swap] = switch_flies(obj1,obj2,obj_1,obj_2,f01,f02,img,chamber(ic),params,scale,ind_twoobj);

   
    %% B.9 Head / Tail and Wing Determination
    
    %%% Fly 1
    if ind_twoobj || (~ind_twoobj && ~ind_swap),
        obj1 = fly_headtailwings(obj1,obj_1,img2,f01,chamber(ic),params,ind_head,0);
    end

    %%% Fly 2
    if (ind_twoobj || (~ind_twoobj && ind_swap)) && params.bool_nfly,
        obj2 = fly_headtailwings(obj2,obj_2,img2,f02,chamber(ic),params,ind_head,1);
    end
    
    if ~ind_twoobj,
        if ~ind_swap,
            obj2 = obj_2;            
        else
            obj1 = obj_1;
        end
    end
    
%% B.10 No Fly Detected
% If no fly is detected, save the previous measurement result into the 
% current measurement result and return the objects. 
else
    obj1 = obj_1;
    obj2 = obj_2;
end
end


function [d] = dist(x,y)
    d = sqrt( sum( (x-y').^2 ) );
end
