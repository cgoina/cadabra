classdef gvPlugin < handle
    %GVPLUGIN Abstract class with many useful functions for a 
    %   gVision Plugin
properties (Abstract)
    vidFileName
    vidFilePtr
    fig
    subROI
end
    
methods (Static, Abstract)
    % This function must return the name of your plugin
    [cName dispName description] = gvPlugName
end

    

methods
    %When writing the fullframe (no sub-ROIs)
    function writeFullROI(plug,frame,time,metadata)
        gui = get(plug.fig,'userdata');
        data = get(gui.vi,'userdata');

%         [frame time metadata] = getdata(gui.vi,1);
        if data.forceGray & (size(frame,3)>1) %colapse the image to grayscale
            frame = sum(frame/size(frame,3),3,'native');
        end
        data.fc = data.fc+1;

        if data.timestamp  %write timestamp onto frame if enabled
            txt = [datestr(metadata.AbsTime),' - ',num2str(metadata.FrameNumber)];
            txtimg = plugins.text2im(txt);
            for i = 1:size(frame,3)
                frame(1:size(txtimg,1),1:size(txtimg,2),i) = txtimg*255;
            end
        end

        %Write frame
        [dn fn ext]=fileparts(data.file.fname);
        switch ext
            case '.avi'
                writeVideo(data.logger.aviobj,frame)
            case '.tif'
                if data.fc~=1
                    %Must create a directory within the TIFF File for
                    %each image
                    data.logger.ts.writeDirectory();
                    data.logger.ts.setTag(data.tagStruct);
                end
                data.logger.ts.write(frame);
            case '.fmf'
%                     frame = rot90(frame,-1);
                fwrite(data.logger.fid,time,'double');
                fwrite(data.logger.fid,frame(:),'uint8');
        end


        preTime = get(gui.RecStatus,'userdata');
        set(gui.RecStatus,'string',sprintf('nFrames: %u\neTime: %.2f s\nfps: %.2f Hz',data.fc,time,1/(time-preTime)),'userdata',time)
        set(gui.vi,'userdata',data)
    end

    %When chopping a frame into individual subregions
    function writeSubROI(plug,gFrame,time,metadata)
        gui = get(plug.fig,'userdata');
        data = get(gui.vi,'userdata');
        data.fc = data.fc+1;

%         [gFrame time] = getdata(gui.vi,1);
        if data.forceGray & (size(frame,3)>1) %colapse the image to grayscale
            frame = sum(frame/size(frame,3),3,'native');
        end

        for i = 1:length(data.sROI)
            ROI = data.sROI(i).ROI;
            [dn fn ext]=fileparts(data.file.fname);
            frame = gFrame(ROI(2):(ROI(2)+ROI(4)-1),ROI(1):(ROI(1)+ROI(3)-1),:);
            switch ext
                case '.avi'
                    writeVideo(data.logger(i).aviobj,frame);
                case '.tif'
                    if data.fc~=1
                        data.logger(i).ts.writeDirectory();
                        data.tagStruct.ImageLength = ROI(4);
                        data.tagStruct.ImageWidth = ROI(3);
                        data.logger(i).ts.setTag(data.tagStruct);
                    end
                    data.logger(i).ts.write(frame);
                case '.fmf'
                    %frame = rot90(frame,-1);
                    fwrite(data.logger(i).fid,time,'double');
                    fwrite(data.logger(i).fid,frame(:),'uint8');
            end


        end

        preTime = get(gui.RecStatus,'userdata');
        set(gui.RecStatus,'string',sprintf('nFrames: %u\neTime: %.2f s\nfps: %.2f Hz',data.fc,time,1/(time-preTime)),'userdata',time)
        set(gui.vi,'userdata',data)
    end

    %Stop function for the video input object, wraps up remaining
    %frames and closes video
    function wrapupVid(plug)
        gui = get(plug.fig,'userdata');
        data = get(gui.vi,'userdata');
        
        %clean up remaining frames
        f=gui.vi.framesacquiredfcn;
        while gui.vi.framesavailable>0
            feval(f,plug);
        end

        closeVidFiles(plug)
        stop(gui.vi)
        set(gui.StartStop,'value',0,'string','Start')
        StartStop(plug);
        
        if get(gui.MetaDataClear,'value')
            set(gui.MetaDataTxt,'string','')
            for i = 1:length(gui.sROIControls)
                ui = get(gui.sROIControls(i),'userdata');
                set(ui.notes,'string','')
            end
        end
        set(gui.RecStatus,'backgroundcolor',[1 .7 .7])
        disp([num2str(data.fc),' Frames Logged to Disk in each of ',num2str(length(data.logger)),' Files'])

        activateGvUI(plug);
    end
    %Close the video files opened for this capture
    function closeVidFiles(plug)
        gui = get(plug.fig,'userdata');
        data = get(gui.vi,'userdata');
        for i = 1:length(data.logger)
            switch data.file.fname(end-2:end)
                case 'avi'
                    close(data.logger(i).aviobj);
                case 'tif'
                    data.logger(i).ts.close();
                case 'fmf'
                    fseek(data.logger(i).fid,20,-1);
                    fwrite(data.logger(i).fid,data.fc,'uint64');
                    fclose(data.logger(i).fid);
            end
        end
    end
    
    %Functions for turning the gVision UI on/off
    function deactivateGvUI(plug)
        gui = get(plug.fig,'userdata');
        set([gui.Params, gui.gROI, gui.gROIrst, gui.addSROI, gui.gROIdraw, gui.remSROI,...
        gui.ModeOptions, gui.fileSelect, gui.FPT, gui.TR, gui.FPS, gui.snapShot gui.sROIpopDown],'enable','off')
        set(gui.streamMen(2:5),'enable','off')
    end
    function activateGvUI(plug)
        gui = get(plug.fig,'userdata');
        set([gui.Params, gui.addSROI gui.remSROI, ...
        gui.ModeOptions, gui.fileSelect, gui.FPT, gui.TR, gui.FPS, gui.snapShot],'enable','on')
        if isempty(gui.sROIControls)
            set([gui.gROI, gui.gROIrst, gui.gROIdraw],'enable','on')
        else
            set(gui.sROIpopDown,'enable','on')
        end
        set(gui.streamMen(2:5),'enable','on')
    end
    
    
     %return logger struct of data files
    %this should also log associated notes with files from notes field
    function data = openVidFiles(plug)
        gui = get(plug.fig,'userdata');
        
        file = get(gui.fName,'userdata');
        data.sROI = get(gui.listSROI,'userdata');

        %Log only the capture enabled subROI
        dind = [];
        for i = 1:length(gui.sROIControls)
            ui = get(gui.sROIControls(i),'userdata');
            data.sROI(i).notes = get(ui.notes,'string');
            if get(ui.capture,'value')==0
                dind = [dind,i];
            end
        end
        data.sROI(dind)=[];

        %Capture info about video stream 
        data.res = gui.vi.ROIPosition;
        data.res = data.res(3:4);
        data.ind = 1; %Counter for 2GB AVI limit
        data.fc = 0; %Frame Count
        if strcmpi(get(gui.streamMen(3),'checked'),'on')
            data.timestamp = 1;
        else data.timestamp = 0;
        end
        if strcmpi(get(gui.streamMen(4),'checked'),'on')
            data.forceGray = 1;
        else data.forceGray = 0;
        end


        %Some camera modes don't present a FPS (because it
        %might be data rate or exposure limited in F7 mode for
        %DCAM).  In this case, the user must specify frame rate
        if strcmpi('on',get(gui.FPS,'visible'))
            data.FPS = get(gui.FPS,'userdata');
        else
            data.FPS = str2double(get(gui.vi.Source,'FrameRate'));
        end
        
        %Iterate File Name
        if strcmpi(get(gui.streamMen(2),'checked'),'on')
            file = plug.checkFileName(file);
        end
        data.file = file;
    
        %Write global metadata
        md = get(gui.MetaDataTxt,'string');
        if get(gui.MetaDataChk,'value') && ~isempty(md)
            [gg fname e] = fileparts(file.fname);

            md_name = [fname,'.txt'];
            a = fopen(md_name,'w');
            fprintf(a,'%s\r\n',datestr(now));
            fprintf(a,'%s\r\n\r\n',file.fname);
            for fi = 1:size(md,1)
                fprintf(a,'%s\r\n',md(fi,:));
            end
            fclose(a);

        end

        % Determine compression method selected from list
        ch = get(gui.streamMen(5),'children');
        for i=1:length(ch)
            if strcmpi(get(ch(i),'checked'),'on')
                comp = get(ch(i),'Label');
                compID = get(ch(i),'userdata');
            end
        end
        
        [p n ftype] = fileparts(file.fname);

        switch ftype
            case '.avi' %Uncompressed AVIs - 2GB Limit (old)                
                %Single AVI
                if isempty(data.sROI)
                    ff = fullfile(file.pname,file.fname);
                    logger.aviobj = VideoWriter(ff,comp);
                    logger.aviobj.FrameRate = data.FPS;
                    open(logger.aviobj);
                else
                    %Multiple AVIs sROI
                    for i = 1:length(data.sROI)
                        f = file;
                        [gg fname e] = fileparts(file.fname);
                        f.fname = sprintf('%s_%s%s',fname,data.sROI(i).name,e);
                        %Write individual sub-ROI metadata
                        if ~isempty(data.sROI(i).notes)
                            notes_name = sprintf('%s_%s%s',fname,data.sROI(i).name,'.txt');
                            a = fopen(notes_name,'w');
                            fprintf(a,'%s\r\n',f.fname);
                            fprintf(a,'%s\r\n\r\n',datestr(now));
                            for fi = 1:size(data.sROI(i).notes,1)
                                fprintf(a,'%s\r\n',data.sROI(i).notes(fi,:));
                            end
                            fclose(a);
                        end

                        ff = fullfile(file.pname,f.fname);
                        data.sROI(i).fname = ff;
                        logger(i).aviobj = VideoWriter(ff,'Uncompressed AVI');
                        logger(i).aviobj.FrameRate = data.FPS;
                        open(logger(i).aviobj)
                        data.sROIsize(i) = data.sROI(i).ROI(3)*data.sROI(i).ROI(4);
                    end
                end
            case '.fmf' %Dickinson Lab "Fly Movie Format" version 1
                if isempty(data.sROI)
                    %single FMF ROI
                    ff = fullfile(file.pname,file.fname);
                    logger.fid = fopen(ff,'w');
                    fwrite(logger.fid,1,'uint32');
                    fwrite(logger.fid,data.res(1),'uint32');
                    fwrite(logger.fid,data.res(2),'uint32');
                    fwrite(logger.fid,data.res(1)*data.res(2)+8,'uint64');
                    fwrite(logger.fid,0,'uint64');
                else
                    %SubROIs
                    for i = 1:length(data.sROI)
                        ROI = data.sROI(i).ROI;
                        f=file;
                        [gg fname e] = fileparts(file.fname);
                        f.fname = sprintf('%s_%s%s',fname,data.sROI(i).name,e);
                        %Write individual sub-ROI metadata
                        if ~isempty(data.sROI(i).notes)
                            notes_name = sprintf('%s_%s%s',fname,data.sROI(i).name,'.txt');
                            a = fopen(notes_name,'w');
                            fprintf(a,'%s\r\n',f.fname);
                            fprintf(a,'%s\r\n\r\n',datestr(now));
                            for fi = 1:size(data.sROI(i).notes,1)
                                fprintf(a,'%s\r\n',data.sROI(i).notes(fi,:));
                            end
                            fclose(a);
                        end

                        ff = fullfile(file.pname,f.fname);
                        %Write FMF Header
                        logger(i).fid = fopen(ff,'w');
                        fwrite(logger(i).fid,1,'uint32');
                        fwrite(logger(i).fid,ROI(4),'uint32');
                        fwrite(logger(i).fid,ROI(3),'uint32');
                        fwrite(logger(i).fid,ROI(4)*ROI(3)+8,'uint64');
                        fwrite(logger(i).fid,0,'uint64');
                    end
                end
            case '.tif' %tiff stack for higher bit depth videos
                tagStruct.Artist = 'gVision';
                tagStruct.DateTime = datestr(now,30);
                tagStruct.ImageLength = data.res(2);
                tagStruct.ImageWidth = data.res(1);
                tagStruct.Compression = compID;

                info = imaqhwinfo(gui.vi);
                switch info.NativeDataType
                    case 'uint8'
                        tagStruct.BitsPerSample = 8;
                    case 'uint16'
                        tagStruct.BitsPerSample = 16;
                end
                if (~data.forceGray) & (get(gui.vi,'NumberOfBands')~=1)
                    tagStruct.SamplesPerPixel = 3;
                    tagStruct.Photometric = Tiff.Photometric.RGB;
                else
                    tagStruct.SamplesPerPixel = 1;
                    tagStruct.Photometric = Tiff.Photometric.MinIsBlack; 
                end

                tagStruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;

                if isempty(data.sROI)
                    %single ROI
                    ff = fullfile(file.pname,file.fname);
                    logger.ts = Tiff(ff,'w');
                    logger.ts.setTag(tagStruct);
                else
                    %subROIs
                    for i = 1:length(data.sROI)
                        ROI = data.sROI(i).ROI;
                        f=file;
                        [gg fname e] = fileparts(file.fname);
                        f.fname = sprintf('%s_%s%s',fname,data.sROI(i).name,e);
                        %Write individual sub-ROI metadata
                        if ~isempty(data.sROI(i).notes)
                            notes_name = sprintf('%s_%s%s',fname,data.sROI(i).name,'.txt');
                            a = fopen(notes_name,'w');
                            fprintf(a,'%s\r\n',f.fname);
                            fprintf(a,'%s\r\n\r\n',datestr(now));
                            for fi = 1:size(data.sROI(i).notes,1)
                                fprintf(a,'%s\r\n',data.sROI(i).notes(fi,:));
                            end
                            fclose(a);
                        end

                        ff = fullfile(file.pname,f.fname);                            
                        tagStruct.ImageLength = ROI(4);
                        tagStruct.ImageWidth = ROI(3);
                        logger(i).ts = Tiff(ff,'w');
                        logger(i).ts.setTag(tagStruct);
                    end
                end
                data.tagStruct = tagStruct;
        end
        data.logger = logger;

    end
end

methods (Static)
    % Itterate file name (max num of itterated files = 1000)
    % This function checks for sub-ROI extensions or other file names
    % in the current directory with the same root.  It returns a file
    % name with an interrated identifier
    function file = checkFileName(file)
        a=dir(file.pname);
        [gg fname ext] = fileparts(file.fname);
        maxind = 0;
        for i = 1:length(a)
            [gg f e] = fileparts(a(i).name);
            temps = sscanf(f,[fname,'%s']);
            if ~isempty(temps) && strcmpi(e,ext)
                %if file name is same extension
                % and CONTAINS name base
                %scan for file name w/ extension for sub-region
                n=sscanf(temps,['_%03d%*s',ext]);
                if isempty(n)
                    %see if there is no extension after the itterator
                    n = sscanf(temps,['_%03d',ext]);
                    if isempty(n)
                        %may be a sub-roi file with same root and no itterator
                        if maxind>1; continue; end
                        maxind = 1;
                        sprintf('%s_%03d%s',fname,1,ext);
                        file.fname = sprintf('%s_%03d%s',fname,1,ext);
                    else
                        maxind=n+1;
                        file.fname = sprintf('%s_%03d%s',fname,n+1,ext);
                    end
                else
                    %Already iterating, increment
                    maxind=n+1;
                    file.fname = sprintf('%s_%03d%s',fname,n+1,ext);
                end

            else
                %if file name matches extension and name base exactly
                if strcmpi(f,fname) && strcmpi(e,ext)
                    if maxind>1; continue; end
                    maxind = 1;
                    sprintf('%s_%03d%s',fname,1,ext);
                    file.fname = sprintf('%s_%03d%s',fname,1,ext);
                else
                    %does not match file at all
                end
            end
        end
        if maxind~=0; fprintf('\nNew File Name = %s\n\n',file.fname); end
    end
end
    
end

