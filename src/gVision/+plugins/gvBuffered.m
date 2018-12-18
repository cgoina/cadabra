classdef gvBuffered < plugins.gvPlugin
properties
    vidFileName
    vidFilePtr
    fig
    subROI
    trigButton
    preBuffEdit
    preBuffLabel
end

methods(Static)
    % This function must return the name of your plugin
    function [cName dispName description] = gvPlugName
        cName = 'gvBuffered';
        dispName = 'Buffered';
        description = 'Manual Trigger with PreBuffered Frames';
    end
end


methods 
    %Constructor which initializes your gVision plugin
    function plug=gvBuffered(fig)
        plug.fig = fig; %extract handles to all UI elements in the main GUI
        gui = get(plug.fig,'userdata');
        set(gui.StartStop,'callback',@(obj,event)StartStop(plug));
        triggerconfig(gui.vi,'immediate','none','none');
        
        gui.vi.tag = 'Buffer';
        gui.vi.framespertrigger = inf;
        gui.vi.triggerrepeat = 0;
        
        %Disable triggering options
        set(gui.FPT,'string','inf','userdata',inf,'enable','off')
        set(gui.TR,'string','0','userdata',0,'enable','off')
        
        % USER INTERFACE ELEMENTS
        %preBuffer size
        plug.preBuffLabel = uicontrol(gui.RecPanel,'style','text','string','preBuffer Size:','units','normalized','position',[.05 .39 .3 .06],'horizontalalignment','left',...
            'backgroundcolor',get(gui.RecPanel,'backgroundcolor'),'fontweight','bold');
        plug.preBuffEdit = uicontrol(gui.RecPanel,'style','edit','string','30','units','normalized','backgroundcolor','w',...
            'position',[0.32 0.38 0.2 0.08],'callback',@(obj,event)validateEntry(plug),'userdata',30);
    end
    %Destructor method to cleanup your gVision plugin when changing to
    %another mode.  This should essentially undo everything you did in
    %the constructor
    function delete(plug)
        gui = get(plug.fig,'userdata');
        
        delete(plug.preBuffLabel)
        delete(plug.preBuffEdit)
        
        try delete(plug.trigButton); end
        
        set(gui.FPT,'string','inf','userdata',inf,'enable','on')
        set(gui.TR,'string','0','userdata',0,'enable','on')
        
        set(gui.RecStatus,'backgroundcolor',[1 .7 .7])
    end

    %Callback for the start/stop button that initiates or stops
    %recording
    function StartStop(plug)
        gui = get(plug.fig,'userdata');
        switch get(gui.StartStop,'Value')
            case 0 %Button is up
                stop(gui.vi)
                set(gui.StartStop,'string','Start')
                delete(plug.trigButton)
                set(gui.StartStop,'position',[0.6 0.01 0.35 0.1],'string','Start');
                set(gui.RecStatus,'backgroundcolor',[1 .7 .7])
                                
                activateGvUI(plug);
            case 1 %Button is down
                                
                %Configure buttons for controlling buffering and
                %acquisition
                plug.trigButton = uicontrol(gui.RecPanel,'style','toggle','string','Trigger','units','normalized','fontweight','bold',...
                    'position',[.6 .01 .19 .1],'backgroundcolor',[.7 .7 .8],'fontweight','bold','callback',@(obj,event)startCapture(plug));
                set(gui.StartStop,'position',[0.8 0.01 0.15 0.1],'string','Halt');

                
                logdata = openVidFiles(plug);
               
                %Video Callbacks         
                gui.vi.framesacquiredfcn = @(obj,event)fillBuffer(plug);
                gui.vi.framesacquiredfcncount = 1;
                gui.vi.stopfcn = @(obj,event)closeVidFiles(plug);
                gui.vi.userdata = logdata;

                deactivateGvUI(plug);

                start(gui.vi)
        end
    end
    
    %Buffer the requested number of frames
    function fillBuffer(plug)
        gui = get(plug.fig,'userdata');
        %clear buffer down to max if it's filled beyond the target amount
        preTime = get(gui.RecStatus,'userdata');
        fa = gui.vi.framesavailable;
        nMax = get(plug.preBuffEdit,'userdata');
        data = gui.vi.userdata;
        data.postcount = 0;
       
        %drain buffer and display buffer status
        if fa > nMax
            [frames time]=getdata(gui.vi,fa-nMax); %Drain Buffer
            fps = 1/(time(1)-preTime);
            if fps <=0;  fps = data.fps;  end
            data.fps = fps;
            et = etime(clock,gui.vi.InitialTriggerTime);
            set(gui.RecStatus,'string',sprintf('nFrames: %u\neTime: %.2f s\nfps: %.1f Hz\nBuffer: %u',data.fc,et,data.fps,fa),'userdata',time(end))
        else 
            et = etime(clock,gui.vi.InitialTriggerTime);
            data.fps = gui.vi.framesacquired/et;
            set(gui.RecStatus,'string',sprintf('nFrames: %u\neTime: %.2f s\nfps: %.1f Hz\nBuffer: %u',data.fc,et,data.fps,fa),'userdata',et)
            
        end
        set(gui.vi,'userdata',data)
    end

    %Button callback to trigger capture of frames and whent to stop
    function startCapture(plug)
        gui = get(plug.fig,'userdata');
        switch get(plug.trigButton,'value')
            case 1 %Start Capture of buffer
                set(plug.trigButton,'string','Stop');
                set(gui.RecStatus,'backgroundcolor',[.5 1 .5])
                gui.vi.framesacquiredfcn = @(obj,event)drainBuffer(plug);
            case 0 %Stop Capture of buffer, finish current frames
                set(plug.trigButton,'string','Trigger');
                set(gui.RecStatus,'backgroundcolor',[1 .7 .7])
                gui.vi.framesacquiredfcn = '';
                fa = gui.vi.framesavailable;
                for i = 1:fa
                    drainBuffer(plug)
                end
                gui.vi.framesacquiredfcn = @(obj,event)fillBuffer(plug);
        end
        
    end
    
    %When trigger is pressed, log the frame history
    function drainBuffer(plug)
        gui = get(plug.fig,'userdata');
        data = get(gui.vi,'userdata');

        [frame time metadata] = getdata(gui.vi,1);
        
        if isempty(data.sROI)
            writeFullROI(plug,frame,time,metadata)
        else
            writeSubROI(plug,frame,time,metadata)
        end
        
    end
    
    %Validate Recording Property string entries
    function validateEntry(plug)
        val = round(str2double(get(plug.preBuffEdit,'string')));
        if isnan(val) | val<0
            set(plug.preBuffEdit,'string',num2str(get(plug.preBuffEdit,'userdata')))
        else
            set(plug.preBuffEdit,'string',sprintf('%d',val),'userdata',val)
        end
        
    end 
end
end