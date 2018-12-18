classdef gvTimeLapse < plugins.gvPlugin
    % Carries out time lapse recording
    %  The user specifies a capture interval, and a software timer drives
    %  frame capture at the given interval
properties
    vidFileName
    vidFilePtr
    fig
    subROI
    trigButton
    timeLabel
    timeEdit
    fpsVisFlag
end

methods(Static)
    % This function must return the name of your plugin
    function [cName dispName description] = gvPlugName
        cName = 'gvTimeLapse';
        dispName = 'Time Lapse';
        description = 'Long Interval Capture';
    end
end


methods 
    %Constructor which initializes your gVision plugin
    function plug=gvTimeLapse(fig)
        plug.fig = fig; %extract handles to all UI elements in the main GUI
        gui = get(plug.fig,'userdata');
        set(gui.StartStop,'callback',@(obj,event)StartStop(plug));
        triggerconfig(gui.vi,'manual','none','none');  
        
        %TimeLapse Interval - Should be moved and created in the time-lapse plugin
        plug.timeLabel = uicontrol(gui.RecPanel,'style','text','string','Cap Interval (s):','units','normalized','position',[.05 .49 .3 .06],'horizontalalignment','left',...
            'backgroundcolor',get(gui.RecPanel,'backgroundcolor'),'fontweight','bold');
        plug.timeEdit = uicontrol(gui.RecPanel,'style','edit','string','1','units','normalized','backgroundcolor','w',...
            'position',[0.32 0.48 0.2 0.08],'callback',@(obj,event)validateEntry(plug),'userdata',1);
        
        %Disable triggering options
        set(gui.FPT,'string','1','userdata',1,'enable','on')
        set(gui.TR,'string','inf','userdata',inf,'enable','on')
        gui.vi.framespertrigger = 1;
        gui.vi.triggerrepeat = inf;
        if strcmpi(get(gui.FPS,'visible'),'on')
            plug.fpsVisFlag = 1;
        else
            set([gui.FPStxt gui.FPS],'visible','on')
            plug.fpsVisFlag = 0;
        end
    end
    %Destructor method to cleanup your gVision plugin when changing to
    %another mode.  This should essentially undo everything you did in
    %the constructor
    function delete(plug)
        gui = get(plug.fig,'userdata');
        
        delete(plug.timeLabel)
        delete(plug.timeEdit)
        
        set(gui.StartStop,'callback','')
        if ~plug.fpsVisFlag %Conditionally turn off FPS display if it was off to start with
            set([gui.FPStxt gui.FPS],'visible','off')
        end
        
        %Reset triggering options
        set(gui.FPT,'string','inf','userdata',inf,'enable','on')
        set(gui.TR,'string','0','userdata',0,'enable','on')
        gui.vi.framespertrigger = inf;
        gui.vi.triggerrepeat = 0;
    end
    
    %Callback for the start/stop button that initiates or stops
    %recording
    function StartStop(plug)
        gui = get(plug.fig,'userdata');
        switch get(gui.StartStop,'Value')
            case 0 %Button is up
                stop(gui.vi)
                set(gui.StartStop,'string','Start')
                activateGvUI(plug);
                set(plug.timeEdit,'enable','on');
                set(gui.RecStatus,'backgroundcolor',[1 .7 .7])
                
            case 1 %Button is down
                               
                %Create Button to Execute Manual Trigger
                set(gui.StartStop,'string','Stop');
                
                %File is created at 15fps regardless of capture rate
                logdata = openVidFiles(plug); 
                
                %Video Callbacks         
                gui.vi.framesacquiredfcn = @(obj,event)frameAcquired(plug);
                gui.vi.framesacquiredfcncount = 1;
                gui.vi.stopfcn = @(obj,event)wrapupVid(plug);
                gui.vi.userdata = logdata;

                deactivateGvUI(plug);
                set(plug.timeEdit,'enable','off');
                set(gui.RecStatus,'backgroundcolor',[.5 1 .5])
                start(gui.vi)
                
                %Use timerfcn built into imaq object
                gui.vi.timerfcn = @(obj,event)trigger(gui.vi);
                gui.vi.timerperiod = get(plug.timeEdit,'userdata');
        end
    end
    
        
    %Calls inherited video writing functions
    function frameAcquired(plug)
        gui = get(plug.fig,'userdata');
        [frame time metadata] = getdata(gui.vi,1);
        
        data = get(gui.vi,'userdata');
        if isempty(data.sROI)
            writeFullROI(plug,frame,time,metadata);
        else
            writeSubROI(plug,frame,time,metadata);
        end
        
        %Append nTriggers repeated onto end of status message
        set(gui.RecStatus,'string',char(get(gui.RecStatus,'string'),sprintf('nTriggers: %d',gui.vi.TriggersExecuted)))
        
    end
    
    %Validate Recording Property string entries
    function validateEntry(plug)
        val = round(str2double(get(plug.timeEdit,'string')));
        if isnan(val) | val<0
            set(plug.timeEdit,'string',num2str(get(plug.timeEdit,'userdata')))
        else
            set(plug.timeEdit,'string',sprintf('%d',val),'userdata',val)
        end
        
    end 

end



end