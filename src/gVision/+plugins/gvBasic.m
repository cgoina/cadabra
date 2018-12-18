classdef gvBasic < plugins.gvPlugin
properties
    vidFileName
    vidFilePtr
    fig
    subROI
    hwTrigLabel
    hwTrigMenu
    trigButton
end

methods(Static)
    % This function must return the name of your plugin
    function [cName dispName description] = gvPlugName
        cName = 'gvBasic';
        dispName = 'Basic';
        description = 'Simple Acquisition Mode';
    end
end


methods 
    %Constructor which initializes your gVision plugin
    %Create gui elements here, setup other devices, etc
    function plug=gvBasic(fig)
        plug.fig = fig; %extract handles to all UI elements in the main GUI
        gui = get(plug.fig,'userdata');
        set(gui.StartStop,'callback',@(obj,event)StartStop(plug));
        
        triggerconfig(gui.vi,'immediate','none','none');
        
        tInfo = triggerinfo(gui.vi);
        str = {};
        for i=1:length(tInfo)
            str{i} = sprintf('%s - %s - %s',...
                tInfo(i).TriggerType, tInfo(i).TriggerCondition, tInfo(i).TriggerSource);
        end
        
        %Create Hardware Triggering Options Menu if any
        plug.hwTrigLabel = uicontrol(gui.RecPanel,'style','text','string','Triggering Options:','units','normalized','position',[.05 .55 .3 .06],'horizontalalignment','left',...
            'backgroundcolor',get(gui.RecPanel,'backgroundcolor'),'fontweight','bold');
        plug.hwTrigMenu = uicontrol(gui.RecPanel,'style','popupmenu',...
            'units','normalized','position',[0.38 0.54 0.57 0.08],'backgroundcolor','w',...
            'callback',@(obj,event)changeTrigMode(plug));
        set(plug.hwTrigMenu,'string',str,'userdata',tInfo)
        
        %For Manual Trigger
        plug.trigButton = uicontrol(gui.RecPanel,'style','pushbutton','string','Trigger','units','normalized','fontweight','bold',...
                    'position',[.6 .01 .19 .1],'backgroundcolor',[.7 .7 .8],'fontweight','bold','callback',@(obj,event)manualTrigger(plug),'visible','off');
        
    end
    %Destructor method to cleanup your gVision plugin when changing to
    %another mode.  This should essentially undo everything you did in
    %the constructor
    function delete(plug)
        gui = get(plug.fig,'userdata');
        
        delete(plug.hwTrigMenu)
        delete(plug.hwTrigLabel)
        delete(plug.trigButton)
        triggerconfig(gui.vi,'immediate','none','none');
        set(gui.StartStop,'position',[0.6 0.01 0.35 0.1],'string','Start');
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
                set(plug.hwTrigMenu,'enable','on')
                set(gui.RecStatus,'backgroundcolor',[1 .7 .7])
            case 1 %Button is down
                
                set(gui.StartStop,'string','Stop')
                
                %Open Video Files for Writing, returns a logdata struct
                %that contains information to be used by logging functions.
                % This struct contains things like file name, logging
                % stream pointers, flags about timestamping and forcegray,
                % etc.
                logdata = openVidFiles(plug);
                
                
                %Video Callbacks         
                gui.vi.framesacquiredfcn = @(obj,event)frameAcquired(plug);
                gui.vi.framesacquiredfcncount = 1;
                gui.vi.stopfcn = @(obj,event)wrapupVid(plug); %Built-in, calls frameAcquired to empty frames
                gui.vi.userdata = logdata;

                
                deactivateGvUI(plug);
                set(plug.hwTrigMenu,'enable','off')
                set(gui.RecStatus,'backgroundcolor',[.5 1 .5])
                start(gui.vi)
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
    
    %Switching to an available trigger mode in immediate
    function changeTrigMode(plug)
        gui = get(plug.fig,'userdata');
        val = get(plug.hwTrigMenu,'value');
        tInfo = get(plug.hwTrigMenu,'userdata');
        
        %Set the trigger mode
        triggerconfig(gui.vi,tInfo(val).TriggerType,...
            tInfo(val).TriggerCondition,tInfo(val).TriggerSource)
        
        if strcmpi(tInfo(val).TriggerType,'manual')
            %Turn on trigger button for manual
            set(plug.trigButton,'visible','on')
            set(gui.StartStop,'position',[0.8 0.01 0.15 0.1],'string','Start');
        else
            set(plug.trigButton,'visible','off')
            set(gui.StartStop,'position',[0.6 0.01 0.35 0.1],'string','Start');
        end
    end
    
    %Deliver a manual trigger to the video input
    function manualTrigger(plug)
        gui = get(plug.fig,'userdata');
        try
            trigger(gui.vi)
        catch
            disp('Manual Trigger Failed')
        end
    end
   
end



end