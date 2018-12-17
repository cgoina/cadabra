function OpenAvi;
global strVideoFExt;
global strInVideoFName;
global strInVideoPath;
global findex;
global filechg;
global InPause;

[strInVideoFName, strInVideoPath, findex] = uigetfile({'*.avi','AVI-file (*.avi)';'*.fmf','fmf-file (*.fmf)'}, ...
                                                        'Open movie file','MultiSelect','off');
if findex, 
    strInVideoPath = strInVideoPath(1:end-1);
    strVideoFExt = strInVideoFName(end-2:end);
    strInVideoFName = strInVideoFName(1:end-4);
    InPause = 0;
    filechg = 1; 
end
return
