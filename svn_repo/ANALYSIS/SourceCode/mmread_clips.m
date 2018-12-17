%% Mmread_clips
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
%% Part of plot_clips to provide clips of detected actions

function mmread_clips(FName,str_frm,lab,pos1,timestamp,window,add_name,igen,ngens,cl,ifile,nfiles,params)
global ImageBuffer;
global dframe;
global indx;
global start_frm;
global end_frm;
global position;
global timeindex;
global addName;
global cluster;
global onefly;
global gen;
global h1;

dframe = window;
addName = add_name;
cluster = lab;
position = pos1;
timeindex = timestamp;
ImageBuffer = [];
onefly = params.oneobj;
gen = cl;
indx = 1;

% COMPUTE START/END FRAME INDICIES
start_frm = str_frm;
end_frm = str_frm + dframe-1;

% REMOVE OVERLAPPING CLIPS
% Clips must be separated by at least 1 frame, because
% output routine passed only once trough each movie
old = end_frm(1); ind = 1;
for i=2:length(start_frm),
    if start_frm(i)-old > 0,
        ind = [ind  i];
        if old < end_frm(i), old = end_frm(i); end
    end
end
start_frm = start_frm(ind); 
end_frm = end_frm(ind);
dframe = dframe(ind);
cluster = cluster(ind);
position.xc1 = position.xc1(ind,:);
position.yc1 = position.yc1(ind,:);
position.xc2 = position.xc2(ind,:);
position.yc2 = position.yc2(ind,:);
timeindex = timeindex(ind,:);
in = strfind(FName,'/'); if isempty(in), in = strfind(FName,'\'); end
in = in(end);
filen = FName(in+1:end); in = strfind(filen,'_'); filen(in) = '-';

% Show a waitbar during movie clip output
h1 = waitbar(0,['Genotype ' num2str(igen) '/' num2str(ngens) ', Movie ' num2str(ifile) '/' num2str(nfiles) ': ' filen ' - Frames: ' num2str(min(start_frm)) '-' num2str(max(end_frm))]);

% if ispc, % Removed, "mmread" is supposed to run universally on both platforms
    ProcessActionClips;
    % Clear previous instance of video grabber
    if ispc,
        clear('mexDDGrab.mexw32');
    else
        clear('FFGrab.mexglx');
        clear('FFGrab.mexmaci');
        clear('FFGrab.mexmaci64');
    end
    % CALL VIDEO GRABBER (C++), 
    % WHICH IN TURN CALLS "ProcessActionClips"
    try
        %mmread([FName '.wmv'],min(start_frm):max(end_frm),[],false,true,'ProcessActionClips',false);
        mmread([FName '.avi'],min(start_frm):max(end_frm),[],false,true,'ProcessActionClips',false); % change to read avi EH100423
    catch
        err = lasterror;
        if ~strcmp(err.message(end-14:end),'STOP PROCESSING')
            rethrow(err); % if it was any other type of error, pass it on so that we can see it
        end
    end
    % Clear previous instance of video grabber; not sure if this helps here
    if ispc,
        clear('mexDDGrab.mexw32'); % Clear previous instance of video grabber
    else
        clear('FFGrab.mexglx');
    end
% else
%     strInFName = [FName '.wmv'];
%     strOutFName = strInFName;
%     strProcFrmFcn = 'ProcessActionClips';
%     bolMovOrFrame = 0;                                      % Get an AVI Movie out
%     bolMakeGray = 0;                                        % Leave RGB as is.
%     intQuality = 10000;                                     % Highest Quality
%     intFps = 29.97;                                         % 30 Frame per Second
%     ProcessAVI(strInFName, strOutFName, strProcFrmFcn, bolMovOrFrame, ...
%         bolMakeGray, min(start_frm), max(end_frm)-min(start_frm)+1, intQuality, intFps);
% end

% Close waitbar
close(h1);

end