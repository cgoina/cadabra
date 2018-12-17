%% Write_clip
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
%% Writes out a detected action as animated GIF with time stamp imprint, 
%% receives information from 'ProcessFrameClustWin', part of 'plot_clips'

function [indx,indxg,b,ImageBuffer] = write_clip(InImage,dframe,start_frm,end_frm,position,timeindex,indx,indxg,b,ImageBuffer,addName,cluster,K,onefly,gen)

% AT START FRAME OF A BOUT, INITIALIZE BUFFER
% FIND ROI TO TIGHTLY CLIP THE LOCATION OF ACTION
if (K == start_frm(indx)),
    if onefly,
        ind_posx = position.xc1(indx,:);
        ind_posy = position.yc1(indx,:);
    else
        ind_posx = [position.xc1(indx,:) position.xc2(indx,:)];
        ind_posy = [position.yc1(indx,:) position.yc2(indx,:)];
    end
    if iscell(position.xc1),
        ind_posx = cell2mat(ind_posx); ind_posy = cell2mat(ind_posy);
    end
    % b = ROI (region of interest); this is the area where the action took
    % place, tightly clipped to focus the expert and to save memory
    b.minx = max(ceil(min(ind_posx)),31)-30; b.maxx = min(floor(max(ind_posx)),610)+30;
    b.miny = max(ceil(min(ind_posy)),31)-30; b.maxy = min(floor(max(ind_posy)),450)+30;
    sarr = size(imresize(zeros(b.maxy-b.miny+1,b.maxx-b.minx+1),.5));
    ImageBuffer = zeros(sarr(1),sarr(2),dframe(indx));
end
% BUFFER ALL FRAMES (ONLY ROI) BELONGING TO A BOUT
if (K >= start_frm(indx)),
    img = imresize(InImage(b.miny:b.maxy,b.minx:b.maxx),.5);
    ImageBuffer(:,:,K-start_frm(indx)+1) = img;
end
% IF THE LAST FRAME OF A BOUT IS READ, OUTPUT THE CLIP
if (K == end_frm(indx)),
    if isempty(gen), 
        clname = ['cluster_' num2str(cluster(indx),'%03d')];
    else
        clname = gen; 
    end
    if ispc, slash = '\'; else slash = '/'; end
    ind_tim = cell2mat(timeindex(indx,:));
    for i=1:dframe(indx),
        img = ImageBuffer(:,:,i);
        % TAKE IMAGE AND ADD A TIME STAMP
        timestmp = [num2str(fix(ind_tim(i)/60),'%02g') ':' num2str(mod(ind_tim(i),60),'%05.2f')];
        img = txtImage(img,timestmp);
        % OUTPUT IMAGE INTO ANIMATED GIF
        if (i==1),
            % Initialize file with first frame
            imwrite(img, [addName clname slash 'movie_' num2str(indxg(cluster(indx)),'%04d') '.gif'],'gif','LoopCount',Inf);
        else
            % Append following frames
            imwrite(img, [addName clname slash 'movie_' num2str(indxg(cluster(indx)),'%04d') '.gif'],'gif','WriteMode','append');
        end
    end
%     delete([addName 'cluster_' num2str(cluster(indx),'%03d') slash 'movie_' num2str(indxg(cluster(indx)),'%04d') '.gif']);
    % REMEMBER CLIP INDEX  - SAVED BY "plot_clips"
    % (in case of software crash, user stop, etc.)
    indxg(cluster(indx)) = indxg(cluster(indx)) + 1;
    indx = indx + 1;
end