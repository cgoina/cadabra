%% Auto_cal_circ
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
% Automatic detection of circular arenas; manually adding missed arenas;
% fully manual arena detemination in case the automatic detection fails

function coord = auto_cal_circ(mean_image,image1,diam_pix,Panels,FigureHandle)
global params;
img = uint8(mean_image*255);

%% DETECT CIRCULAR ARENAS
[accum, circen, cirrad] = CircularHough_Grd(img, ...
                            [diam_pix/2/1.2 diam_pix/2*1.2], 1, 50, .3);

% REMOVE OVERLAPPING CIRCLE
 %******************
 %JL05062009 CircularHough_Grd found and output some circular arenas twice, Add this
 %statement to temporarily fix the bugs.
circen = unique(circen, 'rows');
%*******************
 %JL06042009 Comment out the following statements because Eric is using a new 
 % moive capture software. Some movie files only have one chamber. The following
 % codes assumes there are always more than one circles.
% cendist = tril(sqrt(dist2(circen,circen)));
% inds = find(cendist > 0 & cendist <= diam_pix/3);
% inds = unique([mod(inds,length(circen)) ; find(cirrad <= 0)]);
% inds = setdiff((1:length(circen))',inds);
% circen = circen(inds,:);
% cirrad = cirrad(inds);

% COLLECT COORDINATES
coord = [];
for ii=1:size(circen, 1)
    % Flexible diameter
%     coord(ii).x = [circen(ii,1)-cirrad(ii) circen(ii,1)+cirrad(ii)];
%     coord(ii).y = [circen(ii,2)-cirrad(ii) circen(ii,2)+cirrad(ii)];
    % Constant diameter
    coord(ii).x = [circen(ii,1)-diam_pix/2 circen(ii,1)+diam_pix/2];
    coord(ii).y = [circen(ii,2)-diam_pix/2 circen(ii,2)+diam_pix/2];
end

% VISUAL VERIFICATION
axes(Panels.HMenu);
cla(Panels.HMenu);
image(Panels.BlankMenuPanel);
set(Panels.HMenu,'XTick',[],'YTick',[]);

figure(FigureHandle);
axes(Panels.HImage);
image(image1);
set(Panels.HImage,'XTick',[],'YTick',[]);
drawnow; hold on;
plot(circen(:,1), circen(:,2), 'r+');
for ii = 1 : size(circen, 1)
%     rectangle('Position',[coord(ii).x(1), coord(ii).y(1), 2*cirrad(ii), 2*cirrad(ii)],...
%         'Curvature', [1,1], 'edgecolor', 'b', 'linewidth', 1.5);
    rectangle('Position',[coord(ii).x(1), coord(ii).y(1), diam_pix, diam_pix],...
        'Curvature', [1,1], 'edgecolor', 'b', 'linewidth', 1.5);
end

%% DID THE AUTO-CALIBRATION FIND ALL ARENAS?
if isstruct(coord),
    %JL10012009 show dialog depend on batch processing or not
    if params.batchprocess == 1
        button = 'Yes';
    else
        button = questdlg('Correct?','ROIs','Yes','Some Missing','No','Yes');
    end

else
    button = 'No';
end
            
% AUTO-CALIBRATION INACCURATE?
if ~strcmp(button,'Yes') || ~numel(button),
    axes(Panels.HMenu);
    cla(Panels.HMenu);
    image(Panels.BlankMenuPanel);
    set(Panels.HMenu,'XTick',[],'YTick',[]);
    h = text(round(Panels.PanelDimX/2),105,'Click into');
    set(h,'HorizontalAlignment','center','FontSize',8,'color','w');
    h = text(round(Panels.PanelDimX/2),120,'Arena Center');
    set(h,'HorizontalAlignment','center','FontSize',8,'color','w');
    figure(FigureHandle); axes(Panels.HImage); cla(Panels.HImage);
    if strcmp(button,'No'), 
        coord = [];
        figure(FigureHandle);
        axes(Panels.HImage);
        image(image1);
        set(Panels.HImage,'XTick',[],'YTick',[]);
        drawnow; hold on;
    end
    button = [];
    while ~strcmp(button,'Finish'),
        [coord,circen,button] = getroi_circ(image1,diam_pix,Panels,FigureHandle,coord,circen);
    end
end

% SORT ARENAS
[s,is] = sort(circen(:,2));
tmp = coord; tmpc = circen;
for ii=1:length(coord),
    coord(ii).x = tmp(is(ii)).x;
    coord(ii).y = tmp(is(ii)).y;
    circen(ii,:) = tmpc(is(ii),:);
end
diff = circen(2:end,2) - circen(1:end-1,2);
ind = [0 find(diff > mean(diff))' length(coord)];
tmp = coord; tmpc = circen;
for i=1:numel(ind)-1,
    [s,is] = sort(circen(ind(i)+1:ind(i+1),1));
    for ii=1:length(s),
        coord(ind(i)+ii).x = tmp(ind(i)+is(ii)).x;
        coord(ind(i)+ii).y = tmp(ind(i)+is(ii)).y;
        circen(ind(i)+ii,:) = tmpc(ind(i)+is(ii),:);
    end
end


end

function [coord,circen,button] = getroi_circ(img, diam_pix, Panels, FigureHandle, coord, circen) %#ok<INUSL>

%%
% If the coordinate structure is empty, initialize it with 
% empty arrays, and set the length to zero. Otherwise, 
% update the length to reflect the number of coordinates stored

if isempty(coord),
    circen = [];
    coord.x = [];
    coord.y = [];
    ncoord = 0;
else
    ncoord = length(coord);
end

%%
% Execute the following loop until user is satisfied with
% the circle created on the screen.

buttonid = 0;
while ~buttonid,
    image(img); 
    set(Panels.HImage,'XTick',[],'YTick',[]); 
    drawnow; hold on;
    for ii=1:ncoord,
        plot(sum(coord(ii).x)/2, sum(coord(ii).y)/2, 'r+');
        rectangle('Position',[coord(ii).x(1), coord(ii).y(1), diam_pix, diam_pix],...
            'Curvature', [1,1], 'edgecolor', 'b', 'linewidth', 1.5);
    end
    
    [x,y] = ginput(1);
    plot(x, y, 'r+');
    rectangle('Position',[x-diam_pix/2, y-diam_pix/2, diam_pix, diam_pix],...
        'Curvature', [1,1], 'edgecolor', 'b', 'linewidth', 1.5);
    button = questdlg('Correct?','ROI','Finish','Yes','No','Finish');
    if strcmp(button,'Yes') || strcmp(button,'Finish'), 
        buttonid = 1; 
    end
end

%%
% The rectangle corners are then added into the 
% coordinate array as a new region of interest 

circen = [circen ; [x y]];
coord(ncoord+1).x = [ x - diam_pix/2, x + diam_pix/2];
coord(ncoord+1).y = [ y - diam_pix/2, y + diam_pix/2];

end
