%% Detect_tussls
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
%% Detect tussling

function [fly_feat,tussl] = detect_tussls(fly_feat,params,FileN,ifile,tussl)

if (nargin > 2),
    if ~exist('tussl'),
        tussl = [];
    end
    if ~numel(tussl),
        tussl.ind = []; tussl.mov = []; tussl.t = []; tussl.tim = [];
        tussl.x = []; tussl.y = []; tussl.number = []; tussl.len = [];
    end
end

NFrms = length(fly_feat.frame);
% Time diffrence between two frames
dt = fly_feat(1).time(2:NFrms(1)) - fly_feat(1).time(1:NFrms(1)-1);
dt = median(dt);

% COMPUTE SOME PARAMETERS FOR TUSSLING DETECTION
fly_feat.tussl = zeros(1,NFrms);
diff1 = sqrt((fly_feat.obj1.pos_x(2:end)-fly_feat.obj2.pos_x(1:end-1)).^2 + ...
    (fly_feat.obj1.pos_y(2:end)-fly_feat.obj2.pos_y(1:end-1)).^2);
diff2 = sqrt((fly_feat.obj2.pos_x(2:end)-fly_feat.obj1.pos_x(1:end-1)).^2 + ...
    (fly_feat.obj2.pos_y(2:end)-fly_feat.obj1.pos_y(1:end-1)).^2);
dd = [1 abs(diff1-diff2)]; dd = conv2(dd, [1 1 1], 'same') / 3;
ori = 180/pi*(atan2((fly_feat.obj1.pos_y-fly_feat.obj2.pos_y), ...
    (fly_feat.obj1.pos_x-fly_feat.obj2.pos_x)));
dori1 = abs(ori-fly_feat.obj1.headdir);
ind = find(dori1 > 180); if numel(ind), dori1(ind) = dori1(ind) - 180; end
dori2 = abs(ori-fly_feat.obj2.headdir);
ind = find(dori2 > 180); if numel(ind), dori2(ind) = dori2(ind) - 180; end

% EXCLUDE ARENA BORDER AREA
if ~params.border,
    delta = params.border_width;
else
    delta = 0;
end
[ROI,scale] = load_ROI(FileN);
width.x = (ROI.cols(end)-ROI.cols(1))*scale.x;
width.y = (ROI.rows(end)-ROI.rows(1))*scale.y;
if (width.x > width.y), tmp = width.x; width.x = width.y; width.y = tmp; end
if min(fly_feat.obj1.pos_x)<(-1),
    minx = -width.x/2+delta; maxx = width.x/2-delta;
    miny = -width.y/2+delta; maxy = width.y/2-delta;
else
    minx = delta; maxx = width.x - delta;
    miny = delta; maxy = width.y - delta;
end

% DETECT TUSSLING FRAMES
ind = find((fly_feat.obj1.acc > 80) & (fly_feat.obj2.acc > 80) & ...
    (fly_feat.obj1.vel > 10) & (fly_feat.obj2.vel > 10) & ...
    fly_feat.distc < 1.7 & dd < 1. & ...
    (dori1>60 & dori1<120) & (dori2>60 & dori2<120) & ...
    (fly_feat.obj1.pos_change < 3) & (fly_feat.obj2.pos_change < 3) & ...
    (fly_feat.obj1.pos_x > minx) & (fly_feat.obj1.pos_x < maxx) & ...
    (fly_feat.obj2.pos_x > minx) & (fly_feat.obj2.pos_x < maxx) & ...
    (fly_feat.obj1.pos_y > miny) & (fly_feat.obj1.pos_y < maxy) & ...
    (fly_feat.obj2.pos_y > miny) & (fly_feat.obj2.pos_y < maxy));

% DETECT BOUTS OF TUSSLING AND SETUP DATA  STRUCTURE
ntussl = numel(ind);
if ntussl,
    ntussl = 0;
    if ind(1) == 1, ind(1) = 2; end
    % Bouts with gaps ? 'max_gap' are considered as one bout
    % Minimum bout length = min_len
    % A minimum acceleration is expected, corresponding to observed thrusts
    thres = 0.3; max_gap = round(0.3/dt); min_len = round(0.3/dt);
    for i=1:numel(ind),
        bout = double(fly_feat.obj1.acc > max(fly_feat.obj1.acc(ind(i)-1:ind(i)+1))*thres & ...
            fly_feat.obj2.acc > max(fly_feat.obj2.acc(ind(i)-1:ind(i)+1))*thres);
        bout = conv2(bout,ones(1,max_gap+1),'same'); bout = bout > 0;
        bout = bwlabel(bout);
        val = bout(ind(i)-1:ind(i)+1); val = val(logical(val > 0));
        if val,
            in = find(bout == val(1));
            if numel(in) > min_len,
                ntussl = ntussl + 1;
                fly_feat.tussl(in) = i;
            end
        end
    end
end

if (nargin > 3),
    if ntussl,
        clear bout; bout.len = []; bout.tim = []; bout.ind = [];
        for i=1:ntussl,
            in = find(fly_feat.tussl == i);
            if numel(in) > min_len,
                bout.tim = [bout.tim fly_feat.time(in(1))];
                bout.len = [bout.len length(in)];
                bout.ind = [bout.ind in];
            end
        end
        ntussl = numel(bout.len);

        tussl.t = [tussl.t bout.tim];
        tussl.tim = [tussl.tim fly_feat.time(bout.ind)];
        tussl.len = [tussl.len bout.len];
        tussl.number = [tussl.number ntussl];
        tussl.ind = [tussl.ind fly_feat.lab(bout.ind)];
        tussl.x = [tussl.x fly_feat.obj1.pos_x(bout.ind)]; tussl.y = [tussl.y fly_feat.obj1.pos_y(bout.ind)];
    else
        tussl.number = [tussl.number 0];
    end
    tussl.mov(ifile) = ifile;
end

end