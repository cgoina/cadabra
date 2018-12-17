%% Wing_ext
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
%% Detect wing extension

% WING EXTENSION DETECTION
function wing = wing_ext(fly_feat,params,FileN,min_b_len,max_gap,wing)

% WING EXTENSION RANGE LIMITS
% min_phi = 60; max_phi = 90; min_wingl = 1.15; max_wingl = 2.5; min_flylen = 1.2;
min_phi = 60; max_phi = 88; min_wingl = 1.15; max_wingl = 1.9; min_flylen = 1.2;


% Initialize "wing" structure, in case there is no structure yet
if ~exist('wing'),
    wing.ext.r = []; wing.ext.l = []; wing.ext.b = [];
    wing.fli.r = []; wing.fli.l = []; wing.fli.b = [];
end
% Time difference between two successing frames
NFrms = length(fly_feat.frame);
dt = median(fly_feat(1).time(2:NFrms(1)) - fly_feat(1).time(1:NFrms(1)-1));

% EXCLUDE BORDER AREA OF ARENA?
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

min_b_len1.long = min_b_len.long/dt; min_b_len1.short = min_b_len.short/dt;
max_gap1.long = max_gap.long/dt; max_gap1.short = max_gap.short/dt;

if max(fly_feat.obj1.FLength) < .1, fly_feat.obj1.FLength = fly_feat.obj1.length*2; end
if max(fly_feat.obj2.FLength) < .1, fly_feat.obj2.FLength = fly_feat.obj2.length*2; end

% DETECT FRAMES WITH PROBABLE WING EXTENSIONS RIGHT, LEFT, FLY 1 & 2
if params.radius > 0
    r1 = sqrt(fly_feat.obj1.pos_x.^2 + fly_feat.obj1.pos_y.^2);
    r2 = sqrt(fly_feat.obj2.pos_x.^2 + fly_feat.obj2.pos_y.^2);
    ind1r = (fly_feat.obj1.phir > min_phi) & (fly_feat.obj1.phir < max_phi) & ...
        (fly_feat.obj1.winglr > min_wingl) & (fly_feat.obj1.winglr < max_wingl) & ...
        (fly_feat.obj1.winglr < fly_feat.obj1.FLength) & (fly_feat.obj1.length > min_flylen) & ...
        (r1 < params.radius);
    ind2r = (fly_feat.obj2.phir > min_phi) & (fly_feat.obj2.phir < max_phi) & ...
        (fly_feat.obj2.winglr > min_wingl) & (fly_feat.obj2.winglr < max_wingl) & ...
        (fly_feat.obj2.winglr < fly_feat.obj2.FLength) & (fly_feat.obj2.length > min_flylen) & ...
        (r2 < params.radius);
    ind1l = (fly_feat.obj1.phil > min_phi) & (fly_feat.obj1.phil < max_phi) & ...
        (fly_feat.obj1.wingll > min_wingl) & (fly_feat.obj1.wingll < max_wingl) & ...
        (fly_feat.obj1.wingll < fly_feat.obj1.FLength) & (fly_feat.obj1.length > min_flylen) & ...
        (r1 < params.radius);
    ind2l = (fly_feat.obj2.phil > min_phi) & (fly_feat.obj2.phil < max_phi) & ...
        (fly_feat.obj2.wingll > min_wingl) & (fly_feat.obj2.wingll < max_wingl) & ...
        (fly_feat.obj2.wingll < fly_feat.obj2.FLength) & (fly_feat.obj2.length > min_flylen) & ...
        (r2 < params.radius);    
else
    ind1r = (fly_feat.obj1.phir > min_phi) & (fly_feat.obj1.phir < max_phi) & ...
        (fly_feat.obj1.winglr > min_wingl) & (fly_feat.obj1.winglr < max_wingl) & ...
        (fly_feat.obj1.winglr < fly_feat.obj1.FLength) & (fly_feat.obj1.length > min_flylen) & ...
        (fly_feat.obj1.pos_x > minx) & (fly_feat.obj1.pos_x < maxx) & ...
        (fly_feat.obj2.pos_x > minx) & (fly_feat.obj2.pos_x < maxx) & ...
        (fly_feat.obj1.pos_y > miny) & (fly_feat.obj1.pos_y < maxy) & ...
        (fly_feat.obj2.pos_y > miny) & (fly_feat.obj2.pos_y < maxy);
    ind2r = (fly_feat.obj2.phir > min_phi) & (fly_feat.obj2.phir < max_phi) & ...
        (fly_feat.obj2.winglr > min_wingl) & (fly_feat.obj2.winglr < max_wingl) & ...
        (fly_feat.obj2.winglr < fly_feat.obj2.FLength) & (fly_feat.obj2.length > min_flylen) & ...
        (fly_feat.obj1.pos_x > minx) & (fly_feat.obj1.pos_x < maxx) & ...
        (fly_feat.obj2.pos_x > minx) & (fly_feat.obj2.pos_x < maxx) & ...
        (fly_feat.obj1.pos_y > miny) & (fly_feat.obj1.pos_y < maxy) & ...
        (fly_feat.obj2.pos_y > miny) & (fly_feat.obj2.pos_y < maxy);
    ind1l = (fly_feat.obj1.phil > min_phi) & (fly_feat.obj1.phil < max_phi) & ...
        (fly_feat.obj1.wingll > min_wingl) & (fly_feat.obj1.wingll < max_wingl) & ...
        (fly_feat.obj1.wingll < fly_feat.obj1.FLength) & (fly_feat.obj1.length > min_flylen) & ...
        (fly_feat.obj1.pos_x > minx) & (fly_feat.obj1.pos_x < maxx) & ...
        (fly_feat.obj2.pos_x > minx) & (fly_feat.obj2.pos_x < maxx) & ...
        (fly_feat.obj1.pos_y > miny) & (fly_feat.obj1.pos_y < maxy) & ...
        (fly_feat.obj2.pos_y > miny) & (fly_feat.obj2.pos_y < maxy);
    ind2l = (fly_feat.obj2.phil > min_phi) & (fly_feat.obj2.phil < max_phi) & ...
        (fly_feat.obj2.wingll > min_wingl) & (fly_feat.obj2.wingll < max_wingl) & ...
        (fly_feat.obj2.wingll < fly_feat.obj2.FLength) & (fly_feat.obj2.length > min_flylen) & ...
        (fly_feat.obj1.pos_x > minx) & (fly_feat.obj1.pos_x < maxx) & ...
        (fly_feat.obj2.pos_x > minx) & (fly_feat.obj2.pos_x < maxx) & ...
        (fly_feat.obj1.pos_y > miny) & (fly_feat.obj1.pos_y < maxy) & ...
        (fly_feat.obj2.pos_y > miny) & (fly_feat.obj2.pos_y < maxy);
end

% FILLING 1-FRAME GAPS AND FIND BOTH WINGS EXTENDED (ind1, ind2)
smo = 3;
ind1r = double(bwmorph(ind1r,'clean'));
ind2r = double(bwmorph(ind2r,'clean'));
ind1l = double(bwmorph(ind1l,'clean'));
ind2l = double(bwmorph(ind2l,'clean'));
ind1 = conv2(double(ind1r & ind1l),ones(1,smo+1),'same')>0;
ind2 = conv2(double(ind2r & ind2l),ones(1,smo+1),'same')>0;

% DETECT BOUTS OF RIGHT WING EXTENSIONS
[wing.ext.r,i1,i2] = detect_bouts(fly_feat,ind1r,ind2r,min_b_len1.long,max_gap1.long,wing.ext.r);
% DETECT BOUTS OF RIGHT WING FLICKS (ALL WING EXT. THAT WERE TOO SHORT)
ind1r(i1) = 0; ind2r(i2) = 0;
wing.fli.r = detect_bouts(fly_feat,ind1r,ind2r,min_b_len1.short,max_gap1.short,wing.fli.r);

% DETECT BOUTS OF LEFT WING EXTENSIONS
[wing.ext.l,i1,i2] = detect_bouts(fly_feat,ind1l,ind2l,min_b_len1.long,max_gap1.long,wing.ext.l);
% DETECT BOUTS OF LEFT WING FLICKS (ALL WING EXT. THAT WERE TOO SHORT)
ind1l(i1) = 0; ind2l(i2) = 0;
wing.fli.l = detect_bouts(fly_feat,ind1l,ind2l,min_b_len1.short,max_gap1.short,wing.fli.l);

% DETECT BOUTS OF BOTH WING EXTENSIONS
[wing.ext.b,i1,i2] = detect_bouts(fly_feat,ind1,ind2,min_b_len1.long,0,wing.ext.b);
% DETECT BOUTS OF BOTH WING FLICKS (ALL WING EXT. THAT WERE TOO SHORT)
ind1(i1) = 0; ind2(i2) = 0;
wing.fli.b = detect_bouts(fly_feat,ind1,ind2,min_b_len1.short,0,wing.fli.b);

end