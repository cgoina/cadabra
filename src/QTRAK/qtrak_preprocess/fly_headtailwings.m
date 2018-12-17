%% Fly_headtailwings
% Copyright (C) 2008 Heiko Dankert, California Institute of Technology

% This file is part of QTRAK
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
% * Written by Heiko Dankert 2006/2007
% * Documentation by Edwin Soedarmadji 10/2008
%
%% Head / Tail and Wing Determination
% After we determine the pixel indices, areas, representative
% ellipses (from which the direction of the flies are determined) for
% the two flies, we still have to resolve one more ambiguity: the
% direction to which the head is pointing. To do this, we use the fact
% that the head contains more dark pixels compared to the tail area.
%
% First, we use the grayscale values from img2 instead of img (1a).
% The values exhibit more contrast between the head and the tail area,
% which is useful for distinguishing the two (but at the same time
% less useful for identifying the overall extent of the flies).
% From the previous measurement, the changes in position and direction
% are calculated.

function obj1 = fly_headtailwings(obj1,obj_1,img2,f01,chamber,params,ind_head,bool_fly2)

if bool_fly2,
    ind_head1 = ind_head;
else
    ind_head1 = 1;
end

% ...................(1a)
img_f01 = img2(f01);
obj1.mean = mean(img_f01);

% ...................(1b)
mv1.x = obj1.xc - obj_1.xc;
mv1.y = obj1.yc - obj_1.yc;
obj1.mv_d = sqrt(mv1.x^2+mv1.y^2);
obj1.mv_phi = atan2(mv1.y,mv1.x);

%%
% Next, based on their location with respect to the ellipses'
% centroid, the fly pixels are split into two categories (roughly
% front or back). Note that at this point, it is still not clear
% whether the ellipse's front corresponds to the fly's head (2a).
% For each pixel, the angle with respect to the ellipse's
% orientation is calculated. The two pixel groups are stored in
% _ar_a_ and _ar_b_.
%
% Brightness histograms are then computed for each of the two
% pixel groups (2b). Each histogram has eleven bins, and the
% populations of the first four bins of each group are averaged
% and stored into _a_ and _b_.

% ...................(2a)
y1 = chamber.y_arr(f01) - obj1.eyc;
x1 = chamber.x_arr(f01) - obj1.exc;
aty1x1 = atan2(y1,x1);
th = aty1x1 - obj1.phi + params.hpi;
i = find(th > pi);
th(i) = th(i) - params.twopi;
i = find(th < -pi);
th(i) = th(i) + params.twopi;
ar_a = numel(find(th<0));
ar_b = numel(find(th>0));

% ...................(2b)
a = img_f01;
b = a;
a = fhist(a((th < -pi/4) & (th > -3*pi/4)),0:0.1:1)';
b = fhist(b((th > pi/4) & (th < 3*pi/4)),0:0.1:1)';
a = mean(a(1:4));
b = mean(b(1:4));

%%
% The main concept here is that the front part of a fly contains
% both more pixels (as indicated by _ar_a_ > _ar_b_) and more
% darker pixels (as indicated by _a_ > _b_ ). If these two
% conditions are simultaneously satisfied, then the front of the
% ellipse coincides with the fly's head. Otherwise, the direction
% is flipped (2c).

% ...................(2c)
if (a > b) || (ar_a > ar_b),
    obj1.head = real(obj1.phi);
else
    obj1.head = real(obj1.phi) + pi;
    if (obj1.head > pi),
        obj1.head = obj1.head - params.twopi;
    end
end

%%
% Finally, the change in the fly' head direction is calculated.
% The change is calculated with respect to the direction recorded
% in the previous frame. In addition, the angles are kept to
% within -Pi and Pi. If the change in the head direction is so
% abrupt (approaching Pi), it also indicates that the fly's
% direction is reversed.

if obj_1.head > (-99) && ind_head1,
    if ind_head,% || params.dist2 < 3.2,
        dhead1 = abs(obj1.head - obj_1.head);
        if (dhead1 > pi),
            dhead1 = dhead1 - params.twopi;
        end
        if ((dhead1 > 14.5*pi/18) || (dhead1 < -14.5*pi/18)),
            obj1.head = obj1.head + pi;
        end
        if (obj1.head > pi),
            obj1.head = obj1.head - params.twopi;
        end
    end
    dhead1 = abs(obj1.head - obj1.mv_phi);
    if (dhead1 > pi),
        dhead1 = dhead1 - params.twopi;
    end
    if (obj1.mv_d > 0.8) && ((dhead1 > 14.5*pi/18) || ...
            (dhead1 < -14.5*pi/18)),
        obj1.head = obj1.head + pi;
    end
    if (obj1.head > pi),
        obj1.head = obj1.head - params.twopi;
    end
end

%% B.10 Wing Angle Determination
% One of the last few tasks to do is to determine the wing angle
% relative to the body direction. To accomplish this, the arctan
% array used in (2a) is reused for a different purpose. This time,
% it is used to calculate the angle for each pixels, which is bound
% to a range between -Pi and Pi (3a).
%
% Next, all pixels within a certain angle to the ``left'' of the
% body axis are grouped together, and the distance of these pixels
% to the ellipse center is measured (3b). Finally, the pixel _ind_
% that lies the farthest from the ellipse center (in the prescribed
% angular range) is identified (3c). The left wing angle is then
% easily calculated from _ind_. This process is repeated in (3d)
% to obtain the right wing angle.

% ...................(3a)
if (obj1.head ~= obj1.phi),
    th = aty1x1 - obj1.head + params.hpi;
    i = find(th > pi);
    th(i) = th(i) - params.twopi;
    i = find(th < -pi);
    th(i) = th(i) + params.twopi;
end

% ...................(3b)
indxy = find((th > -params.hpi*7/9) & (th < 0));
obj1.l = sqrt(y1(indxy).^2 + x1(indxy).^2);

% ...................(3c)
ind = find(obj1.l == max(max(obj1.l)));
if ind,
    %% Method 1 (chosen) - Use wing pixel 
    %% (Pixel with farthest distance to fly body center) to 
    %% Compute Wing Length & Angle
    ind = ind(1);
    obj1.phil = params.hpi-abs(th(indxy(ind)));
    obj1.l = obj1.l(ind);
    %% Method 2 - weighted by wing-pixel distance to the fourth power
%     obj1.phil = params.hpi-abs(sum(th(indxy).*obj1.l.^4)./sum(obj1.l.^4));
%     obj1.l = sum(obj1.l.^4)/sum(obj1.l.^3);
    %% Method 3 - Fit ellipse to all wing pixels to 
    %% Compute Wing Length & Angle
%     [obj1.l,obj1.phil] = wing_ellipsefit(x1(indxy),y1(indxy),obj1.head,params);
%     obj1.phil = params.hpi-abs(obj1.phil);

    if obj1.l <= 1.15 || obj1.l > 2.5, obj1.l = obj_1.l; obj1.phil = obj_1.phil; end
else
    obj1.l = 0;
    obj1.phil = 0;
end

% ...................(3d)
indxy = find((th > -pi) & (th < -params.hpi*11/9));
obj1.r = sqrt(y1(indxy).^2 + x1(indxy).^2);
ind = find(obj1.r == max(max(obj1.r)));
if ind,
    %% Method 1 (chosen) - Use wing-pixel 
    %% (Pixel with farthest distance to fly body center) to 
    %% Compute Wing Length & Angle
    ind = ind(1);
    obj1.phir = abs(th(indxy(ind)))-params.hpi;
    obj1.r = obj1.r(ind);
    %% Method 2 - weighted by wing-pixel distance to the fourth power
%     obj1.phir = abs(sum(th(indxy).*obj1.r.^4)./sum(obj1.r.^4))-params.hpi;
%     obj1.r = sum(obj1.r.^4)/sum(obj1.r.^3);
    %% Method 3 - Fit ellipse to all wing pixels to 
    %% Compute Wing Length & Angle
%     [obj1.r,obj1.phir] = wing_ellipsefit(x1(indxy),y1(indxy),obj1.head,params);
%     obj1.phir = abs(obj1.phir)-params.hpi;

    if obj1.r <= 1.15 || obj1.r > 2.5, obj1.r = obj_1.r; obj1.phir = obj_1.phir; end
else
    obj1.r = 0;
    obj1.phir = 0;
end

%%
% Finally, the head and tail locations are calculated (4a) from
% the center of the ellipse, the angular direction of the head,
% and the ellipse's major axis (which is equivalent to half the
% ellipse's length). In addition, the fly is fit into a circle,
% resulting into three relevant parameters: the coordinate of the
% center of the circle, as well as its radius (4b).

% ...................(4a)
obj1.xh = real(obj1.xc + cos(obj1.head) * obj1.A);
obj1.yh = real(obj1.yc + sin(obj1.head) * obj1.A);
obj1.xt = real(obj1.xc - cos(obj1.head) * obj1.A);
obj1.yt = real(obj1.yc - sin(obj1.head) * obj1.A);

% obj1.xh = real(obj1.xc + cos(obj1.head) * obj1.length/2);
% obj1.yh = real(obj1.yc + sin(obj1.head) * obj1.length/2);
% obj1.xt = real(obj1.xc - cos(obj1.head) * obj1.length/2);
% obj1.yt = real(obj1.yc - sin(obj1.head) * obj1.length/2);

% ...................(4b)
[Xcc1,Ycc1,R1] = circfit(chamber.x_arr(f01),chamber.y_arr(f01));
obj1.xcc = Xcc1;
obj1.ycc = Ycc1;
obj1.rcc = R1;
end

%% C.2 CircFit 
function [xc,yc,R] = circfit( x , y )
   x=x(:); y=y(:);
   a=[x y ones(size(x))]\(-(x.^2+y.^2));
   xc = -.5*a(1);
   yc = -.5*a(2);
   R  =  sqrt((a(1)^2+a(2)^2)/4-a(3));
end

%% C.3 Fit Ellipse to wing pixel
function [l,phi] = wing_ellipsefit(x,y,phi_head,params)
    [Xce1,Yce1,A01,B01,phi,P] = ellipsefit(x,y);
    l = 2*A01;
    phi = phi-phi_head+params.hpi;
    if phi > pi, phi = phi - params.twopi; end
    if phi < -pi, phi = phi + params.twopi; end
end