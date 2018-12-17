%% Switch_flies
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
%% Fly Switchover
% There are variations in the flies' apparent size in the images
% captured by the camera. These variations can result in a switchover
% in the flies' identities (fly 1 becomes fly 2 and vice versa). Left
% uncorrected, the motion parameter time series would appear erratic.
%
% To remedy this situation, the program computes and compares the
% distance between the current and the previous fly positions (which
% can be obtained from the input argument to _meas_dist_) (1a).
% If switching the flies result in smaller movements, the flies will
% be swapped (1b). Switching is not performed if courtship events are
% expected due to the close proximity of the two flies.
function [obj1,obj2,f01,f02,ind_swap] = switch_flies(obj1,obj2,obj_1,obj_2,f01,f02,img,chamber,params,scale,ind_twoobj)

ind_swap = 0;
if ~params.bool_court,

    % ...................(1a)
    if ind_twoobj,
        r_11 = dist([obj1.xc obj1.yc],[obj_1.xc obj_1.yc]');
        r_12 = dist([obj2.xc obj2.yc],[obj_1.xc obj_1.yc]');
    else
        r_11 = dist([obj1.xc obj1.yc],[obj_1.xc obj_1.yc]');
        r_12 = dist([obj1.xc obj1.yc],[obj_2.xc obj_2.yc]');
    end

    % ...................(1b)
    if (1*r_12 < r_11) && params.bool_nfly,
        tmp = obj1;
        obj1 = obj2;
        obj2 = tmp;
        tmp4 = f01;
        f01 = f02;
        f02 = tmp4;
        ind_swap = 1;
    end
end

%%
% Alternatively, instead of using the previous and current positions,
% the decision to switch the flies can also be based on the parameters
% of the ellipses representing the two flies. The objective is to
% adjust the fly identities such that the interframe changes of these
% parameters are minimized (2).

if (obj1.B/obj2.B < .8) && (obj1.A/obj2.A < .95) && ...
        (obj1.A/obj2.A > .7) && params.bool_court,
    r_12 = dist([obj1.xc obj1.yc],[obj_2.xc obj_2.yc]');
    r_21 = dist([obj2.xc obj2.yc],[obj_1.xc obj_1.yc]');

    % ...................(2)
    if ((r_12 < 1) && (r_21 < 1)) || ((obj_1.B/obj_2.B < .8) && ...
            (obj_1.A/obj_2.A < .95) && (obj_1.A/obj_2.A > .7)),
        tmp = obj1;
        obj1 = obj2;
        obj2 = tmp;
        tmp4 = f01;
        f01 = f02;
        f02 = tmp4;
        if ind_swap,
            ind_swap = 0;
        else
            ind_swap = 1;
        end
    end
end

%%
% Resolve the flies' identities by using the identifying dot 
% if one of the flies has one, and the user explicitly
% specifies this during initialization.

if ind_twoobj && params.bool_dot && params.bool_nfly,
    
    nel = 20;

    % get (x,y) coordinates of along flies' long axis
    j = (1:nel)'*2/nel-1; dr = 1.1*[obj1.A obj2.A] / scale.r;
    x = [obj1.yc obj2.yc] / scale.y;
    y = [obj1.xc obj2.xc] / scale.x;
    dx = sin([obj1.phi obj2.phi]) .* dr;
    dy = cos([obj1.phi obj2.phi]) .* dr;
    tmpx = j*dx; tmpy = j*dy-1;

    % calculate distance between min/max brightness value and 
    % location of minimum brightness
    gr = zeros(2,1); gmi = gr; meagrv = gr;
    for i=1:2,
        indx = round(tmpx(:,i) + x(i));
        indy = round(tmpy(:,i) + y(i));
        in = find(indx <= chamber.nrows & ...
                  indy <= chamber.ncols);
        grv = img(indx(in) + indy(in)*chamber.nrows);
        mgrv = min(grv); meagrv(i) = mean(grv);
        gr(i) = max(grv) - mgrv;
        in = abs(find(grv == mgrv)-nel/2); 
        gmi(i) = min(in); 
    end
    if (gr(1) == 0), gr(1) = 0.0001; end

    % switch flies at following conditions:
    if (abs(gr(2)/gr(1)) > 1.2) && (gmi(2) <= gmi(1)+1) && ...
            (abs(1-meagrv(1)/meagrv(2)) < 0.3) && (gr(2) > 0),
        tmp = obj1;
        obj1 = obj2;
        obj2 = tmp;
        tmp4 = f01;
        f01 = f02;
        f02 = tmp4;
        if (obj_1.head == -99) && (obj_2.head == -99),
            obj_1.xc = obj1.xc;
            obj_1.yc = obj1.yc;
            obj_2.xc = obj2.xc;
            obj_2.yc = obj2.yc;
        end
        if ind_swap,
            ind_swap = 0;
        else
            ind_swap = 1;
        end
    end
end
end

function [d] = dist(x,y)
    d = sqrt( sum( (x-y').^2 ) );
end

function plt_ind(tmpx,tmpy,x,y,img,chamber)
figure(100); clf; hold on;
for i=1:2,
    indx = round(tmpx(:,i) + x(i)); indy = round(tmpy(:,i) + y(i));
    if i==1, col = 'b'; else col = 'r'; end;
    plot(img(indx + indy*chamber.nrows),'Color',col);
end; hold off;
figure(200);
for i=1:1,
    indx = round(tmpx(:,i) + x(i)); indy = round(tmpy(:,i) + y(i));
    a = img * 0; a(indx + indy*chamber.nrows) = 1;
    imagesc(a);
end
figure(201); imagesc(img);
end