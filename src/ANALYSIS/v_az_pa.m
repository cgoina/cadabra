%% V_az_pa
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
%% Compute parallel and azimuthal fly velocity for all frames

function [v_az,v_pa,phi] = v_az_pa(x,y,chamber,oneobj)
x1 = x(1,:); x2 = x(2,:); y1 = y(1,:); y2 = y(2,:);
v1 = sqrt((x1(2:end)-x1(1:end-1)).^2+(y1(2:end)-y1(1:end-1)).^2); v1 = [v1(1) v1];
th1= atan2(y1(2:end)-y1(1:end-1),x1(2:end)-x1(1:end-1)); th1 = [th1(1) th1];
if oneobj,
    phi = atan2((chamber.height/2-y1),(chamber.width/2-x1)) * 180/pi;
else
    phi = atan2(y2-y1,x2-x1) * 180/pi;
    v2 = sqrt((x2(2:end)-x2(1:end-1)).^2+(y2(2:end)-y2(1:end-1)).^2); v2 = [v2(1) v2];
    th2= atan2(y2(2:end)-y2(1:end-1),x2(2:end)-x2(1:end-1)); th2 = [th2(1) th2];
end
% VELOCITY AZIMUTHAL AND PARALLEL TO THE FLY'S MOVING DIRECTION
% WITH RESPECT TO THE RELATIVE POSITION OF THE OPPONENT
v_az = v1 .* sin(phi*pi/180-th1);
v_pa = v1 .* cos(phi*pi/180-th1);
if ~oneobj,
    v_az = [v_az ; v2 .* sin(phi*pi/180-th2)];
    v_pa = [v_pa ; v2 .* cos(phi*pi/180-th2)];
end
end