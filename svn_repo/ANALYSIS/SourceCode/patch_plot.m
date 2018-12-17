%% Patch_plot
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
%% Plot single bars as part of redisigned barplot routine

% PLOT BARS (my replacement for Matlabs 'barplot' function)
function patch_plot(x,y,params)
if ~params.oneobj && (size(y,2) == 3),
    d1 = -params.bardispl; d2 = -d1; w = params.twobarswidth/2;
    for i=1:size(y,1),
        patch([x(i)+d1-w x(i)+d1+w x(i)+d1+w x(i)+d1-w x(i)+d1-w],...
            [0 0 y(i,1) y(i,1) 0],params.flycol.win); hold on;
        patch([x(i)-w x(i)+w x(i)+w x(i)-w x(i)-w],...
            [0 0 y(i,2) y(i,2) 0],params.flycol.los);
        patch([x(i)+d2-w x(i)+d2+w x(i)+d2+w x(i)+d2-w x(i)+d2-w],...
            [0 0 y(i,3) y(i,3) 0],[0.5 .8 .76]);
    end
elseif ~params.oneobj && (size(y,2) == 2),
    d1 = -params.bardispl; d2 = -d1; w = params.twobarswidth;
    for i=1:size(y,1),
        patch([x(i)+d1-w x(i)+d1+w x(i)+d1+w x(i)+d1-w x(i)+d1-w],...
            [0 0 y(i,1) y(i,1) 0],params.flycol.win); hold on;
        patch([x(i)+d2-w x(i)+d2+w x(i)+d2+w x(i)+d2-w x(i)+d2-w],...
            [0 0 y(i,2) y(i,2) 0],params.flycol.los);
    end
else
    d1 = 0; w = params.barwidth;
    for i=1:size(y,1),
        patch([x(i)+d1-w x(i)+d1+w x(i)+d1+w x(i)+d1-w x(i)+d1-w],...
            [0 0 y(i,1) y(i,1) 0],params.flycol.win); hold on;
    end
end