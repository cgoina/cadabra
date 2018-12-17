%% Smooth_angle
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
%% Smooth orientations, part of 'analyze_wing_mov'

% SMOOTHES ORIENTATIONS AND COMPUTES ANGULAR VELOCITIES
function [data,phi_velo] = smooth_angle(data,maxit)
% Swap orientations by analyzing the angular velocities
for z=1:maxit,
    phi_velo = filter([1 -1],1,data);
    phi_velo = phi_velo(2:end-1); data = data(2:end-1);
    % Number of positive/negative angles
    if z==1, np = numel(find(data>=0)); nn = numel(find(data<0)); end
    if z==2, tmp=nn; nn=np; np=tmp; end
    if np>nn,
        in = find((phi_velo > -220) & (phi_velo < -165)); if numel(in), data(in) = data(in) + 180; end
    else
        in = find((phi_velo > 165) & (phi_velo < 220)); if numel(in), data(in) = data(in) - 180; end
    end
end
% Smooth orientations and compute angular velocities
data = filter([1 1 1]/3,1,data); phi_velo = filter([1 -1],1,data);
% Remove jumps in angular velocities
in = find(phi_velo < -300); if phi_velo(in), phi_velo(in) = phi_velo(in) + 360; end
in = find(phi_velo > 300); if phi_velo(in), phi_velo(in) = phi_velo(in) - 360; end
end
