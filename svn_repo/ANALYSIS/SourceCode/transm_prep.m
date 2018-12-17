%% Transm_prep
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
%% Preprocessing actions for computing a transition matrix, part of
%% 'transitions'

% COLLECT ACTION START TIME, END TIME, CODE
function [X,I] = transm_prep(datat, datan, datal, cop, dt, valI, X, I, j, trans_t, iloop)
    ind = [0 cumsum(datan)];
    if datan(j) ~= 0,
        indd = ind(j)+1:ind(j+1);
        if cop.pre(j) > 0,
            % In case of copulation, only take period before begin of
            % copulation into account
            ind = find(datat(indd)+datal(indd)*dt < cop.pre(j)*dt);
            if numel(ind), indd = indd(ind); end
        end
        ind = find((datat(indd) >= (iloop-1)*trans_t) & (datat(indd) < iloop*trans_t));
        if numel(ind) > 1,
            indd = indd(ind);
            X.A = [X.A datat(indd)]; % Start time vector
            X.E = [X.E datat(indd)+datal(indd)*dt]; % End time vector
            I = [I ; ones(1,indd(end)-indd(1)+1)'*valI]; % Action code vector
        end
    end
end
