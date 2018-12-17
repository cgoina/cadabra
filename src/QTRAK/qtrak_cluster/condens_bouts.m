%% Condens_bouts
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
%% Merge bouts that are apart by less than max_gap

% BOUT MERGING
function behavior = condens_bouts(behavior,params)

data1i= [0 cumsum(behavior.len)];
if length(data1i)>1,    
    data1g = behavior.lab(data1i(2:end-1)+1)-behavior.lab(data1i(2:end-1));

    % FIND BOUTS (SEPARATED BY LARGER GAP)
    ind = [1 find(data1g >= params.mf*params.max_gap)+1 length(data1i)];
    behavior.t = behavior.t(ind(1:end-1));
    behavior.len = behavior.len(ind(1:end-1));
    behavior.phir = behavior.phir(ind(1:end-1));
    behavior.phil = behavior.phil(ind(1:end-1));

    % MERGE EVENTS
    data1i = data1i(ind);
    ind = []; ind0 = []; behavior.number = 0; behavior.len = [];
    for i=1:length(data1i)-1,
        dat = data1i(i)+1:data1i(i+1);
        if length(dat) >= params.min_len,
            ind = [ind dat]; ind0 = [ind0 i];
            behavior.len = [behavior.len length(dat)];
            behavior.number = behavior.number + 1;
        end
    end
    behavior.t = behavior.t(ind0);
    behavior.phir = behavior.phir(ind0);
    behavior.phil = behavior.phil(ind0);
    behavior.ind = behavior.ind(ind);
    behavior.lab = behavior.lab(ind);
    behavior.frm = behavior.frm(ind);
    behavior.tim = behavior.tim(ind);
    behavior.do1 = behavior.do1(ind);
    behavior.do2 = behavior.do2(ind);
    behavior.x1 = behavior.x1(ind);
    behavior.x2 = behavior.x2(ind);
    behavior.y1 = behavior.y1(ind);
    behavior.y2 = behavior.y2(ind);
end