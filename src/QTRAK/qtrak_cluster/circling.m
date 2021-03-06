%% Circling
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
%% Analyze CIRCLING action

function court = circling(fly_feat,b_len,maxg,court)

if ~exist('court'),
    court = [];
end

NFrms = length(fly_feat.frame);
dt = fly_feat(1).time(2:NFrms(1)) - fly_feat(1).time(1:NFrms(1)-1);
dt = median(dt);

if isstruct(b_len), 
    b_len.min = round(b_len.min/dt); b_len.max = round(b_len.max/dt); 
else
    b_len = round(b_len/dt);
end
maxg = round(maxg/dt);
ind1 = zeros(1,NFrms); ind1(fly_feat.obj1.court > 0) = 1;
ind2 = zeros(1,NFrms); ind2(fly_feat.obj2.court > 0) = 1;
% DETECT BOUTS OF ACTION AND SETUP DATA STRUCTURE
[court,i1,i2] = detect_bouts(fly_feat,ind1,ind2,b_len,maxg,court,3);

end