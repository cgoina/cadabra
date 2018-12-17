%% Two2oneflymatr
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
%% Merges transition matrices of two flies to one

% MERGE TRANSITION MATRICES OF TWO FLIES INTO ONE
function [h,lab] = two2oneflymatr(h,lab,flyno)

% FIND TRANSITIONS FLY 1-->1, FLY 2-->2
indx1 = []; indx2 = [];
for i=1:length(lab),
    if strfind(lab{i},'f1'), indx1 = [indx1 i]; end
    if strfind(lab{i},'f2'), indx2 = [indx2 i]; end
end
n = length(indx1); [indxb,ib] = intersect(indx1,indx2);
if ~flyno,
    % ADD TRANSITIONS FLY 1->1 + FLY 2->2
    h0 = h(indx1,indx1) + h(indx2,indx2);
    % Remove double-counted transitions
    h0(ib,ib) = h0(ib,ib) - h(indxb,indxb); 
    h = h0;
elseif flyno == 1221,
    % ADD REACTIONS FLY 1->2 + FLY 2->1
    h = h(indx1,indx2) + h(indx2,indx1);
    h(n,:) = 0; h(:,n) = 0; % exclude from/to 'no behavior'
elseif flyno == 1,
    % ONLY TRANSITIONS FLY 1->1
    h = h(indx1,indx1);
elseif flyno == 2
    % ONLY TRANSITIONS FLY 2->2
    h = h(indx2,indx2);
elseif flyno == 12
    % ADD REACTIONS FLY 1->2
    h = h(indx1,indx2);
    h(n,:) = 0; h(:,n) = 0; % exclude from/to 'no behavior'
elseif flyno == 21
    % ADD REACTIONS FLY 2->1
    h = h(indx2,indx1);
    h(n,:) = 0; h(:,n) = 0; % exclude from/to 'no behavior'
end
lab = lab(indx1);

