%% Detect_copul
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
%% Detect copulation

% COPULATION DETECTION
function [pre_cop,post_cop] = detect_copul(fly_feat,dt,pre_cop,post_cop)
min_duration = 300; %[s] min. copulation duration

if nargin < 3,
    pre_cop = [];
    post_cop = [];
end
a = filtfilt(ones(1,100)/100,1,fly_feat.distc);
window = 250;
std_a = zeros(1,floor(length(a)-window)); mea_a = std_a;
for j=1:floor(length(a)-window),
    mea_a(j) = mean(a(j:j+window));
    std_a(j) = std(a(j:j+window));
end
ind1 = [];
% FIND FRAMES WHERE MEAN/STDEV OF RELATIVE FLY DISTANCE 
% IS SMALL
ind = find((std_a < 0.3) & (mea_a < 2.));
if numel(ind) == 0,
    ind1(1) = -99; ind1(2) = -99;
else
    b = a * 0; b1 = b; 
    b(ind) = 1;    
    % Special case of zeros acceleration of fly 2 - not sure why I had 
    % to add this, catches those frames as probable copulation phase
    b1(fly_feat.obj2.acc == 0) = 1; 
    c = bwlabel(b1);
    maxc = max(c); arr = zeros(1,maxc);
    for i=1:maxc, arr(i) = numel(find(c == i)); end
    inda = find(arr > 10); indd = [];
    for i=1:numel(inda), indd = [indd find(c == inda(i))]; end
    b(indd) = 1;
    
    % TRY ISOLATING THE LARGEST CONNECTED AREAS
    c = bwlabel(conv2(b,ones(1,window*2+1),'same'));
    maxc = max(c); arr = zeros(1,maxc);
    for i=1:maxc,
        arr(i) = numel(find(c == i));
    end
    maxind = find(arr == max(arr));
    if arr(maxind)*dt < min_duration, 
        ind1(1) = -99; ind1(2) = -99; 
    else
        % IF DURATION OF LARGEST CONNECTED AREA IS LONG
        % ENOUGH FIND START/END FRAME = COPULATION START/END
        ind = find(c == maxind);
        ind1 = fly_feat.frame([ind(1)+window ind(end)-window]);
    end
end

pre_cop = [pre_cop ind1(1)]; post_cop = [post_cop ind1(2)];
end