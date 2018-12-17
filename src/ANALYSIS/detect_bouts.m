%% Detect_bouts
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
%% Detect bouts of actions

% DETECT BOUTS OF ACTIONS
function [bouts,ind1,ind2] = detect_bouts(fly_feat,ind1,ind2,b_len,max_gap,bouts,mf,cut,scat_part)
% Prepare data structure
if isstruct(b_len), % min/max length of bout
    min_len = b_len.min; max_len = b_len.max;
else
    min_len = b_len; max_len = 1e6;
end
if nargin < 8, cut = 0; end % remove front/end part of bout from smoothing filter mask
if nargin < 9, scat_part = 0; end % minimum percentage/100 of detections within a bout
if nargin < 7, mf = 3.; end % merge events with gap < mf*max_gap
params.min_len = min_len; params.max_gap = max_gap; params.mf = mf;

% DATA STRUCTURE
% Used for: Wing, Chase, Court (=Circling), Tussl
% ind = frame index after denoising (see 'analyze_data.m' calls 'den_feat.m')
% lab = absolute frame index of original movie (prior to denoising)
% frm = frame number
% mov = movie (fly pair) number
% obj = fly identity
% t = start time of an action [seconds]
% tim = time stamps over the length of an action [seconds]
% number = number of lunges per genotype
% obj1.x1 = fly 1 x position [mm]
% obj1.y1 = fly 1 y position [mm]
% obj1.x2 = fly 2 x position [mm]
% obj1.y2 = fly 2 y position [mm]
% obj1.do1 = fly 1 head direction [deg]
% obj1.do2 = fly 2 head direction at time when fly 1 performs an action [deg]
% phir = right wing angle [deg]
% phil = left wing angle [deg]
% len = length an action [frames]
% len_unf = sum of frames that were scored positive for an action, 
% before filtering (see below)
% number = number of actions per pair

if ~numel(bouts),
        bouts.obj1.ind = []; bouts.obj1.frm = []; bouts.obj1.lab = []; 
        bouts.obj1.t = []; bouts.obj1.tim = [];
        bouts.obj1.do1 = []; bouts.obj1.do2 = []; 
        bouts.obj2.do1 = []; bouts.obj2.do2 = [];
        bouts.obj1.x1 = []; bouts.obj1.y1 = []; 
        bouts.obj1.x2 = []; bouts.obj1.y2 = []; 
        bouts.obj1.len = []; bouts.obj1.phir = []; bouts.obj1.phil = [];
        bouts.obj1.len_unfilt = []; bouts.obj1.number = [];
        bouts.obj2 = bouts.obj1;
end
bouts0.obj1.t = []; bouts0.obj1.len = [];
bouts0.obj1.phir = []; bouts0.obj1.phil = [];
bouts0.obj1.len_unfilt = []; bouts0.obj1.number = [];
bouts0.obj2 = bouts0.obj1;
if ~isstruct(fly_feat), tim = fly_feat; else tim = fly_feat.time; end

% Length of unfiltered bouts
bouts0.obj1.len_unfilt = numel(find(ind1));
bouts0.obj2.len_unfilt = numel(find(ind2));
ind10 = ind1>0; ind20 = ind2>0;

% CLOSE GAPS
if max_gap, 
    ind1 = conv2(ind1,ones(1,round(max_gap+1)),'same');
    ind2 = conv2(ind2,ones(1,round(max_gap+1)),'same');
end

% DETECT BOUTS AND SETUP DATA STRUCTURE
% FLY 1
lab = bwlabel(ind1); ind1 = lab * 0; indt1 = [];
for i=1:max(lab),
    inda = find(lab == i); 
    % Remove the convolution artifact from bouts (max_gap/2 at start/end of bout)
    inda = inda(1+floor(max_gap/2):end-floor(max_gap/2));
    scat = sum(ind10(inda))/numel(inda);
    if numel(inda) >= round(min_len) && numel(inda) <= round(max_len) && ...
            scat >= scat_part,
        indt1 = [indt1 tim(inda(round(end/2)))];
        ind1(inda) = 1;
        bouts0.obj1.len = [bouts0.obj1.len numel(inda)];
        if isstruct(fly_feat),
            bouts0.obj1.phir = [bouts0.obj1.phir quantile(fly_feat.obj1.phir(inda),.95)];
            bouts0.obj1.phil = [bouts0.obj1.phil quantile(fly_feat.obj1.phil(inda),.95)];        
        end
    end
end
ind10 = ind1;
ind1 = find(ind1>0);

% FLY 2
lab = bwlabel(ind2); ind2 = lab * 0; indt2 = [];
for i=1:max(lab),
    inda = find(lab == i); 
    inda = inda(1+floor(max_gap/2):end-floor(max_gap/2));
    scat = sum(ind20(inda))/numel(inda);
    if numel(inda) >= round(min_len) && numel(inda) <= round(max_len) && ...
            scat >= scat_part,
        if sum(ind10(inda))/numel(inda) < .9, % same beh. already detected for fly 1?
            indt2 = [indt2 tim(inda(round(end/2)))];
            ind2(inda) = 1;
            bouts0.obj2.len = [bouts0.obj2.len numel(inda)];
            if isstruct(fly_feat),
                bouts0.obj2.phir = [bouts0.obj2.phir quantile(fly_feat.obj2.phir(inda),.95)];
                bouts0.obj2.phil = [bouts0.obj2.phil quantile(fly_feat.obj2.phil(inda),.95)];
            end
        end
    end
end
ind2 = find(ind2>0);

bouts0.obj1.t = indt1; bouts0.obj2.t = indt2;
bouts0.obj1.number = numel(indt1); bouts0.obj2.number = numel(indt2);

% GATHER POSITION, TIME, ORIENTATION INFORMATION FROM 
% GIVEN BOUT INDICES
if isstruct(fly_feat),
    % FLY 1
    if numel(ind1),
        bouts0.obj1.ind = []; bouts0.obj1.lab = []; bouts0.obj1.tim = [];
        bouts0.obj1.do1 = []; bouts0.obj1.do2 = [];
        bouts0.obj1.x1 = []; bouts0.obj1.y1 = []; bouts0.obj1.x2 = []; bouts0.obj1.y2 = [];

        bouts0.obj1.ind = ind1; bouts0.obj1.lab = fly_feat.lab(ind1); bouts0.obj1.frm = fly_feat.frame(ind1);
        bouts0.obj1.x1 = fly_feat.obj1.pos_x(ind1); bouts0.obj1.y1 = fly_feat.obj1.pos_y(ind1);
        bouts0.obj1.x2 = fly_feat.obj2.pos_x(ind1); bouts0.obj1.y2 = fly_feat.obj2.pos_y(ind1);
        bouts0.obj1.tim = fly_feat.time(ind1);
        bouts0.obj1.do1 = fly_feat.obj1.headdir(ind1); bouts0.obj1.do2 = fly_feat.obj2.headdir(ind1);

        % CONDENSE BOUT DATA
        bouts0.obj1 = condens_bouts(bouts0.obj1,params);
    end
    % FLY 2
    if numel(ind2),
        bouts0.obj2.ind = []; bouts0.obj2.lab = []; bouts0.obj2.tim = [];
        bouts0.obj2.do1 = []; bouts0.obj2.do2 = [];
        bouts0.obj2.x1 = []; bouts0.obj2.y1 = []; bouts0.obj2.x2 = []; bouts0.obj2.y2 = [];

        bouts0.obj2.ind = ind2; bouts0.obj2.lab = fly_feat.lab(ind2); bouts0.obj2.frm = fly_feat.frame(ind2);
        bouts0.obj2.x1 = fly_feat.obj2.pos_x(ind2); bouts0.obj2.y1 = fly_feat.obj2.pos_y(ind2);
        bouts0.obj2.x2 = fly_feat.obj1.pos_x(ind2); bouts0.obj2.y2 = fly_feat.obj1.pos_y(ind2);
        bouts0.obj2.tim = fly_feat.time(ind2);
        bouts0.obj2.do1 = fly_feat.obj2.headdir(ind2); bouts0.obj2.do2 = fly_feat.obj1.headdir(ind2);
        
        % CONDENSE BOUT DATA
        bouts0.obj2 = condens_bouts(bouts0.obj2,params);
    end
end

% ADD BOUT STRUCTURE TO GIVEN STRUCTURE FROM OTHER MOVIES
% FLY 1
fieldn = fieldnames(bouts0.obj1);
for i=1:numel(fieldn),
    eval(['bouts.obj1.' fieldn{i} ' = [bouts.obj1.' fieldn{i} ' bouts0.obj1.' fieldn{i} '];']);
end
% FLY 2
fieldn = fieldnames(bouts0.obj2);
for i=1:numel(fieldn),
    eval(['bouts.obj2.' fieldn{i} ' = [bouts.obj2.' fieldn{i} ' bouts0.obj2.' fieldn{i} '];']);
end

end