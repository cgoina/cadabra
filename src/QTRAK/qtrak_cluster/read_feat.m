%% Read_feat
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
%% Read feature-text files from movies processed the tracking software 
%% (QTrak), detect a selection of actions including lunges (more actions 
%% are detected by 'analyze data'), part of 'analyze_data'

function [fly_feat,intNFrms] = read_feat(FeatureFileName,params)
% intNFrms = number of frames
% params = parameters, defined in GUI ('analysis.m')
%JL0831 add the global variables ind1_count and ind2_count to
%output of the aggress_court.mat.

global ind1_count ind2_count;

if nargin < 2
    params.border = 0;
    params.oneobj = 0;
end

% CHOSEN OPERATION POINT FOR LUNGE DETECTION (DEFAULT = 0.7)
decision_thres = params.tune.lunging.thresh;

% IMPORT FEATURE TABLE (.FEAT) INTO "ALL"
all = importdata(FeatureFileName);
intNFrms = size(all.data,1);
all.data = single(all.data);

% DECLARE ARRAYS
% Binary arrays containing occurring events (1 = event; 0 = no event)
% Lunge events
fly_feat.ind1_count = zeros(1,intNFrms); fly_feat.ind2_count = zeros(1,intNFrms);
ind1_count = zeros(1,intNFrms); ind2_count = zeros(1,intNFrms);
% Wing extensions (wingl = left wing; wingr = right wing; wing = both wings)
fly_feat.obj1.wingl = zeros(1,intNFrms); fly_feat.obj1.wingr = zeros(1,intNFrms);
fly_feat.obj1.wing  = zeros(1,intNFrms); fly_feat.obj2.wingl = zeros(1,intNFrms);
fly_feat.obj2.wingr = zeros(1,intNFrms); fly_feat.obj2.wing  = zeros(1,intNFrms);
fly_feat.obj1.lunge = zeros(1,intNFrms); fly_feat.obj2.lunge = zeros(1,intNFrms);
% Wing threats
fly_feat.obj1.wingthr = zeros(1,intNFrms); fly_feat.obj2.wingthr = zeros(1,intNFrms);
% Chasing
fly_feat.obj1.chase = zeros(1,intNFrms);fly_feat.obj2.chase = zeros(1,intNFrms);
% Jumping
fly_feat.obj1.jump = zeros(1,intNFrms); fly_feat.obj2.jump = zeros(1,intNFrms);
% Circling
fly_feat.obj1.court = zeros(1,intNFrms); fly_feat.obj2.court = zeros(1,intNFrms);

% ASSIGN READ DATA TO STRUCTURE
% objX = fly X
% fly_feat.frame = frames; time = time [seconds];
% headdir = fly head direction [deg]; movedir = fly move direction [deg];
% orient = fly orientation (with 180 degree ambiguity) [deg];
% obj_dirdiff = head direction difference between both flies [deg]
% obj_mvdirdiff = moving direction difference between both flies [deg]
% objXtoYmvdirdiff = difference between moving direction of fly X and the vector from fly X to Y [deg]
% mea = average brightness of fly body
% length = 2*major axis of ellipse [mm], fitted to fly body
% area = area of ellipse [mm2], fitted to fly body
% vel = fly velocity [mm/s]; acc = fly acceleration [mm/s2]
% r = distance between fly and arena center [mm]
% distc = distance between both flies (center of ellipses) [mm]
% disthXtY = distance between head of fly X and abdomen of fly Y [mm]
% disth = head distance between both flies [mm]
% distt = abdomen distance between both flies [mm]
% der_distc,der_disth,der_distt = frame-to-frame change of distc,disth,distt [mm]
% FLength = fly length [mm]; FArea = fly area (measured from body pixels) [mm2]
% pos_change = frame-to-frame fly position change [mm]
% pos_x = x-position [mm]; pos_y = y-position [mm]
% phir = right wing angle [deg]; phil = left wing angle [deg]
% winglr = right wing length [mm]; wingll = left wing length [mm]
% der_disthcXY = frame-to-frame change of distance between head of fly X and center of fly 2 [mm]

fly_feat.frame = all.data(:,1)'; % frames
fly_feat.time = all.data(:,2)'; % time [seconds]
fly_feat.obj1.headdir = all.data(:,3)'; fly_feat.obj2.headdir = all.data(:,4)';
fly_feat.obj1.movedir = all.data(:,5)'; fly_feat.obj2.movedir = all.data(:,6)';
fly_feat.obj1.orient = all.data(:,7)'; fly_feat.obj2.orient = all.data(:,8)'; 
fly_feat.obj_dirdiff = all.data(:,9)'; fly_feat.obj_mvdirdiff = all.data(:,10)';
fly_feat.obj1to2mvdirdiff = all.data(:,11)'; fly_feat.obj2to1mvdirdiff = all.data(:,12)';
fly_feat.obj1.mea = all.data(:,13)'; fly_feat.obj2.mea = all.data(:,14)';
fly_feat.obj1.length = all.data(:,15)'; fly_feat.obj2.length = all.data(:,16)';
fly_feat.obj1.area = all.data(:,17)'; fly_feat.obj2.area = all.data(:,18)';
fly_feat.obj1.vel = all.data(:,19)'; fly_feat.obj2.vel = all.data(:,20)';
fly_feat.obj1.acc = all.data(:,21)'; fly_feat.obj2.acc = all.data(:,22)';
fly_feat.obj1.r = all.data(:,23)'; fly_feat.obj2.r = all.data(:,24)';
fly_feat.distc = all.data(:,25)'; fly_feat.disth1t2 = all.data(:,26)';
fly_feat.disth2t1 = all.data(:,27)'; fly_feat.disth = all.data(:,28)';
fly_feat.distt = all.data(:,29)'; fly_feat.der_distc = all.data(:,30)';
fly_feat.der_disth = all.data(:,31)'; fly_feat.der_distt = all.data(:,32)';
fly_feat.obj1.FArea = all.data(:,33)'; fly_feat.obj2.FArea = all.data(:,34)';
fly_feat.obj1.FLength = all.data(:,35)'; fly_feat.obj2.FLength = all.data(:,36)';
fly_feat.obj1.pos_change = all.data(:,47)'; fly_feat.obj2.pos_change = all.data(:,48)';
fly_feat.obj1.pos_x = all.data(:,49)'; fly_feat.obj1.pos_y = all.data(:,50)';
fly_feat.obj2.pos_x = all.data(:,51)'; fly_feat.obj2.pos_y = all.data(:,52)';
fly_feat.obj1.phir = all.data(:,53)'; fly_feat.obj1.phil = all.data(:,54)';
fly_feat.obj1.winglr = all.data(:,55)'; fly_feat.obj1.wingll = all.data(:,56)';
fly_feat.obj2.phir = all.data(:,57)'; fly_feat.obj2.phil = all.data(:,58)';
fly_feat.obj2.winglr = all.data(:,59)'; fly_feat.obj2.wingll = all.data(:,60)';
fly_feat.der_disthc12 = [0 (fly_feat.disth1t2(2:end)+fly_feat.disth(2:end))/2 - (fly_feat.disth1t2(1:end-1)+fly_feat.disth(1:end-1))/2];
fly_feat.der_disthc21 = [0 (fly_feat.disth2t1(2:end)+fly_feat.disth(2:end))/2 - (fly_feat.disth2t1(1:end-1)+fly_feat.disth(1:end-1))/2];
% REMOVE "ALL" STRUCTURE FROM MEMORY
clear('all');

% FIND AND CORRECT SWAPPED FLY HEAD DIRECTIONS
% 0.5 = 50% weight on fly velocity angle
if params.tune.correct_orient,
    fly_feat = viterbi_vel_ori(fly_feat,0.5);
end

% FIND AND CORRECT SWAPPED FLY IDENTITIES
if params.tune.correct_positions,
    fly_feat = viterbi_pos(fly_feat);
end

% READ ROI AND SCALE INFORMATION
[ROI,scale] = load_ROI(FeatureFileName);
width.x = (ROI.cols(end)-ROI.cols(1))*scale.x;
width.y = (ROI.rows(end)-ROI.rows(1))*scale.y;

% TRANSPOSE RECTANGULAR CHAMBER IN CHASE OF DIFFERENT ORIENTATION
if (width.x > width.y) && (width.y/width.x < .85),
    tmp = fly_feat.obj1.pos_x; fly_feat.obj1.pos_x = fly_feat.obj1.pos_y; fly_feat.obj1.pos_y = tmp;
    
    tmp = 90 - fly_feat.obj1.headdir; ind = find(tmp > 180);
    if numel(ind), tmp(ind) = tmp(ind) - 360; end
    fly_feat.obj1.headdir = tmp;
    
    tmp = 90 - fly_feat.obj1.movedir; ind = find(tmp > 180);
    if numel(ind), tmp(ind) = tmp(ind) - 360; end
    fly_feat.obj1.movedir = tmp;
    
    tmp = 90 - fly_feat.obj1.orient; ind = find(tmp < 0);
    if numel(ind), tmp(ind) = tmp(ind) + 180; end
    fly_feat.obj1.orient = tmp;
    
    tmp = fly_feat.obj2.pos_x; fly_feat.obj2.pos_x = fly_feat.obj2.pos_y; fly_feat.obj2.pos_y = tmp;
    
    tmp = 90 - fly_feat.obj2.headdir; ind = find(tmp > 180);
    if numel(ind), tmp(ind) = tmp(ind) - 360; end
    fly_feat.obj2.headdir = tmp;
    
    tmp = 90 - fly_feat.obj2.movedir; ind = find(tmp > 180);
    if numel(ind), tmp(ind) = tmp(ind) - 360; end
    fly_feat.obj2.movedir = tmp;
    
    tmp = 90 - fly_feat.obj2.orient; ind = find(tmp < 0);
    if numel(ind), tmp(ind) = tmp(ind) + 180; end
    fly_feat.obj2.orient = tmp;
end
if (width.x > width.y), tmp = width.x; width.x = width.y; width.y = tmp; end

% IN CASE OF ONLY ONE FLY RELATE DISTANCES TO ARENA CENTER
if params.oneobj,
    fly_feat.distc = sqrt((fly_feat.obj1.pos_x-width.x/2).^2+(fly_feat.obj1.pos_y.^2-width.y/2));
    fly_feat.der_distc = [0 fly_feat.distc(2:end) - fly_feat.distc(1:end-1)];
end

% COMPUTE AZIMUTHAL AND PARALLEL FLY VELOCITIES
x = [fly_feat.obj1.pos_x ; fly_feat.obj2.pos_x];
y = [fly_feat.obj1.pos_y ; fly_feat.obj2.pos_y];
chamber.height = width.y; chamber.width = width.x;
[v_az,v_pa,phi] = v_az_pa(x,y,chamber,params.oneobj);


% DETECT ACTIONS AT EACH FRAME
for FrameNumber=5:intNFrms-1,
    % If "border" is dectivated, only scan for actions within the arena
    % excluding the "border_width", otherwise take whole ROI into account
    if params.radius > 0,
        r1 = sqrt((fly_feat.obj1.pos_x(FrameNumber)-width.x/2).^2 + ...
                   (fly_feat.obj1.pos_y(FrameNumber)-width.y/2).^2);
        r2 = sqrt((fly_feat.obj2.pos_x(FrameNumber)-width.x/2).^2 + ...
                   (fly_feat.obj2.pos_y(FrameNumber)-width.y/2).^2);
    else
        r1 = 0; r2 = 0;
    end
    if params.border || ...
            (params.radius < 0 && ...
            (fly_feat.obj1.pos_x(FrameNumber) > params.border_width && ...
            fly_feat.obj1.pos_x(FrameNumber) < width.x-params.border_width && ...
            fly_feat.obj1.pos_y(FrameNumber) > params.border_width && ...
            fly_feat.obj1.pos_y(FrameNumber) < width.y-params.border_width) && ...
            (params.oneobj || ...
            fly_feat.obj2.pos_x(FrameNumber) > params.border_width && ...
            fly_feat.obj2.pos_x(FrameNumber) < width.x-params.border_width && ...
            fly_feat.obj2.pos_y(FrameNumber) > params.border_width && ...
            fly_feat.obj2.pos_y(FrameNumber) < width.y-params.border_width)) || ...
            (params.radius > 0 && (r1 < params.radius || r2 < params.radius)),
%             (r1 < params.radius && (params.oneobj || r2 < params.radius)),
        
        % Wing threat ranges
        min_phi_threat = 30; max_phi_threat = 80; min_wingl = 1.1; max_wingl = 1.9; max_v = 5;
        
        % Compute fly head direction difference
        headdirdiff = abs(fly_feat.obj1.headdir(FrameNumber) - fly_feat.obj2.headdir(FrameNumber));
        if headdirdiff > 180, headdirdiff = 360 - headdirdiff; end
        if headdirdiff > 90, headdirdiff = 180 - headdirdiff; end
        
        % Detect frames with a probable wing threat for fly 1
        if ((r1 < params.radius && params.radius > 0) || params.radius < 0) && ...
                (fly_feat.obj1.phir(FrameNumber-1) > min_phi_threat) && (fly_feat.obj1.phir(FrameNumber) > min_phi_threat) && ...
                (fly_feat.obj1.phir(FrameNumber-1) < max_phi_threat) && (fly_feat.obj1.phir(FrameNumber) < max_phi_threat) && ...
                (fly_feat.obj1.phil(FrameNumber-1) > min_phi_threat) && (fly_feat.obj1.phil(FrameNumber) > min_phi_threat) && ...
                (fly_feat.obj1.phil(FrameNumber-1) < max_phi_threat) && (fly_feat.obj1.phil(FrameNumber) < max_phi_threat) && ...
                abs(fly_feat.obj1.phil(FrameNumber)-fly_feat.obj1.phil(FrameNumber-1)) > 0 && ...
                abs(fly_feat.obj1.phir(FrameNumber)-fly_feat.obj1.phir(FrameNumber-1)) > 0 && ...
                (fly_feat.obj1.winglr(FrameNumber-1) > min_wingl) && (fly_feat.obj1.winglr(FrameNumber) > min_wingl) && ...
                (fly_feat.obj1.wingll(FrameNumber-1) > min_wingl) && (fly_feat.obj1.wingll(FrameNumber) > min_wingl) && ...
                (fly_feat.obj1.winglr(FrameNumber-1) < max_wingl) && (fly_feat.obj1.winglr(FrameNumber) < max_wingl) && ...
                (fly_feat.obj1.wingll(FrameNumber-1) < max_wingl) && (fly_feat.obj1.wingll(FrameNumber) < max_wingl) && ...
                (fly_feat.obj1.length(FrameNumber-1) < 1.8) && (fly_feat.obj1.length(FrameNumber) < 1.8) && ...
                (fly_feat.obj1.vel(FrameNumber-1) > 0.01) && (fly_feat.obj1.vel(FrameNumber) > 0.01) && ...
                (fly_feat.obj1.vel(FrameNumber-1) < max_v) && (fly_feat.obj1.vel(FrameNumber) < max_v) && ...
                (fly_feat.distc(FrameNumber) > 2 || (headdirdiff > 30 && fly_feat.distc(FrameNumber) > 1)) && ...
                (fly_feat.distc(FrameNumber) < 30) && ...
                abs(fly_feat.obj1to2mvdirdiff(FrameNumber)) < 100,
            fly_feat.obj1.wingthr(FrameNumber) = 1;
        else
            fly_feat.obj1.wingthr(FrameNumber) = 0;
        end        
        % Detect frames with a probable wing threat for fly 2
        if ((r2 < params.radius && params.radius > 0) || params.radius < 0) && ...
                (fly_feat.obj2.phir(FrameNumber-1) > min_phi_threat) && (fly_feat.obj2.phir(FrameNumber) > min_phi_threat) && ...
                (fly_feat.obj2.phir(FrameNumber-1) < max_phi_threat) && (fly_feat.obj2.phir(FrameNumber) < max_phi_threat) && ...
                (fly_feat.obj2.phil(FrameNumber-1) > min_phi_threat) && (fly_feat.obj2.phil(FrameNumber) > min_phi_threat) && ...
                (fly_feat.obj2.phil(FrameNumber-1) < max_phi_threat) && (fly_feat.obj2.phil(FrameNumber) < max_phi_threat) && ...
                abs(fly_feat.obj2.phil(FrameNumber)-fly_feat.obj2.phil(FrameNumber-1)) > 0 && ...
                abs(fly_feat.obj2.phir(FrameNumber)-fly_feat.obj2.phir(FrameNumber-1)) > 0 && ...
                (fly_feat.obj2.winglr(FrameNumber-1) > min_wingl) && (fly_feat.obj2.winglr(FrameNumber) > min_wingl) && ...
                (fly_feat.obj2.wingll(FrameNumber-1) > min_wingl) && (fly_feat.obj2.wingll(FrameNumber) > min_wingl) && ...
                (fly_feat.obj2.winglr(FrameNumber-1) < max_wingl) && (fly_feat.obj2.winglr(FrameNumber) < max_wingl) && ...
                (fly_feat.obj2.wingll(FrameNumber-1) < max_wingl) && (fly_feat.obj2.wingll(FrameNumber) < max_wingl) && ...
                (fly_feat.obj2.length(FrameNumber-1) < 1.8) && (fly_feat.obj2.length(FrameNumber) < 1.8) && ...
                (fly_feat.obj2.vel(FrameNumber-1) > 0.01) && (fly_feat.obj2.vel(FrameNumber) > 0.01) && ...
                (fly_feat.obj2.vel(FrameNumber-1) < max_v) && (fly_feat.obj2.vel(FrameNumber) < max_v) && ...
                (fly_feat.distc(FrameNumber) > 2 || (headdirdiff > 30 && fly_feat.distc(FrameNumber) > 1)) && ...
                (fly_feat.distc(FrameNumber) < 30) && ...
                abs(fly_feat.obj2to1mvdirdiff(FrameNumber)) < 100, %125
            fly_feat.obj2.wingthr(FrameNumber) = 1;
        else
            fly_feat.obj2.wingthr(FrameNumber) = 0;
        end
        
        if abs(fly_feat.obj1to2mvdirdiff(FrameNumber-1)) > 145, fly_feat.obj1to2mvdirdiff(FrameNumber-1) = 180 - abs(fly_feat.obj1to2mvdirdiff(FrameNumber-1)); end % just in case the last frame was outside borders
        if abs(fly_feat.obj2to1mvdirdiff(FrameNumber-1)) > 145, fly_feat.obj2to1mvdirdiff(FrameNumber-1) = 180 - abs(fly_feat.obj2to1mvdirdiff(FrameNumber-1)); end
        if abs(fly_feat.obj1to2mvdirdiff(FrameNumber)) > 145, fly_feat.obj1to2mvdirdiff(FrameNumber) = 180 - abs(fly_feat.obj1to2mvdirdiff(FrameNumber)); end
        if abs(fly_feat.obj2to1mvdirdiff(FrameNumber)) > 145, fly_feat.obj2to1mvdirdiff(FrameNumber) = 180 - abs(fly_feat.obj2to1mvdirdiff(FrameNumber)); end
        
        % LUNGE DETECTION PREPARATION, STEP 1.1:   
        % Pre-selection of frames with a probable lunge for fly 1        
        % Version 07-29-2008
        if ((r1 < params.radius && params.radius > 0) || params.radius < 0) && ...
                (fly_feat.obj1.length(FrameNumber-1) <= fly_feat.obj1.length(FrameNumber-2)) && ...
                (fly_feat.obj1.area(FrameNumber) > 0.1) && ...
                (fly_feat.obj1.vel(FrameNumber) < 200) && (fly_feat.obj1.vel(FrameNumber) > 0.5) && ...
                (fly_feat.obj2.vel(FrameNumber) < 20) && ...
                (fly_feat.obj1.acc(FrameNumber) < 2000) && (fly_feat.obj1.acc(FrameNumber) > 15) && ...
                (fly_feat.obj1.length(FrameNumber) < 2.5) && (fly_feat.obj1.length(FrameNumber) > 0.8) && ...
                (fly_feat.obj1.length(FrameNumber-2) < 2.5) && (fly_feat.obj1.length(FrameNumber-2) > 0.8) && ...
                fly_feat.obj1.pos_change(FrameNumber) > 0.05 && fly_feat.obj1.pos_change(FrameNumber) < 5 && ...
                abs(fly_feat.obj1to2mvdirdiff(FrameNumber-1)) < 45 && ...
                (fly_feat.der_distc(FrameNumber) > (-2.1)) && (fly_feat.der_distc(FrameNumber) < -0.15) && ...
                fly_feat.distc(FrameNumber) > 0.9 && fly_feat.distc(FrameNumber) < 4,
            ind1_count(FrameNumber) = FrameNumber;
        end
        % Pre-selection of frames with a probable lunge for fly 2
        if ((r2 < params.radius && params.radius > 0) || params.radius < 0) && ...
                (fly_feat.obj2.length(FrameNumber-1) <= fly_feat.obj2.length(FrameNumber-2)) && ...
                (fly_feat.obj2.area(FrameNumber) > 0.1) && ...
                (fly_feat.obj2.vel(FrameNumber) < 200) && (fly_feat.obj2.vel(FrameNumber) > 0.5) && ...
                (fly_feat.obj1.vel(FrameNumber) < 20) && ...
                (fly_feat.obj2.acc(FrameNumber) < 2000) && (fly_feat.obj2.acc(FrameNumber) > 15) && ...
                (fly_feat.obj2.length(FrameNumber) < 2.5) && (fly_feat.obj2.length(FrameNumber) > 0.8) && ...
                (fly_feat.obj2.length(FrameNumber-2) < 2.5) && (fly_feat.obj2.length(FrameNumber-2) > 0.8) && ...
                fly_feat.obj2.pos_change(FrameNumber) > 0.05 && fly_feat.obj2.pos_change(FrameNumber) < 5 && ...
                abs(fly_feat.obj2to1mvdirdiff(FrameNumber-1)) < 45 && ...
                (fly_feat.der_distc(FrameNumber) > (-2.1)) && (fly_feat.der_distc(FrameNumber) < -0.15) && ...
                fly_feat.distc(FrameNumber) > 0.9 && fly_feat.distc(FrameNumber) < 4,
            ind2_count(FrameNumber) = FrameNumber;
        end
        
        % TEST Version 12-2008
%         if ((r1 < params.radius && params.radius > 0) || params.radius < 0) && ...
%                 (fly_feat.obj1.length(FrameNumber-1) <= fly_feat.obj1.length(FrameNumber-2)) && ... % <
%                 (fly_feat.obj1.area(FrameNumber) > 0.1) && ...
%                 (fly_feat.obj1.vel(FrameNumber) < 200) && (fly_feat.obj1.vel(FrameNumber) > 0.5) && ...
%                 (fly_feat.obj2.vel(FrameNumber-4) < 25) && ... % remove!
%                 (fly_feat.obj2.vel(FrameNumber) < 45) && ... % < 20!!!
%                 (fly_feat.obj1.acc(FrameNumber) < 2000) && (fly_feat.obj1.acc(FrameNumber) > 15) && ...
%                 (fly_feat.obj1.length(FrameNumber) < 5) && (fly_feat.obj1.length(FrameNumber) > 0.7) && ... % > 0.8 < 2.5!!!
%                 (fly_feat.obj1.length(FrameNumber-2) < 2.5) && (fly_feat.obj1.length(FrameNumber-2) > 0.8) && ...
%                 fly_feat.obj1.pos_change(FrameNumber) > 0.05 && fly_feat.obj1.pos_change(FrameNumber) < 5 && ...
%                 abs(fly_feat.obj1to2mvdirdiff(FrameNumber-1)) < 45 && ...
%                 (fly_feat.der_distc(FrameNumber) > (-2.1)) && (fly_feat.der_distc(FrameNumber) < -0.15) && ... % -2.1 -- -0.3
%                 fly_feat.distc(FrameNumber) > 0.9 && fly_feat.distc(FrameNumber) < 4,
%             ind1_count(FrameNumber) = FrameNumber;
%         end
%         if ((r2 < params.radius && params.radius > 0) || params.radius < 0) && ...
%                 (fly_feat.obj2.length(FrameNumber-1) <= fly_feat.obj2.length(FrameNumber-2)) && ... % <
%                 (fly_feat.obj2.area(FrameNumber) > 0.1) && ...
%                 (fly_feat.obj2.vel(FrameNumber) < 200) && (fly_feat.obj2.vel(FrameNumber) > 0.5) && ...
%                 (fly_feat.obj1.vel(FrameNumber-4) < 25) && ... % remove!
%                 (fly_feat.obj1.vel(FrameNumber) < 45) && ... % < 20!!!
%                 (fly_feat.obj2.acc(FrameNumber) < 2000) && (fly_feat.obj2.acc(FrameNumber) > 15) && ...
%                 (fly_feat.obj2.length(FrameNumber) < 5) && (fly_feat.obj2.length(FrameNumber) > 0.7) && ... % > 0.8 < 2.5!!!
%                 (fly_feat.obj2.length(FrameNumber-2) < 2.5) && (fly_feat.obj2.length(FrameNumber-2) > 0.8) && ...
%                 fly_feat.obj2.pos_change(FrameNumber) > 0.05 && fly_feat.obj2.pos_change(FrameNumber) < 5 && ...
%                 abs(fly_feat.obj2to1mvdirdiff(FrameNumber-1)) < 45 && ...
%                 (fly_feat.der_distc(FrameNumber) > (-2.1)) && (fly_feat.der_distc(FrameNumber) < -0.15) && ... % -2.1 -- -0.3
%                 fly_feat.distc(FrameNumber) > 0.9 && fly_feat.distc(FrameNumber) < 4,
%             ind2_count(FrameNumber) = FrameNumber;
%         end
    
        % SEE 'WING_EXT.M' FOR CURRENT WING EXTENSION DETECTION
        % FIRST ATTEMPT (2006) IS LISTED HERE:
        min_phi = 55; max_phi = 90; min_wingl = 1.2; max_wingl = 2.; min_flylen = 1.2;
        
        % Fly 1, right wing
        if ((r1 < params.radius && params.radius > 0) || params.radius < 0) && ...
            (fly_feat.obj1.phir(FrameNumber-1) > min_phi) && (fly_feat.obj1.phir(FrameNumber) > min_phi) && ...
                (fly_feat.obj1.phir(FrameNumber-1) < max_phi) && (fly_feat.obj1.phir(FrameNumber) < max_phi) && ...
                (fly_feat.obj1.winglr(FrameNumber-1) > min_wingl) && (fly_feat.obj1.winglr(FrameNumber) > min_wingl) && ...
                (fly_feat.obj1.winglr(FrameNumber-1) < max_wingl) && (fly_feat.obj1.winglr(FrameNumber) < max_wingl) && ...
                (fly_feat.obj1.length(FrameNumber) > min_flylen),
            fly_feat.obj1.wingr(FrameNumber) = 1;
        else
            fly_feat.obj1.wingr(FrameNumber) = 0;
        end
        % Fly 1, left wing
        if ((r1 < params.radius && params.radius > 0) || params.radius < 0) && ...
            (fly_feat.obj1.phil(FrameNumber-1) > min_phi) && (fly_feat.obj1.phil(FrameNumber) > min_phi) && ...
                (fly_feat.obj1.phil(FrameNumber-1) < max_phi) && (fly_feat.obj1.phil(FrameNumber) < max_phi) && ...
                (fly_feat.obj1.wingll(FrameNumber-1) > min_wingl) && (fly_feat.obj1.wingll(FrameNumber) > min_wingl) && ...
                (fly_feat.obj1.wingll(FrameNumber-1) < max_wingl) && (fly_feat.obj1.wingll(FrameNumber) < max_wingl) && ...
                (fly_feat.obj1.length(FrameNumber) > min_flylen),
            fly_feat.obj1.wingl(FrameNumber) = 1;
        else
            fly_feat.obj1.wingl(FrameNumber) = 0;
        end
        % Fly 1, both wings
        if fly_feat.obj1.wingr(FrameNumber) && fly_feat.obj1.wingl(FrameNumber),
            fly_feat.obj1.wing(FrameNumber) = 1;
        end
        % Fly 2, right wing
        if ((r2 < params.radius && params.radius > 0) || params.radius < 0) && ...
            (fly_feat.obj2.phir(FrameNumber-1) > min_phi) && (fly_feat.obj2.phir(FrameNumber) > min_phi) && ...
                (fly_feat.obj2.phir(FrameNumber-1) < max_phi) && (fly_feat.obj2.phir(FrameNumber) < max_phi) && ...
                (fly_feat.obj2.winglr(FrameNumber-1) > min_wingl) && (fly_feat.obj2.winglr(FrameNumber) > min_wingl) && ...
                (fly_feat.obj2.winglr(FrameNumber-1) < max_wingl) && (fly_feat.obj2.winglr(FrameNumber) < max_wingl) && ...
                (fly_feat.obj2.length(FrameNumber) > min_flylen),
            fly_feat.obj2.wingr(FrameNumber) = 1;
        else
            fly_feat.obj2.wingr(FrameNumber) = 0;
        end
        % Fly 2, left wing
        if ((r2 < params.radius && params.radius > 0) || params.radius < 0) && ...
            (fly_feat.obj2.phil(FrameNumber-1) > min_phi) && (fly_feat.obj2.phil(FrameNumber) > min_phi) && ...
                (fly_feat.obj2.phil(FrameNumber-1) < max_phi) && (fly_feat.obj2.phil(FrameNumber) < max_phi) && ...
                (fly_feat.obj2.wingll(FrameNumber-1) > min_wingl) && (fly_feat.obj2.wingll(FrameNumber) > min_wingl) && ...
                (fly_feat.obj2.wingll(FrameNumber-1) < max_wingl) && (fly_feat.obj2.wingll(FrameNumber) < max_wingl) && ...
                (fly_feat.obj2.length(FrameNumber) > min_flylen),
            fly_feat.obj2.wingl(FrameNumber) = 1;
        else
            fly_feat.obj2.wingl(FrameNumber) = 0;
        end
        % Fly 2, both wings
        if fly_feat.obj2.wingr(FrameNumber) && fly_feat.obj2.wingl(FrameNumber),
            fly_feat.obj2.wing(FrameNumber) = 1;
        end
        
        % Jumping (probable frames), fly 1
        if ((r1 < params.radius && params.radius > 0) || params.radius < 0) && ...
            fly_feat.obj1.vel(FrameNumber) > 100 && fly_feat.obj2.vel(FrameNumber) < 20 && ...
                fly_feat.obj1.pos_change(FrameNumber) > 1 && ...
                fly_feat.obj2.pos_change(FrameNumber) > 0 && fly_feat.obj2.pos_change(FrameNumber) < 1 && ...
                fly_feat.obj1.pos_x(FrameNumber) > 5 && fly_feat.obj1.pos_x(FrameNumber) < width.x-5 && ...
                fly_feat.obj1.pos_y(FrameNumber) > 5 && fly_feat.obj1.pos_y(FrameNumber) < width.y-5,
            fly_feat.obj1.jump(FrameNumber) = 1;
        end
        % Jumping (probable frames), fly 2
        if ((r2 < params.radius && params.radius > 0) || params.radius < 0) && ...
            fly_feat.obj2.vel(FrameNumber) > 100 && fly_feat.obj1.vel(FrameNumber) < 20 && ...
                fly_feat.obj2.pos_change(FrameNumber) > 1 && ...
                fly_feat.obj1.pos_change(FrameNumber) > 0 && fly_feat.obj1.pos_change(FrameNumber) < 1 && ...
                fly_feat.obj2.pos_x(FrameNumber) > 5 && fly_feat.obj2.pos_x(FrameNumber) < width.x-5 && ...
                fly_feat.obj2.pos_y(FrameNumber) > 5 && fly_feat.obj2.pos_y(FrameNumber) < width.y-5,
            fly_feat.obj2.jump(FrameNumber) = 1;
        end
        
        % Chasing (probable frames), fly 1
        if ((r1 < params.radius && params.radius > 0) || params.radius < 0) && ...
            abs(fly_feat.obj1to2mvdirdiff(FrameNumber)) < 45 && ...
                fly_feat.disth1t2(FrameNumber) < fly_feat.disth2t1(FrameNumber) && ...
                fly_feat.distc(FrameNumber) > 3 && fly_feat.distc(FrameNumber) < 10 &&  ...
                fly_feat.obj_mvdirdiff(FrameNumber) < 45 && ...
                (fly_feat.obj1.pos_change(FrameNumber) > 0.5) && (fly_feat.obj2.pos_change(FrameNumber) > 0.5) && ...
                abs(fly_feat.der_disthc12(FrameNumber)) < 2 && ...
                (fly_feat.obj1.vel(FrameNumber) > 5) && (fly_feat.obj2.vel(FrameNumber) > 5),
            fly_feat.obj1.chase(FrameNumber) = 1;
        end
        % Chasing (probable frames), fly 2
        if ((r2 < params.radius && params.radius > 0) || params.radius < 0) && ...
            abs(fly_feat.obj2to1mvdirdiff(FrameNumber)) < 45 && ...
                fly_feat.disth2t1(FrameNumber) < fly_feat.disth1t2(FrameNumber) && ...
                fly_feat.distc(FrameNumber) > 3 && fly_feat.distc(FrameNumber) < 10 && ...
                fly_feat.obj_mvdirdiff(FrameNumber) < 45 && ...
                (fly_feat.obj1.pos_change(FrameNumber) > 0.5) && (fly_feat.obj2.pos_change(FrameNumber) > 0.5) && ...
                abs(fly_feat.der_disthc21(FrameNumber)) < 2 && ...
                (fly_feat.obj2.vel(FrameNumber) > 5) && (fly_feat.obj1.vel(FrameNumber) > 5),
            fly_feat.obj2.chase(FrameNumber) = 1;
        end
        
        if params.oneobj,
            % Relate some features of fly 1 to the arena center, in case there is no fly 2
            fly_feat.distc(FrameNumber) = sqrt((fly_feat.obj1.pos_x(FrameNumber)-width.x/2).^2 + (fly_feat.obj1.pos_y(FrameNumber)-width.y/2).^2);
            obj12dir = 180/pi*atan2(fly_feat.obj1.pos_y(FrameNumber)-width.y/2,fly_feat.obj1.pos_x(FrameNumber)-width.x/2);
            fly_feat.obj1to2mvdirdiff(FrameNumber) = (fly_feat.obj1.headdir(FrameNumber) - obj12dir);
            fly_feat.obj2to1mvdirdiff(FrameNumber) = (fly_feat.obj2.headdir(FrameNumber) - (obj12dir-pi));
            if fly_feat.obj1to2mvdirdiff(FrameNumber)>180, fly_feat.obj1to2mvdirdiff(FrameNumber) = fly_feat.obj1to2mvdirdiff(FrameNumber) - 360; end
            if fly_feat.obj1to2mvdirdiff(FrameNumber)<(-180), fly_feat.obj1to2mvdirdiff(FrameNumber) = fly_feat.obj1to2mvdirdiff(FrameNumber) + 360; end
            if fly_feat.obj2to1mvdirdiff(FrameNumber)>180, fly_feat.obj2to1mvdirdiff(FrameNumber) = fly_feat.obj2to1mvdirdiff(FrameNumber) - 360; end
            if fly_feat.obj2to1mvdirdiff(FrameNumber)<(-180), fly_feat.obj2to1mvdirdiff(FrameNumber) = fly_feat.obj2to1mvdirdiff(FrameNumber) + 360; end
            if abs(fly_feat.obj1to2mvdirdiff(FrameNumber)) > 145, fly_feat.obj1to2mvdirdiff(FrameNumber) = 180 - abs(fly_feat.obj1to2mvdirdiff(FrameNumber)); end
            if abs(fly_feat.obj2to1mvdirdiff(FrameNumber)) > 145, fly_feat.obj2to1mvdirdiff(FrameNumber) = 180 - abs(fly_feat.obj2to1mvdirdiff(FrameNumber)); end
            fly_feat.obj2.vel(FrameNumber) = 2;
        end
        % Circling (probable frames), fly 1
        if ((r1 < params.radius && params.radius > 0) || params.radius < 0) && ...
            abs(fly_feat.der_disthc12(FrameNumber)) < 0.5 && abs(fly_feat.der_disthc12(FrameNumber-1)) < 0.5  && ...
                abs(fly_feat.obj1to2mvdirdiff(FrameNumber)) < 20 && ...
                fly_feat.distc(FrameNumber) > 1 && fly_feat.distc(FrameNumber) < 5 &&  ...
                (fly_feat.obj1.vel(FrameNumber) > .25) && (fly_feat.obj2.vel(FrameNumber) < 5) && ...
                ((abs(v_az(1,FrameNumber)) > 0.25) || ...
                (abs(v_az(1,FrameNumber)) > .05 && (fly_feat.obj1.wingr(FrameNumber) || fly_feat.obj1.wingl(FrameNumber)))),
            fly_feat.obj1.court(FrameNumber) = 1;
        end
        % Circling (probable frames), fly 2
        if ~params.oneobj && ((r2 < params.radius && params.radius > 0) || params.radius < 0),
            if abs(fly_feat.der_disthc21(FrameNumber)) < 0.5 && abs(fly_feat.der_disthc21(FrameNumber-1)) < 0.5  && ...
                    abs(fly_feat.obj2to1mvdirdiff(FrameNumber)) < 20 && ...
                    fly_feat.distc(FrameNumber) > 1 && fly_feat.distc(FrameNumber) < 5 &&  ...
                    (fly_feat.obj2.vel(FrameNumber) > .25) && (fly_feat.obj1.vel(FrameNumber) < 5) && ...
                    ((abs(v_az(2,FrameNumber)) > .25) || ...
                    (abs(v_az(2,FrameNumber)) > .05 && (fly_feat.obj2.wingr(FrameNumber) || fly_feat.obj2.wingl(FrameNumber)))),
                fly_feat.obj2.court(FrameNumber) = 1;
            end
        end
    end
end
% Orientation, head direction, and move direction change feature
fly_feat.der_obj1_orient = pi/180*(conv2(fly_feat.obj1.orient,[1 0 -1],'same'));
fly_feat.der_obj2_orient = pi/180*(conv2(fly_feat.obj2.orient,[1 0 -1],'same'));
fly_feat.der_obj1_headdir = pi/180*(conv2(fly_feat.obj1.headdir,[1 0 -1],'same'));
fly_feat.der_obj2_headdir = pi/180*(conv2(fly_feat.obj2.headdir,[1 0 -1],'same'));
fly_feat.der_obj1_movedir = pi/180*(conv2(fly_feat.obj1.movedir,[1 0 -1],'same'));
fly_feat.der_obj2_movedir = pi/180*(conv2(fly_feat.obj2.movedir,[1 0 -1],'same'));

%Save the probableLunge.mat per Eric request

[featureFilePath, ~] = fileparts(FeatureFileName);
probableLungeFileName = fullfile(featureFilePath, 'probablelunge.mat');
save(probableLungeFileName, 'ind1_count', 'ind2_count');

% ==================================================
% LUNGE DETECTION STEPS 1.2 & 2:
% ==================================================
% STEP 1.2:
% Build up a 11-dimensional feature array, finally it 10 dimension will be used
no_feat = 11;
NumOfSamples = numel(ind1_count);
count_ind = zeros(NumOfSamples,1); y_02r = count_ind; ind_lunge0 = count_ind;
x_0r = zeros(NumOfSamples,no_feat);
count = 0;
for i=1:NumOfSamples,
    if ind1_count(i) > 0,
        ind = ind1_count(i);
        x_a1(1:no_feat) = [fly_feat.obj1.vel(ind), fly_feat.obj2.vel(ind), ...
            fly_feat.obj1.acc(ind-1), fly_feat.obj2.acc(ind-1), ...
            fly_feat.obj1.length(ind-2), fly_feat.obj1.length(ind-1), ...
            fly_feat.obj1.length(ind-2)-fly_feat.obj1.length(ind-1), ...
            pi/180*(abs(fly_feat.obj1to2mvdirdiff(ind-1))), ...
            fly_feat.distc(ind), fly_feat.der_distc(ind), fly_feat.obj1.pos_change(ind)];
        if (~numel(find(x_a1 > 1.3e3)) && ~numel(find(x_a1(5:6) > 4)) && ...
                ~numel(find(x_a1(1:2) > 100)) && (pi/180*(x_a1(7)) < 60)),
            count = count + 1;
            count_ind(count) = ind;
            x_0r(count,:) = x_a1;
            y_02r(count) = 1;
            ind_lunge0(count) = fly_feat.frame(ind);
        end
    end
    if ind2_count(i) > 0,
        ind = ind2_count(i);
        x_a1(1:no_feat) = [fly_feat.obj2.vel(ind), fly_feat.obj1.vel(ind), ...
            fly_feat.obj2.acc(ind-1), fly_feat.obj1.acc(ind-1), ...
            fly_feat.obj2.length(ind-2), fly_feat.obj2.length(ind-1), ...
            fly_feat.obj2.length(ind-2)-fly_feat.obj2.length(ind-1), ...
            pi/180*(abs(fly_feat.obj2to1mvdirdiff(ind-1))), ...
            fly_feat.distc(ind), fly_feat.der_distc(ind), fly_feat.obj2.pos_change(ind)];
        if (~numel(find(x_a1 > 1.3e3)) && ~numel(find(x_a1(5:6) > 4)) && ...
                ~numel(find(x_a1(1:2) > 100)) && (pi/180*(x_a1(7)) < 60)),
            count = count + 1;
            count_ind(count) = ind;
            x_0r(count,:) = x_a1;
            y_02r(count) = 2;
            ind_lunge0(count) = fly_feat.frame(ind);
        end
    end
end
count_ind = count_ind(1:count); y_02r = y_02r(1:count); ind_lunge0 = ind_lunge0(1:count);
x_0r = x_0r(1:count,:);

% LUNGE DETECTION, STEP 2:
% EXAMPLE-BASED CLASSIFICATION

% Load the examples data base
bool_kNN = 1;
path2 = '';
DATA_FILE = [path2 'tta1.mat'];

load(DATA_FILE);

% Chose 10 features
x_1 = [x_1(:,1:6) x_1(:,8:11)];
x_0 = [x_0(:,1:6) x_0(:,8:11)];
x_0r = [x_0r(:,1:6) x_0r(:,8:11)];

NN_pos = size(x_1,1);
NN_neg = size(x_0,1);
NN_feat = size(x_1,2);

% Load indices of negative examples
load([path2 'ttr1'], 'rrn');
NN_neg = floor(NN_neg / 7);
x_0 = x_0(rrn(1:NN_neg),:);

% Centering data
for feat=1:NN_feat,
    fmed = median(x_1(:,feat));
    fstd = std(x_1(:,feat));
    x_1(:,feat) = (x_1(:,feat)-fmed) / fstd;
    x_0(:,feat) = (x_0(:,feat)-fmed) / fstd;
    x_0r(:,feat) = (x_0r(:,feat)-fmed) / fstd;
end;

groups = [ones(NN_pos,1) ; zeros(NN_neg,1)];

if (bool_kNN),
    % Sphere Data
    load([path2 'ttC']);
    [U,S,V] = svd(C);
    x_1 = x_1 * U * inv(sqrt(S));
    x_0 = x_0 * U * inv(sqrt(S));
    x_0r = x_0r * U * inv(sqrt(S));
    data_train = [x_1 ; x_0];
    % k-NN classifier preparation
    if nargin<4,
        K = [3 5 9 15 21 27]; %% Number of nearest neighbors considered
        kn = K(4);
    end
    colr = ['r';'g';'b';'k';'c';'m']; N_colr = length(colr);
else
    data_train = [x_1 ; x_0];
    % SVM-Classifier (alternative)
    svmStruct = svmtrain(data_train,groups,'showplot',false,'Kernel_Function', 'polynomial','Polyorder', 5); %% svm
    classes = svmclassify(svmStruct,x_0r); %% classify
end

% K-NN (NEAREST NEIGHBOR) CLASSIFICATION OF LUNGES
if nargin<4, kn = 15; end
fly_feat.obj1.lunge = zeros(1,intNFrms);
fly_feat.obj2.lunge = zeros(1,intNFrms);
N_TEST = length(y_02r); N_TRAIN = length(data_train);
score = zeros(1,N_TEST); classes = score; nolung = score;
% Test all pre-selected (suspicious) frames from part 1
for ii = 1 : N_TEST
    cp = x_0r(ii,:);
    [dd,best_index] = sort(dist2(cp,data_train));
    dd = dd(1:kn);
    best_index = best_index(1:kn);
    cl_label = groups(best_index);
    score(ii) = sum(cl_label); if dd(1) < 1e-4, score(ii) = kn; end
    % Lunge is detected if "score ? decision_thres * kn"
    if (score(ii) >= decision_thres * kn) || dd(1) < 1e-4,
        classes(ii) = 1;
    else
        classes(ii) = 0;
    end
    
    if classes(ii) && (y_02r(ii) == 1),
        % Lunge of fly 1
        fly_feat.obj1.lunge(count_ind(ii)) = 1;
    elseif classes(ii) && (y_02r(ii) == 2),
        % Lunge of fly 2
        fly_feat.obj2.lunge(count_ind(ii)) = 1;
    else
        % No lunge
        nolung(ii) = 1;
    end
end

% DELETE DOUBLE-ASSIGMENTS (BOTH FLIES LUNGE) BY COMPARING SCORES
ind = intersect(find(fly_feat.obj1.lunge),find(fly_feat.obj2.lunge));
if numel(ind),
    for i=1:numel(ind),
        in = find(ind_lunge0 == fly_feat.frame(ind(i)));
        [msc,isc] = max(score(in));
        if y_02r(in(isc)) == 1,
            fly_feat.obj2.lunge(ind(i)) = 0;
        else
            fly_feat.obj1.lunge(ind(i)) = 0;
        end
    end
end
fly_feat.lab = 1:intNFrms;
%%
end