%% Viterbi_pos
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
%% Try to correct fly position switches
% we will set the angle to theta_t = phi_t + s_t * pi
% we want to choose s_t to minimize
% \sum_t cost(s_t|s_{t-1})
% cost(s_t|s_{t-1}) = [(1 - w(||v_t||^2))*d(\theta_t,\theta_{t-1}) +
%                      w(||v_t||^2)*d(\theta_t,angle(v_t))]
% where w(||v_t||^2) = \min{1, c*||v_t||^2}
% we will find the most likely states s_t using the recursion
% cost_t(s_t) = min_{s_{t-1}} { cost_{t-1}(s_{t-1}) + cost(s_t|s_{t-1})

% FIND AND CORRECT SWAPPED FLY IDENTITIES
function fly_feat = viterbi_pos(fly_feat)

verbose = 0;

NFrms = length(fly_feat.frame);

old_fly_feat = fly_feat;

% Allocate space for storing the optimal path
stateprev = zeros(NFrms,2);

% Allocate space for computing costs
tmpcost = zeros(1,2);

% Initialize first frame
costprev = zeros(1,2);

id{1} = ['obj' num2str(1)];
id{2} = ['obj' num2str(2)];

% Compute iteratively
for t=2:NFrms,
    
    % Compute for both possible states
    for scurr=1:2,
        if scurr == 1,
            pos_x(1).curr = fly_feat.(id{1}).pos_x(t);
            pos_y(1).curr = fly_feat.(id{1}).pos_y(t);            
            pos_x(2).curr = fly_feat.(id{2}).pos_x(t);
            pos_y(2).curr = fly_feat.(id{2}).pos_y(t);            
        else
            pos_x(1).curr = fly_feat.(id{2}).pos_x(t);
            pos_y(1).curr = fly_feat.(id{2}).pos_y(t);            
            pos_x(2).curr = fly_feat.(id{1}).pos_x(t);
            pos_y(2).curr = fly_feat.(id{1}).pos_y(t);            
        end
        
        % Try both previous states
        for sprev=1:2,
            if sprev == 1,
                pos_x(1).prev = fly_feat.(id{1}).pos_x(t-1);
                pos_y(1).prev = fly_feat.(id{1}).pos_y(t-1);
                pos_x(2).prev = fly_feat.(id{2}).pos_x(t-1);
                pos_y(2).prev = fly_feat.(id{2}).pos_y(t-1);                
            else
                pos_x(1).prev = fly_feat.(id{2}).pos_x(t-1);
                pos_y(1).prev = fly_feat.(id{2}).pos_y(t-1);
                pos_x(2).prev = fly_feat.(id{1}).pos_x(t-1);
                pos_y(2).prev = fly_feat.(id{1}).pos_y(t-1);                
            end
            
            costcurr = min(sqrt((pos_x(1).curr-pos_x(1).prev).^2+(pos_y(1).curr-pos_y(1).prev).^2), ...
                sqrt((pos_x(2).curr-pos_x(2).prev).^2+(pos_y(2).curr-pos_y(2).prev).^2));
%             costcurr = max([sqrt((pos_x(1).curr-pos_x(1).prev).^2+(pos_y(1).curr-pos_y(1).prev).^2), ...
%                 sqrt((pos_x(2).curr-pos_x(2).prev).^2+(pos_y(2).curr-pos_y(2).prev).^2)]);
%             costcurr = sqrt((pos_x(1).curr-pos_x(1).prev).^2+(pos_y(1).curr-pos_y(1).prev).^2);
%             costcurr = sqrt((pos_x(2).curr-pos_x(2).prev).^2+(pos_y(2).curr-pos_y(2).prev).^2);
            
            tmpcost(sprev) = costprev(sprev) + costcurr;
        end
        
        % Choose argmin
        [m,sprev] = min(tmpcost);
        
        % Set pointer for path
        stateprev(t-1,scurr) = sprev;
        
        % Set cost
        costprev(scurr) = tmpcost(sprev);
    end
    
    % Choose the best last state
    [m,scurr] = min(costprev);
    if scurr == 2,
        fly_feat = swap_pos(fly_feat,id,NFrms-1);
    end
end

% Choose the best states
iswap = zeros(NFrms,1);
for t=NFrms-1:-1:1,
    scurr = stateprev(t,scurr);
    if scurr == 2,
        iswap(t) = 1;
    end
end

% SWAP FLIES
iswap = logical(iswap);
if sum(iswap)>0, 
    fly_feat = swap_pos(fly_feat,id,iswap); 
end

% SWAP DATA STRUCTURE FOR FRAMES WITH SWAPPED FLIES
fieldn = fieldnames(fly_feat.obj1);
for i=1:numel(fieldn),
    if ~numel(strfind(cell2mat(fieldn(i)),'pos')),
        tmp = fly_feat.(id{1}).(fieldn{i})(iswap);
        fly_feat.(id{1}).(fieldn{i})(iswap) = fly_feat.(id{2}).(fieldn{i})(iswap);
        fly_feat.(id{2}).(fieldn{i})(iswap) = tmp;
    end
end
tmp = fly_feat.disth1t2(iswap); 
fly_feat.disth1t2(iswap) = fly_feat.disth2t1(iswap);
fly_feat.disth2t1(iswap) = tmp;

% IDENTIFY FLY 1 BY COMPARING THE SUM OF THE SQUARED DIFFERENCES OF
% NEW FLY POSITIONS AND THE ORIGINAL POSITIONS
diff1 = sum((fly_feat.(id{1}).pos_x - old_fly_feat.(id{1}).pos_x).^2 + ...
             (fly_feat.(id{1}).pos_y - old_fly_feat.(id{1}).pos_y).^2);
diff2 = sum((fly_feat.(id{2}).pos_x - old_fly_feat.(id{1}).pos_x).^2 + ...
             (fly_feat.(id{2}).pos_y - old_fly_feat.(id{1}).pos_y).^2);
if diff2 < diff1,
    tmp = fly_feat.obj1;
    fly_feat.obj1 = fly_feat.obj2;
    fly_feat.obj2 = tmp;
end

% RE-COMPUTE AFFECTED FEATURES
fly_feat = compute_features(fly_feat,id);

% VISUAL CHECK OF CHANGES IN FLIES X-POSITIONS
if verbose,
    strfrm = 1; endfrm = 30000;
    figure(100); clf;
    ax(1) = subplot(2,1,1); hold on;
    plot(old_fly_feat.obj1.pos_x(strfrm:endfrm),'b');
    plot(old_fly_feat.obj2.pos_x(strfrm:endfrm),'r'); hold off;
    ax(2) = subplot(2,1,2); hold on; 
    plot(fly_feat.obj1.pos_x(strfrm:endfrm),'Color',[0 .2 1]);
    plot(fly_feat.obj2.pos_x(strfrm:endfrm),'Color',[1 .2 0]); hold off;
    linkaxes(ax,'xy');
    drawnow;
end

end


function fly_feat = swap_pos(fly_feat,id,frame)

tmp_x = fly_feat.(id{1}).pos_x(frame); tmp_y = fly_feat.(id{1}).pos_y(frame);
fly_feat.(id{1}).pos_x(frame) = fly_feat.(id{2}).pos_x(frame);
fly_feat.(id{1}).pos_y(frame) = fly_feat.(id{2}).pos_y(frame);
fly_feat.(id{2}).pos_x(frame) = tmp_x;
fly_feat.(id{2}).pos_y(frame) = tmp_y;

end

function fly_feat = compute_features(fly_feat,id)

dt = median(fly_feat.time(2:end) - fly_feat.time(1:end-1));
for i=1:2,
    mx = [0 fly_feat.(id{i}).pos_x(2:end)-fly_feat.(id{i}).pos_x(1:end-1)];
    my = [0 fly_feat.(id{i}).pos_y(2:end)-fly_feat.(id{i}).pos_y(1:end-1)];
    vx = conv2(mx,[.25 .5 .25],'same')/dt;
    vy = conv2(my,[.25 .5 .25],'same')/dt;
    fly_feat.(id{i}).vel = sqrt(vx.^2 + vy.^2);
    fly_feat.(id{i}).acc = sqrt((conv2(vx,[-1 0 1],'same')/2/dt).^2 + ...
                                (conv2(vy,[-1 0 1],'same')/2/dt).^2);
    fly_feat.(id{i}).pos_change = sqrt(mx.^2+my.^2);
    fly_feat.(id{i}).movedir = atan2(my,mx)*180/pi;
end
fly_feat.obj_mvdirdiff = anglemod(abs(fly_feat.(id{1}).movedir - fly_feat.(id{2}).movedir));
obj12dir = atan2(fly_feat.(id{2}).pos_y - fly_feat.(id{1}).pos_y, ...
                 fly_feat.(id{2}).pos_x - fly_feat.(id{1}).pos_x);
fly_feat.obj1to2mvdirdiff = anglemod(180/pi*(fly_feat.(id{1}).headdir*pi/180  - obj12dir));
fly_feat.obj2to1mvdirdiff = anglemod(180/pi*(fly_feat.(id{2}).headdir*pi/180 - (obj12dir-pi)));
fly_feat.der_disthc12 = [0 (fly_feat.disth1t2(2:end)+fly_feat.disth(2:end))/2 - ...
                           (fly_feat.disth1t2(1:end-1)+fly_feat.disth(1:end-1))/2];
fly_feat.der_disthc21 = [0 (fly_feat.disth2t1(2:end)+fly_feat.disth(2:end))/2 - ...
                           (fly_feat.disth2t1(1:end-1)+fly_feat.disth(1:end-1))/2];

end

function angle = anglemod(angle)

if angle > 180, angle = angle - 360; end
if angle <= (-180), angle = angle + 360; end

end
