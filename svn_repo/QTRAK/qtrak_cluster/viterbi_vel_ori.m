%% Viterbi_vel_ori
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
%% Try to correct head direction switches
% we will set the angle to theta_t = phi_t + s_t * pi
% we want to choose s_t to minimize
% \sum_t cost(s_t|s_{t-1})
% cost(s_t|s_{t-1}) = [(1 - w(||v_t||^2))*d(\theta_t,\theta_{t-1}) +
%                      w(||v_t||^2)*d(\theta_t,angle(v_t))]
% where w(||v_t||^2) = \min{1, c*||v_t||^2}
% we will find the most likely states s_t using the recursion
% cost_t(s_t) = min_{s_{t-1}} { cost_{t-1}(s_{t-1}) + cost(s_t|s_{t-1})

% HEAD DIRECTION CORRECTION
function fly_feat = viterbi_vel_ori(fly_feat,vel_weight)
% vel_weight .. weight (0..1) of velocity (fly move) direction

verbose = 0;
if verbose, old_fly_feat = fly_feat; end

NFrms = length(fly_feat.frame);

iswap = zeros(NFrms,2);

for ident=1:2,
    
    % Fly ID
    id = ['obj' num2str(ident)];    
    
    % ALLOCATE SPACE FOR STORING THE OPTIMAL PATH
    stateprev = zeros(NFrms,2);
    
    % ALLOCATE SPACE FOR COMPUTING COSTS
    tmpcost = zeros(1,2);
    
    % INITIALIZE FIRST FRAME
    costprev = zeros(1,2);
    
    % COMPUTE ITERATIVELY
    for t=2:NFrms,
        
        % COMPUTE FOR BOTH POSSIBLE STATES
        for scurr=1:2,
            theta.curr = fly_feat.(id).headdir(t) * pi/180 + pi*(scurr-1);
            x.curr = fly_feat.(id).pos_x(t);
            y.curr = fly_feat.(id).pos_y(t);
            
            % TRY BOTH PREVIOUS STATES
            for sprev=1:2,
                theta.prev = fly_feat.(id).headdir(t-1) * pi/180 + pi*(sprev-1);
                x.prev = fly_feat.(id).pos_x(t-1);
                y.prev = fly_feat.(id).pos_y(t-1);
                
                vx = x.curr - x.prev;
                vy = y.curr - y.prev;
                
                w = min(1,vel_weight*sqrt(vx.^2 + vy.^2));
                velocityangle = atan2(vy,vx);
                
                costcurr = (1.-w) * angledist(theta.prev,theta.curr) + ...
                            w * angledist(theta.curr,velocityangle);
                
                tmpcost(sprev) = costprev(sprev) + costcurr;
            end
            
            % CHOOSE ARGMIN
            [m,sprev] = min(tmpcost);
            
            % SET POINTER FOR PATH
            stateprev(t-1,scurr) = sprev;
            
            % SET COST
            costprev(scurr) = tmpcost(sprev);
        end
        
        % CHOOSE THE BEST LAST STATE
        [m,scurr] = min(costprev);
        if scurr == 2,
            fly_feat.(id).headdir(NFrms-1) = fly_feat.(id).headdir(NFrms-1) + 180;
            fly_feat.(id).headdir(NFrms-1) = anglemod(fly_feat.(id).headdir(NFrms-1));
        end
    end
    
    % CHOOSE THE BEST STATES
    for t=NFrms-1:-1:1,
        scurr = stateprev(t,scurr);
        if scurr == 2,
            iswap(t,ident) = 1;            
        end
    end
    
    ind = logical(iswap(:,ident));
    
    % SWITCH HEAD DIRECTION
    fly_feat.(id).headdir(ind) = fly_feat.(id).headdir(ind) + 180;
    fly_feat.(id).headdir(ind) = anglemod(fly_feat.(id).headdir(ind));
    % ALSO SWITCH LEFT/RIGHT WING
    tmp_phir = fly_feat.(id).phir(ind); tmp_winglr = fly_feat.(id).winglr(ind);
    fly_feat.(id).phir(ind) = fly_feat.(id).phil(ind);
    fly_feat.(id).winglr(ind) = fly_feat.(id).wingll(ind);
    fly_feat.(id).phil(ind) = tmp_phir;
    fly_feat.(id).wingll(ind) = tmp_winglr;

    if verbose,
        figure(100+ident);
        subplot(2,1,1); plot(old_fly_feat.(id).headdir(1:5000),'r');
        subplot(2,1,2); plot(fly_feat.(id).headdir(1:5000),'g');
        drawnow;
    end
end

% REARRANGE DISTANCES THAT DEPEND ON HEAD-TAIL LOCATION
% Case 1: fly 1 swapped
ind = logical(iswap(:,1) == 1 & iswap(:,2) == 0);
tmp = fly_feat.disth(ind); 
fly_feat.disth(ind) = fly_feat.disth2t1(ind);
fly_feat.disth2t1(ind) = tmp;
tmp = fly_feat.distt(ind); 
fly_feat.distt(ind) = fly_feat.disth1t2(ind);
fly_feat.disth1t2(ind) = tmp;

% Case 2: fly 2 swapped
ind = logical(iswap(:,1) == 0 & iswap(:,2) == 1);
tmp = fly_feat.disth(ind); 
fly_feat.disth(ind) = fly_feat.disth1t2(ind);
fly_feat.disth1t2(ind) = tmp;
tmp = fly_feat.distt(ind); 
fly_feat.distt(ind) = fly_feat.disth2t1(ind);
fly_feat.disth2t1(ind) = tmp;

% Case 3: both flies swapped
ind = logical(iswap(:,1) == 1 & iswap(:,2) == 1);
tmp = fly_feat.disth(ind); 
fly_feat.disth(ind) = fly_feat.distt(ind);
fly_feat.distt(ind) = tmp;

end


function angle = angledist(angle1,angle2)

% angle = dist(angle1,angle2);
angle = abs(angle1-angle2);
angle = abs(anglemod(angle*180/pi))*pi/180;

end

function angle = anglemod(angle)

ind = logical(angle > 180); angle(ind) = angle(ind) - 360;
ind = logical(angle <= (-180)); angle(ind) = angle(ind) + 360;

end
