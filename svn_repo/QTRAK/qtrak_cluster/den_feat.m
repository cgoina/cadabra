%% Den_feat
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
%% Remove outliers - frames with measured high velocities (jumping, flying)

% DENOISE TRACKING DATA BY REMOVING OUTLIERS, SUCH AS
% JUMPING AND FLYING
function [fly_feat, NFrms] = den_feat(fly_feat,quant,MARGIN,dt,params,width)

if nargin < 6,
    params.border = 1;
end
bord = params.border_width;

Z1 = fly_feat.obj1.pos_change; Z2 = fly_feat.obj2.pos_change;
if sum(Z2) > 1 && ~params.oneobj,
    % MIN/MAX ALLOWED DISPLACEMENTS
    fminv = min(quantile([fly_feat.obj1.vel fly_feat.obj2.vel],quant)); % >0 mm
    % fmaxv = max(quantile([fly_feat.obj1.vel fly_feat.obj2.vel],1-quant));
    fmaxv = 100; % mm
    
    % MIN/MAX ALLOWED VELOCITIES
    fmin = min(quantile([Z1 Z2],quant)); % >0 mm/s
    % fmax = max(quantile([Z1 Z2],1-quant));
    fmax = fmaxv * dt; % < ~30 mm/s

    fdel = fmax-fmin; fdelv = fmaxv-fminv;
    
    if params.border,
        I_GOOD = find((Z1 >= fmin-MARGIN*fdel & Z1 < fmax+MARGIN*fdel) & ...
            (Z2 >= fmin-MARGIN*fdel & Z2 < fmax+MARGIN*fdel) & ...
            (fly_feat.obj1.vel >= fminv-MARGIN*fdelv & fly_feat.obj1.vel < fmaxv+MARGIN*fdelv) & ...
            (fly_feat.obj2.vel >= fminv-MARGIN*fdelv & fly_feat.obj2.vel < fmaxv+MARGIN*fdelv));
    else
        if params.radius > 0,
            r1 = sqrt(fly_feat.obj1.pos_x.^2 + fly_feat.obj1.pos_y.^2);
            r2 = sqrt(fly_feat.obj2.pos_x.^2 + fly_feat.obj2.pos_y.^2);
            I_GOOD = find((Z1 >= fmin-MARGIN*fdel & Z1 < fmax+MARGIN*fdel) & ...
                (Z2 >= fmin-MARGIN*fdel & Z2 < fmax+MARGIN*fdel) & ...
                (fly_feat.obj1.vel >= fminv-MARGIN*fdelv & fly_feat.obj1.vel < fmaxv+MARGIN*fdelv) & ...
                (fly_feat.obj2.vel >= fminv-MARGIN*fdelv & fly_feat.obj2.vel < fmaxv+MARGIN*fdelv) & ...
                (r1 < params.radius) & (r2 < params.radius));
        else
            I_GOOD = find((Z1 >= fmin-MARGIN*fdel & Z1 < fmax+MARGIN*fdel) & ...
                (Z2 >= fmin-MARGIN*fdel & Z2 < fmax+MARGIN*fdel) & ...
                (fly_feat.obj1.vel >= fminv-MARGIN*fdelv & fly_feat.obj1.vel < fmaxv+MARGIN*fdelv) & ...
                (fly_feat.obj2.vel >= fminv-MARGIN*fdelv & fly_feat.obj2.vel < fmaxv+MARGIN*fdelv) & ...
                fly_feat.obj1.pos_x > -width.x/2+bord & fly_feat.obj1.pos_x < width.x/2-bord & ...
                fly_feat.obj1.pos_y > -width.y/2+bord & fly_feat.obj1.pos_y < width.y/2-bord & ...
                fly_feat.obj2.pos_x > -width.x/2+bord & fly_feat.obj2.pos_x < width.x/2-bord & ...
                fly_feat.obj2.pos_y > -width.y/2+bord & fly_feat.obj2.pos_y < width.y/2-bord);
        end
    end
else
    fmin = min(quantile(Z1,quant));
    % fmax = max(quantile(Z1,1-quant));
    fmax = 10;

    fminv = min(quantile(fly_feat.obj1.vel,quant));
    % fmaxv = max(quantile(fly_feat.obj1.vel,1-quant));
    fmaxv = 100;
    fdel = fmax-fmin; fdelv = fmaxv-fminv;

    if params.border,
        I_GOOD = find((Z1 >= fmin-MARGIN*fdel & Z1 < fmax+MARGIN*fdel) & ...
            (fly_feat.obj1.vel >= fminv-MARGIN*fdelv & fly_feat.obj1.vel < fmaxv+MARGIN*fdelv));
    else
        if params.radius > 0,
            r1 = sqrt((fly_feat.obj1.pos_x-width.x/2).^2 + ...
                (fly_feat.obj1.pos_y-width.y/2).^2);
            I_GOOD = find((Z1 >= fmin-MARGIN*fdel & Z1 < fmax+MARGIN*fdel) & ...
                (fly_feat.obj1.vel >= fminv-MARGIN*fdelv & fly_feat.obj1.vel < fmaxv+MARGIN*fdelv) & ...
                (r1 < params.radius));
        else
            I_GOOD = find((Z1 >= fmin-MARGIN*fdel & Z1 < fmax+MARGIN*fdel) & ...
                (fly_feat.obj1.vel >= fminv-MARGIN*fdelv & fly_feat.obj1.vel < fmaxv+MARGIN*fdelv) & ...
                fly_feat.obj1.pos_x > -width.x/2+bord & fly_feat.obj1.pos_x < width.x/2-bord & ...
                fly_feat.obj1.pos_y > -width.y/2+bord & fly_feat.obj1.pos_y < width.y/2-bord);
        end
    end
end

% REMOVE OUTLIERS FROM DATA STRUCTURE
fieldn = fieldnames(fly_feat);
for i=1:numel(fieldn),
    if ~strcmp(fieldn(i),'obj1') && ~strcmp(fieldn(i),'obj2'),
        eval(['fly_feat.' fieldn{i} '=fly_feat.' fieldn{i} '(I_GOOD);']);
    end
end
fieldn = fieldnames(fly_feat.obj1);
for i=1:numel(fieldn),
    eval(['fly_feat.obj1.' fieldn{i} '=fly_feat.obj1.' fieldn{i} '(I_GOOD);']);
    eval(['fly_feat.obj2.' fieldn{i} '=fly_feat.obj2.' fieldn{i} '(I_GOOD);']);
end
NFrms = length(fly_feat.frame);

end