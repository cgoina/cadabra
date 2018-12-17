%% Plot_events
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
%% Part of 'plot_clips', receives information from 'read_clips' and passes
%% the processed information to 'mmread_clips'

% FURTHER PROCESSING OF DATA FROM "PLOT_CLIPS" BEFORE CALLING
% ACTION CLIP OUTPUT
function plot_events(fly_feat,feature,nevents,clnumber,ncl,cl,ifile,nfiles,min_len,feat_filen,path,add_path,add,params)

if ispc, slash = '\'; else slash = '/'; end
frames = zeros(1,nevents); nframe = frames; frind = cell(nevents,1);
eval(['data = fly_feat.' feature ';']); nfr = length(data); cnt = 1;
for i=1:nevents,
    % EXTRACT START FRAME NUMBER AND LENGTH OF EACH BOUT
    in = find(data == i);
    if numel(in) >= min_len,
        if exist('add'),  %#ok<EXIST>
            in = [in(1)-add:in(1)-1 in in(end)+1:in(end)+add];  %#ok<AGROW>
            in = in((in > 0) & (in <= nfr));
        end
        if strcmp(feature,'copulation'),
            % Add frames for copulation to watch attemps,
            % mounting, unmounting
            frms = 600;
            ina = in(1)-frms/2+1:in(1)+frms/2;
            ina = ina((ina > 0) & (ina <= nfr));
            ine = in(end)-frms/2+1:in(end)+frms/2;
            ine = ine((ine > 0) & (ine <= nfr));
            frames(cnt) = fly_feat.frame(ina(1)); frind{cnt} = ina; nframe(cnt) = length(ina);
            frames(cnt+1) = fly_feat.frame(ine(1)); frind{cnt+1} = ine; nframe(cnt+1) = length(ine);
            cnt = cnt + 2;
        else
            frames(i) = fly_feat.frame(in(1)); frind{i} = in; nframe(i) = length(in);
        end
    end
end
in = find(nframe >= min_len); frames = frames(in); frind = frind(in); nframe = nframe(in);
% Load ROI and scale information
[ROI,scale,fname] = load_ROI(feat_filen);
labels = frames * 0 + clnumber;

% COMPUTE POSITION AND TIME INFORMATION TO LOCALIZE THE AREA
% WHERE THE ACTIONS TOOK PLACE
if numel(frames),
    width.x = (ROI.cols(end)-ROI.cols(1))*scale.x;
    width.y = (ROI.rows(end)-ROI.rows(1))*scale.y;
    if (width.x > width.y) && (width.y/width.x < .85), % transpose chamber in chase of different orientation
        tmp = fly_feat.obj1.pos_x; fly_feat.obj1.pos_x = fly_feat.obj1.pos_y; fly_feat.obj1.pos_y = tmp;
        tmp = fly_feat.obj2.pos_x; fly_feat.obj2.pos_x = fly_feat.obj2.pos_y; fly_feat.obj2.pos_y = tmp;
    end
    pos.xc1 = cell(length(frind),1); pos.yc1 = pos.xc1; pos.xc2 = pos.xc1; pos.yc2 = pos.xc1;
    timestamp = pos.xc1;
    for i=1:length(frind),
        pos.xc1{i} = fly_feat.obj1.pos_x(frind{i}) / scale.x + ROI.cols(1);
        pos.yc1{i} = fly_feat.obj1.pos_y(frind{i}) / scale.y + ROI.rows(1);
        pos.xc2{i} = fly_feat.obj2.pos_x(frind{i}) / scale.x + ROI.cols(1);
        pos.yc2{i} = fly_feat.obj2.pos_y(frind{i}) / scale.y + ROI.rows(1);
        timestamp{i} = fly_feat.time(frind{i});
    end
    
    if exist('add_path'), path = [path add_path slash]; [s,mess,messid] = mkdir(path); end

    if isempty(cl),
        [s,mess,messid] = mkdir([path 'cluster_' num2str(clnumber,'%03d')]);
    else
        [s,mess,messid] = mkdir([path cl]);
    end

    % CALL ACTION CLIP OUTPUT ROUTINE
    mmread_clips(fname, frames, labels, pos, timestamp, nframe, path, clnumber, ncl, cl, ifile, nfiles, params);
end

end