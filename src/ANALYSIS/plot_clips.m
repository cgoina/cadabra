%% Plot_clips
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
%% MAIN ROUTINE FOR PROCESSING MOVIE CLIPS OF DETECTED ACTIONS, 
%% display them as rows of animated GIFs with time stamp on a webpage,
%% start the default web browser showing links to the action movie clips.
%% Each clip is shown with half the original resolution for memory and
%% CPU-time saving reasons.

function plot_clips(path,sorted_names,gen_ident,genotypes,params)
global indxg;

if params.plots.movieclips,
    % Always re-initialize clip output (=1) or continue a lastest known position (=0)
    % (=0): Routine would continuously save the output position (movie, clip
    % index) to continue after, i.e. a program reset, break, etc.
    params.init = 1;
    
    % Which action?
    fieldn = fieldnames(params.plots.stat);
    for i=1:numel(fieldn),
        if numel(params.plots.stat.(fieldn{i})) == 1,
            if params.plots.stat.(fieldn{i}), action = fieldn{i}; end
        end
    end
    
    add_frames = 0; % add frames before/after each detected bout;
    add_path = 'clips'; % sub-folder for the video-clips to be stored
    min_len = 20; % minmal length of events in frames

    % Minimum clip length and add-frames exceptions for lunging and wing threat
    if strcmp(action,'lunging'), min_len = 1; if add_frames == 0, add_frames = 1; end; end
    if strcmp(action,'wingthreat'), min_len = 1; if add_frames == 0, add_frames = 5; end; end

    % PREPARE A "fly_feat.action" VECTOR WITH ALL DETECTED BOUTS 
    % OF THE ACTION OF INTEREST
    ngens = max(genotypes);
    add_path_base = add_path; add_path = [add_path '/' add_path '_' action];
    if ~ispc, slash = '/'; else slash = '\'; end
    % FOR ALL GENOTYPES
    igen = 1; gen = '';
    while igen <= ngens,
        gen_ind = find(genotypes == igen); nfiles = length(gen_ind);
        if iscell(gen_ident), gen = gen_ident{igen}; end
        
        gen = gen(setdiff(1:length(gen),strfind(gen,'_')));
        % Which was the latest processed genotype?
        load([path 'analysis' slash gen '_aggress_court.mat']);
        
        if strcmp(action,'wingthreat'),
            num1 = [0 cumsum(wing.threat.obj1.number)]+1;
            num2 = [0 cumsum(wing.threat.obj2.number)]+1;
            ind1 = [0 cumsum(wing.threat.obj1.len)];
            ind2 = [0 cumsum(wing.threat.obj2.len)];
        elseif strcmp(action,'tussl'),
            num1 = [0 cumsum(tussl.number)]+1;
            ind1 = [0 cumsum(tussl.len)];
        elseif strcmp(action,'onewing'),
            num1r = [0 cumsum(wing.ext.r.obj1.number)]+1;
            num2r = [0 cumsum(wing.ext.r.obj2.number)]+1;
            num1l = [0 cumsum(wing.ext.l.obj1.number)]+1;
            num2l = [0 cumsum(wing.ext.l.obj2.number)]+1;
            ind1r = [0 cumsum(wing.ext.r.obj1.len)];
            ind2r = [0 cumsum(wing.ext.r.obj2.len)];
            ind1l = [0 cumsum(wing.ext.l.obj1.len)];
            ind2l = [0 cumsum(wing.ext.l.obj2.len)];
        elseif strcmp(action,'circling'),
            num1 = [0 cumsum(court.obj1.number)]+1;
            num2 = [0 cumsum(court.obj2.number)]+1;
            ind1 = [0 cumsum(court.obj1.len)];
            ind2 = [0 cumsum(court.obj2.len)];
        elseif strcmp(action,'chasing'),
            num1 = [0 cumsum(chase.obj1.number)]+1;
            num2 = [0 cumsum(chase.obj2.number)]+1;
            ind1 = [0 cumsum(chase.obj1.len)];
            ind2 = [0 cumsum(chase.obj2.len)];
        end

        % PROCESS EACH MOVIE
        ifile = 1;
        while ifile <= nfiles,
            
            % Keep track of the already processed movie clips
            fid = fopen([path add_path slash 'sav_indxg.mat'],'r');
            if (fid >= 0) && (~params.init),
                % Which was the latest processed genotype, movie, clip index?
                load([path add_path slash 'sav_indxg.mat']);
                gen_ind = find(genotypes == igen);
            else
                % Reset
                indxg = ones(1,ngens);                
                [s, mess, messid] = mkdir([path add_path_base]);
                [s, mess, messid] = mkdir([path add_path]);
                if isempty(gen),
                    if ~isempty(dir([path add_path slash 'cluster_0*'])),
                        rmdir([path add_path slash 'cluster_0*'],'s');
                    end
                else
                    if ~isempty(dir([path add_path slash gen])),
                        rmdir([path add_path slash gen],'s');
                    end
                end
                params.init = 0;
            end
            if (fid >= 0), fclose(fid); end
            if iscell(sorted_names), FileN = sorted_names{gen_ind(ifile)}; else FileN = sorted_names; end

            % Read some addtional features like positions
            fid = fopen([FileN(1:end-5) '_feat.mat'],'r');
            if (fid < 0) || params.feat_read_new,
                % Process feature file
                [fly_feat,NFrms] = read_feat(FileN,params);
                save([FileN(1:end-5) '_feat.mat'],'fly_feat','NFrms');
            else
                % Load feature file
                load([FileN(1:end-5) '_feat.mat']);
                fclose(fid);
            end
            
            % EXTRACT BOUTS OF ACTIONS
            if strcmp(action,'lunging'),
                if (params.courtship && cop.pre(ifile) > 0),
                    I_GOOD = find(fly_feat.frame < cop.pre(ifile) | fly_feat.frame > cop.post(ifile));
                else
                    I_GOOD = find(fly_feat.frame <= params.max_frames/dt(ifile));
                end
                tmp = fly_feat.obj1.lunge(I_GOOD) + fly_feat.obj2.lunge(I_GOOD);
                tmp = tmp > 0; tmp = bwlabel(tmp);
                eval(['fly_feat.' action ' = tmp;']);
            elseif strcmp(action,'wingthreat'),
                tmp = zeros(1,NFrms);
                for ii=num1(ifile):num1(ifile+1)-1,
                    tmp(wing.threat.obj1.lab(ind1(ii)+1):wing.threat.obj1.lab(ind1(ii+1))) = 1;
                end
                for ii=num2(ifile):num2(ifile+1)-1,
                    tmp(wing.threat.obj2.lab(ind2(ii)+1):wing.threat.obj2.lab(ind2(ii+1))) = 1;
                end
                tmp = tmp > 0; tmp = bwlabel(tmp);
                eval(['fly_feat.' action ' = tmp;']);
            elseif strcmp(action,'tussl'),
                tmp = zeros(1,NFrms);
                for ii=num1(ifile):num1(ifile+1)-1,
                    tmp(tussl.ind(ind1(ii)+1):tussl.ind(ind1(ii+1))) = 1;
                end
                tmp = tmp > 0; tmp = bwlabel(tmp);
                eval(['fly_feat.' action ' = tmp;']);
            elseif strcmp(action,'onewing'),
                tmp = zeros(1,NFrms);
                for ii=num1r(ifile):num1r(ifile+1)-1,
                    tmp(wing.ext.r.obj1.lab(ind1r(ii)+1):wing.ext.r.obj1.lab(ind1r(ii+1))) = 1;
                end
                for ii=num1l(ifile):num1l(ifile+1)-1,
                    tmp(wing.ext.l.obj1.lab(ind1l(ii)+1):wing.ext.l.obj1.lab(ind1l(ii+1))) = 1;
                end
                for ii=num2r(ifile):num2r(ifile+1)-1,
                    tmp(wing.ext.r.obj2.lab(ind2r(ii)+1):wing.ext.r.obj2.lab(ind2r(ii+1))) = 1;
                end
                for ii=num2l(ifile):num2l(ifile+1)-1,
                    tmp(wing.ext.l.obj2.lab(ind2l(ii)+1):wing.ext.l.obj2.lab(ind2l(ii+1))) = 1;
                end
                tmp = tmp > 0; tmp = bwlabel(tmp);
                eval(['fly_feat.' action ' = tmp;']);
            elseif strcmp(action,'circling'),
                tmp = zeros(1,NFrms);
                for ii=num1(ifile):num1(ifile+1)-1,
                    tmp(court.obj1.lab(ind1(ii)+1):court.obj1.lab(ind1(ii+1))) = 1;
                end
                for ii=num2(ifile):num2(ifile+1)-1,
                    tmp(court.obj2.lab(ind2(ii)+1):court.obj2.lab(ind2(ii+1))) = 1;
                end
                tmp = tmp > 0; tmp = bwlabel(tmp);
                eval(['fly_feat.' action ' = tmp;']);
            elseif strcmp(action,'chasing'),
                tmp = zeros(1,NFrms);
                for ii=num1(ifile):num1(ifile+1)-1,
                    tmp(chase.obj1.lab(ind1(ii)+1):chase.obj1.lab(ind1(ii+1))) = 1;
                end
                for ii=num2(ifile):num2(ifile+1)-1,
                    tmp(chase.obj2.lab(ind2(ii)+1):chase.obj2.lab(ind2(ii+1))) = 1;
                end
                tmp = tmp > 0; tmp = bwlabel(tmp);
                eval(['fly_feat.' action ' = tmp;']);
            elseif strcmp(action,'copulation'),
                dt = median(fly_feat(1).time(2:NFrms(1)) - fly_feat(1).time(1:NFrms(1)-1));
                [pre,post] = detect_copul(fly_feat,dt);
                tmp = fly_feat.obj1.chase * 0; icnt = 0;
                for ii=1:numel(pre),
                    if (pre(ii) > 0),
                        icnt = icnt + 1;
                        if post(ii) < 0, post(ii) = length(tmp); end
                        pre = find(fly_feat.frame == pre);
                        post = find(fly_feat.frame == post);
                        tmp(pre(ii):post(ii)) = icnt;
                    end
                end
                eval(['fly_feat.' action ' = tmp;']);
            end
            nevents = max(tmp);

            % CALL "PLOT_EVENTS" FUNCTION
            if nevents > 0,
                plot_events(fly_feat,action,nevents,igen,ngens,gen,ifile,nfiles,min_len,FileN,path,add_path,add_frames,params);
            end
            
            ifile = ifile + 1;
            % Save last position (movie, clip index)
            save([path add_path slash 'sav_indxg.mat'],'indxg','igen','ifile');
        end
        ifile = 1; igen = igen + 1;
        % Save last position (movie, clip index)
        save([path add_path slash 'sav_indxg.mat'],'indxg','igen','ifile');
    end
    % COMPUTE A WEBSITE (HTML) TO WATCH THE ACTION CLIPS
    if params.plots.showhtml,
        write_html(path);
    end
end

end

