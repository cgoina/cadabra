%% Analyze_data
% Copyright (C) 2009 Heiko Dankert, California Institute of Technology
%                    Jinyang Liu, Howard Huges Mecial Institute

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
% * Implementation by Heiko Dankert, Jinyang Liu
%
%% Analyzed data provided by 'read_feat', detect behaviors, store results

function analyze_data(path,name,params)

% PARAMETERS
% min_b_len = minimum bout length [s]
% gap = time between two bouts
% Action bouts with gaps < 'max_gap' are considered as one bout
% Wing extensions

%JL0831 add the global variables ind1_count and ind2_count to
%output of the aggress_court.mat.
global ind1_count ind2_count;

    min_b_len.long = 1; max_gap.long = .6; % [s] (wing extension)
    min_b_len.short = 0.01; max_gap.short = 0.03; % [s] (wing flick!)
    % Circling
    min_b_len_circ = .69; max_gap_circ = .6; % [s]
    % Chasing
    min_b_len_chase = 1; max_gap_chase = 0.6; % [s]
    % Wing threats
    b_len_threat.min = .3; b_len_threat.max = 30;
    max_gap_threat = .5; % [s]

    % ANALYZE DATA
    % DECLARE VARIABLES USED FOR PLOT FUNCTIONS 
    % CONSISTING OF DATA FOR ALL GENOTYPES
    ngens = 1;
    nmovs = zeros(ngens,1); dts = cell(ngens,1); nframes = dts;
    obj1 = cell(ngens,1); obj2 = cell(ngens,1); 
    lunges = cell(ngens,1); tussls = cell(ngens,1);
    jumps = cell(ngens,1); chases = cell(ngens,1); charges = cell(ngens,1); 
    courts = cell(ngens,1); 
    wings.ext.r = cell(ngens,1); wings.ext.l = cell(ngens,1); wings.ext.b = cell(ngens,1);
    wings.fli.r = cell(ngens,1); wings.fli.l = cell(ngens,1); wings.fli.b = cell(ngens,1);
    wings.threat = cell(ngens,1);
    wings.min_bout.long = 0; wings.min_bout.short = 0; wings.max_gap.long = 0; wings.max_gap.short = 0;
    wings.bout_threat = 0; wings.max_gap_threat = 0;
    copu.pre = dts; copu.post = dts; copu.int = dts;
    dists = cell(ngens,1); %ddists = cell(ngens,1);
    proxi.times = cell(ngens,1); proxi.indi = cell(ngens,1);
    addname = ['analysis' params.slash];
    [s,mess,messid] = mkdir([path addname]);
    
    fid = fopen([path addname name '_analysis.mat'],'r');
    if (fid < 0) || params.analyze_new
        Y1.t = []; Y1.x = []; Y1.y = [];
        Y1.vx = []; Y1.vy = []; Y1.v = []; Y1.vs = []; Y1.do = [];
        Y1.am = []; Y1.as = []; Y1.lm = []; Y1.ls = [];
        Y2 = Y1;
        lunge.lab = []; lunge.mov = []; lunge.obj = []; lunge.t = []; lunge.number = [];
        lunge.obj1.x1 = []; lunge.obj1.y1 = []; lunge.obj1.x2 = []; lunge.obj1.y2 = [];
        lunge.obj1.t = []; lunge.obj1.len = []; lunge.obj1.number = [];
        lunge.obj2 = lunge.obj1;

        wing.ext.r = []; wing.ext.l = []; wing.ext.b = [];
        wing.fli.r = []; wing.fli.l = []; wing.fli.b = [];
        wing.threat = [];
        jump = []; charge = []; chase = []; court = []; tussl = [];
        dist_c = []; ddist_c = [];
        dist_h = []; dist_t = []; dist_h1t2 = []; dist_h2t1 = [];
        food_dist.obj1 = []; food_dist.obj2 = [];
        food_time.obj1.near = []; food_time.obj1.mid = []; food_time.obj1.far = [];
        food_time.obj2.near = []; food_time.obj2.mid = []; food_time.obj2.far = [];
        proxi_time.near = []; proxi_time.mid = []; proxi_time.far = [];
        proxi_ind = proxi_time;
        proxi_vel.obj1.near = []; proxi_vel.obj1.mid = []; proxi_vel.obj1.far = [];
        proxi_vel.obj2.near = []; proxi_vel.obj2.mid = []; proxi_vel.obj2.far = [];
        proxi_mvdirdiff = proxi_vel; movdist.obj1 = []; movdist.obj2 = [];
        movdirdiff.obj1 = []; movdirdiff.obj2 = [];
        NFrms = []; fly_feat = []; dt = []; nfrms = []; cop.pre = []; cop.post = [];
        nmov = params.nchambers;
        ifile = 0;
        
        for imov=1:nmov
            % READ MOVIE FEATURE FILE
            %if iscell(name), FileN = name{gen_ind(imov)}; else FileN = name; end
            FileN = [path name '_' num2str(imov) '.feat'];

            fid = fopen(FileN,'r');
            if (fid < 0) || params.analyze_new
                % READ OUT FEATURE TEXT FILE
                [fly_feat,NFrms] = read_feat(FileN,params);
                save([FileN(1:end-5) '_feat.mat'],'fly_feat','NFrms');
            else
                % READ BINARY
                load([FileN '_feat.mat']);
                fclose(fid);
            end
            
            % READ ROI AND SCALE DATA
            [ROI,scale] = load_ROI(FileN);
            width.x = (ROI.cols(end)-ROI.cols(1))*scale.x;
            width.y = (ROI.rows(end)-ROI.rows(1))*scale.y;
            if (width.x > width.y), tmp = width.x; width.x = width.y; width.y = tmp; end
            
            % CENTER THE ROI COORDINATES
            fly_feat.obj1.pos_x = fly_feat.obj1.pos_x - width.x/2;
            fly_feat.obj1.pos_y = fly_feat.obj1.pos_y - width.y/2;
            fly_feat.obj2.pos_x = fly_feat.obj2.pos_x - width.x/2;
            fly_feat.obj2.pos_y = fly_feat.obj2.pos_y - width.y/2;

            % TIME DIFFERENCE BETWEEN TWO FRAMES
            dtf = fly_feat(1).time(2:NFrms(1)) - fly_feat(1).time(1:NFrms(1)-1);

            % Test whether both flies are within the borders
            test = den_feat(fly_feat,0.01,0,median(dtf),params,width);
            % Required min. data length after filtering
            min_data_len = 1000; %[frames]
            
            % FILE NOT EMPTY?
            if numel(test.time) >= min_data_len && sum(fly_feat.obj1.pos_change+fly_feat.obj2.pos_change) > 1,
                ifile = ifile + 1;
                
                dt(ifile) = median(dtf); %#ok<AGROW>   
                                                
                % JUMP DETECTION (in preparation)
                max_gap_jump = 0.6; %[s]
                jump = jumping(fly_feat,max_gap_jump,jump);
                
                fly_feat.lab = 1:NFrms;
                
                % PREOVIDE SMOOTHED VELOCITY VECTORS FOR VELOCITY TIME
                % SERIES PLOTS
                window_size = 5; % seconds
                bin_size = ceil(window_size/dt(ifile));

                fly_feat.obj1.vels = filtfilt(ones(1,bin_size)/bin_size,1,fly_feat.obj1.vel);
                fly_feat.obj2.vels = filtfilt(ones(1,bin_size)/bin_size,1,fly_feat.obj2.vel);
                
                % DENOISE DATA (REMOVE JUMPS, HIGH VELOCITIES)
                [fly_feat, NFrms] = den_feat(fly_feat,0.01,0,dt(ifile),params,width); % features, quantile, margin
                nfrms(ifile) = NFrms; %#ok<AGROW>

                % Horizontal adjustment of ROI-position
                diff = min(min([fly_feat.obj1.pos_x fly_feat.obj2.pos_x])) - 0.5; 
                if  diff > 0.1,
                    fly_feat.obj1.pos_x = fly_feat.obj1.pos_x - diff; fly_feat.obj2.pos_x = fly_feat.obj2.pos_x - diff;
                end

                % COPULATION DETECTION
                if params.courtship && ~params.oneobj,
                    % pre = start frame; post = end frame
                    [cop.pre,cop.post] = detect_copul(fly_feat,dt(ifile),cop.pre,cop.post);
                else
                    cop.pre = [cop.pre 0]; cop.post = [cop.post NFrms];
                end

                % TIME LIMIT (in case of data-sets with different lengths
                if fly_feat.frame(end) > params.max_frames || (cop.pre(ifile) > 0),
                    if (params.courtship && cop.pre(ifile) > 0),
                        I_GOOD = find((fly_feat.frame < cop.pre(ifile) | fly_feat.frame > cop.post(ifile)) & ...
                                       fly_feat.frame <= params.max_frames); %#ok<NASGU>
                    else
                        I_GOOD = find(fly_feat.frame <= params.max_frames); %#ok<NASGU>
                    end
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

                % Position change from frame to frame, filter out constant
                % values (flies leaving the arena)
                dpos_ch = conv2(fly_feat.obj1.pos_change,[-1 0 1],'same');
                ind = find((dpos_ch == 0) | (fly_feat.obj1.pos_change > 2));
                if numel(ind), fly_feat.obj1.pos_change(ind) = 0; end
                dpos_ch = conv2(fly_feat.obj2.pos_change,[-1 0 1],'same');
                ind = find((dpos_ch == 0) | (fly_feat.obj2.pos_change > 2));
                if numel(ind), fly_feat.obj2.pos_change(ind) = 0; end

                % WALKING DISTANCE
                movd1 = single(cumsum(fly_feat.obj1.pos_change));
                movd2 = single(cumsum(fly_feat.obj2.pos_change));
                
                % Save only every 2nd value a higher frame rates 
                % to safe memory
                if dt(ifile) < 1/15,
                    indsc = 1:2:length(fly_feat.obj1.pos_x);
                else
                    indsc = 1:length(fly_feat.obj1.pos_x);
                end
                
                % FLY POSITIONS, VELOCITIES, TIME
                x1 = single(fly_feat.obj1.pos_x(indsc)); y1 = single(fly_feat.obj1.pos_y(indsc));
                x2 = single(fly_feat.obj2.pos_x(indsc)); y2 = single(fly_feat.obj2.pos_y(indsc));
                Y1.x = [Y1.x x1]; Y1.y = [Y1.y y1];
                Y2.x = [Y2.x x2]; Y2.y = [Y2.y y2];
                Y1.vx = [Y1.vx single(fly_feat.obj1.vel(indsc) .* cos(fly_feat.obj1.headdir(indsc)*pi/180))];
                Y1.vy = [Y1.vy single(fly_feat.obj1.vel(indsc) .* sin(fly_feat.obj1.headdir(indsc)*pi/180))];
                Y2.vx = [Y2.vx single(fly_feat.obj2.vel(indsc) .* cos(fly_feat.obj2.headdir(indsc)*pi/180))];
                Y2.vy = [Y2.vy single(fly_feat.obj2.vel(indsc) .* sin(fly_feat.obj2.headdir(indsc)*pi/180))];
                Y1.t = [Y1.t single(fly_feat.time(indsc))]; Y2.t = [Y2.t single(fly_feat.time(indsc))];
                Y1.v = [Y1.v single(fly_feat.obj1.vel(indsc))]; Y2.v = [Y2.v single(fly_feat.obj2.vel(indsc))];
                Y1.vs = [Y1.vs single(fly_feat.obj1.vels(indsc))]; Y2.vs = [Y2.vs single(fly_feat.obj2.vels(indsc))];
                Y1.do = [Y1.do single(fly_feat.obj1.headdir(indsc))]; Y2.do = [Y2.do single(fly_feat.obj2.headdir(indsc))];

                % FLY BODY AREA, LENGTH
                if mean(fly_feat.obj1.FArea) > 0,
                    a = single(fly_feat.obj1.FArea); l = single(fly_feat.obj1.FLength); l(l == 0) = 0.0001;
                    ind = find((l > quantile(l,.2)) & (l < quantile(l,.9)) & ...
                        (a > quantile(a,.2)) & (a./l.^2 > .1)); %a/l^2 > .25
                    Y1.am = [Y1.am mean(a(ind))]; Y1.as = [Y1.as std(a(ind))];
                    Y1.lm = [Y1.lm mean(l(ind))]; Y1.ls = [Y1.ls std(l(ind))];
                    a = fly_feat.obj2.FArea; l = fly_feat.obj2.FLength; l(l == 0) = 0.0001;
                    ind = find((l > quantile(l,.2)) & (l < quantile(l,.9)) & ...
                        (a > quantile(a,.2)) & (a./l.^2 > .1)); %a/l^2 > .25
                    Y2.am = [Y2.am mean(a(ind))]; Y2.as = [Y2.as std(a(ind))];
                    Y2.lm = [Y2.lm mean(l(ind))]; Y2.ls = [Y2.ls std(l(ind))];
                else
                    Y1.am = [Y1.am 0]; Y2.am = [Y2.am 0]; Y1.as = [Y1.as 0]; Y2.as = [Y2.as 0];
                    Y1.lm = [Y1.lm 0]; Y2.lm = [Y2.lm 0]; Y1.ls = [Y1.ls 0]; Y2.ls = [Y2.ls 0];
                end

                % RELATIVE FLY DISTANCE AND CHANGE OF DISTANCE
                % c = center, h = head, t = abdomen
                dist_c = [dist_c single(fly_feat.distc(indsc))]; %#ok<AGROW>
                ddist_c = [ddist_c single(fly_feat.der_distc(indsc))]; %#ok<AGROW>
                dist_h = [dist_h single(fly_feat.disth(indsc))]; %#ok<AGROW>
                dist_t = [dist_t single(fly_feat.distt(indsc))]; %#ok<AGROW>
                dist_h1t2 = [dist_h1t2 single(fly_feat.disth1t2(indsc))]; %#ok<AGROW>
                dist_h2t1 = [dist_h2t1 single(fly_feat.disth2t1(indsc))]; %#ok<AGROW>

                % FLY VELOCITIES AT DIFFERENT RELATIVE DISTANCES
                % Near (< 5 mm)
                ind = find(fly_feat.distc < 5); a = single(fly_feat.time(ind));
                proxi_time.near = [proxi_time.near (numel(a)-1) * dt(ifile)];
                if numel(ind) > 1000, ind = ind(1:2:numel(ind)); end % try saving some memory
                proxi_ind.near = [proxi_ind.near numel(a(1:2:numel(ind)))];
                proxi_vel.obj1.near = [proxi_vel.obj1.near single(fly_feat.obj1.vels(ind))];
                proxi_vel.obj2.near = [proxi_vel.obj2.near single(fly_feat.obj2.vels(ind))];
                proxi_mvdirdiff.obj1.near = [proxi_mvdirdiff.obj1.near single(fly_feat.obj1to2mvdirdiff(ind))];
                proxi_mvdirdiff.obj2.near = [proxi_mvdirdiff.obj2.near single(fly_feat.obj2to1mvdirdiff(ind))];
                % Mid-range (5 mm - 10 mm)
                ind = find((fly_feat.distc >= 5) & (fly_feat.distc < 10)); a = single(fly_feat.time(ind));
                proxi_time.mid = [proxi_time.mid (numel(a)-1) * dt(ifile)];
                if numel(ind) > 1000, ind = ind(1:2:numel(ind)); end % try saving some memory
                proxi_ind.mid = [proxi_ind.mid numel(a(1:2:numel(ind)))];
                proxi_vel.obj1.mid = [proxi_vel.obj1.mid single(fly_feat.obj1.vels(ind))];
                proxi_vel.obj2.mid = [proxi_vel.obj2.mid single(fly_feat.obj2.vels(ind))];
                proxi_mvdirdiff.obj1.mid = [proxi_mvdirdiff.obj1.mid single(fly_feat.obj1to2mvdirdiff(ind))];
                proxi_mvdirdiff.obj2.mid = [proxi_mvdirdiff.obj2.mid single(fly_feat.obj2to1mvdirdiff(ind))];
                % Far (>= 10 mm)
                ind = find(fly_feat.distc >= 10); a = single(fly_feat.time(ind));
                proxi_time.far = [proxi_time.far (numel(a)-1) * dt(ifile)];
                if numel(ind) > 1000, ind = ind(1:2:numel(ind)); end % try saving some memory
                proxi_ind.far = [proxi_ind.far numel(a(1:2:numel(ind)))];
                proxi_vel.obj1.far = [proxi_vel.obj1.far single(fly_feat.obj1.vels(ind))];
                proxi_vel.obj2.far = [proxi_vel.obj2.far single(fly_feat.obj2.vels(ind))];
                proxi_mvdirdiff.obj1.far = [proxi_mvdirdiff.obj1.far single(fly_feat.obj1to2mvdirdiff(ind))];
                proxi_mvdirdiff.obj2.far = [proxi_mvdirdiff.obj2.far single(fly_feat.obj2to1mvdirdiff(ind))];

                % TIME SPEND ON/OFF FOOD FLY 1
                % Near (< 5 mm)
                ind = logical(fly_feat.obj1.r < 5); a = single(fly_feat.time(ind));
                food_time.obj1.near = [food_time.obj1.near (numel(a)-1) * dt(ifile)];
                % Mid-range (5 mm - 10 mm)
                ind = logical((fly_feat.obj1.r >= 5) & (fly_feat.obj1.r < 10)); a = single(fly_feat.time(ind));
                food_time.obj1.mid = [food_time.obj1.mid (numel(a)-1) * dt(ifile)];
                % Far (>= 10 mm)
                ind = logical(fly_feat.obj1.r >= 10); a = single(fly_feat.time(ind));
                food_time.obj1.far = [food_time.obj1.far (numel(a)-1) * dt(ifile)];

                % TIME SPEND ON/OFF FOOD FLY 2
                % Near (< 5 mm)
                ind = logical(fly_feat.obj2.r < 5); a = single(fly_feat.time(ind));
                food_time.obj2.near = [food_time.obj2.near (numel(a)-1) * dt(ifile)];
                % Mid-range (5 mm - 10 mm)
                ind = logical((fly_feat.obj2.r >= 5) & (fly_feat.obj2.r < 10)); a = single(fly_feat.time(ind));
                food_time.obj2.mid = [food_time.obj2.mid (numel(a)-1) * dt(ifile)];
                % Far (>= 10 mm)
                ind = logical(fly_feat.obj2.r >= 10); a = single(fly_feat.time(ind));
                food_time.obj2.far = [food_time.obj2.far (numel(a)-1) * dt(ifile)];

                % DISTANCE FROM FOOD
                food_dist.obj1 = [food_dist.obj1 single(fly_feat.obj1.r(indsc))];
                food_dist.obj2 = [food_dist.obj2 single(fly_feat.obj2.r(indsc))];

                % DIFFERENCE BETWEEN MOVING DIRECTION OF FLY X AND THE
                % VECTOR FROM FLY X TO Y [deg]
                movdirdiff.obj1 = [movdirdiff.obj1 single(fly_feat.obj1to2mvdirdiff(indsc))];
                movdirdiff.obj2 = [movdirdiff.obj2 single(fly_feat.obj2to1mvdirdiff(indsc))];
                
                % WALKING DISTANCE (cumulative)
                movdist.obj1 = [movdist.obj1 movd1(indsc)];
                movdist.obj2 = [movdist.obj2 movd2(indsc)];                                

                % LUNGES (detection in 'read_feat')
                ind = zeros(1,NFrms); ind(logical(fly_feat.obj1.lunge)) = 1;
                ind(logical(fly_feat.obj2.lunge)) = 2; ind1 = ind(logical(ind>0)); ind = find(ind>0);
                lunge.mov(ifile) = ifile;
                lunge.obj1.number = [lunge.obj1.number numel(find(fly_feat.obj1.lunge))];
                lunge.obj2.number = [lunge.obj2.number numel(find(fly_feat.obj2.lunge))];
                lunge.lab = [lunge.lab fly_feat.lab(ind)]; lunge.obj = [lunge.obj ind1];
                lunge.t = [lunge.t fly_feat.time(ind)];
                ind1 = find(fly_feat.obj1.lunge); ind2 = find(fly_feat.obj2.lunge);
                lunge.number = [lunge.number numel(fly_feat.time(ind))];
                lunge.obj1.x1 = [lunge.obj1.x1 fly_feat.obj1.pos_x(ind1)];
                lunge.obj1.y1 = [lunge.obj1.y1 fly_feat.obj1.pos_y(ind1)];
                lunge.obj1.x2 = [lunge.obj1.x2 fly_feat.obj2.pos_x(ind1)];
                lunge.obj1.y2 = [lunge.obj1.y2 fly_feat.obj2.pos_y(ind1)];
                lunge.obj1.t = [lunge.obj1.t fly_feat.time(ind1)];
                lunge.obj2.x1 = [lunge.obj2.x1 fly_feat.obj2.pos_x(ind2)];
                lunge.obj2.y1 = [lunge.obj2.y1 fly_feat.obj2.pos_y(ind2)];
                lunge.obj2.x2 = [lunge.obj2.x2 fly_feat.obj1.pos_x(ind2)];
                lunge.obj2.y2 = [lunge.obj2.y2 fly_feat.obj1.pos_y(ind2)];
                lunge.obj2.t = [lunge.obj2.t fly_feat.time(ind2)];

                % WING EXTENSIONS
                wing = wing_ext(fly_feat,params,FileN,min_b_len,max_gap,wing);

                % CHASING (detection in 'read_feat')
                chase = chasing(fly_feat,min_b_len_chase,max_gap_chase,chase);

                % CIRCLING (detection in 'read_feat')
                court = circling(fly_feat,min_b_len_circ,max_gap_circ,court);

                % TUSSLING (detection in 'read_feat')
                [fly_feat,tussl] = detect_tussls(fly_feat,params,FileN,ifile,tussl);
                
                % WING THREATS (detection in 'read_feat')
                wing.threat = detect_threat(fly_feat,b_len_threat,max_gap_threat,wing.threat);
            end
        end
        lunge.obj1.len = ones(1,sum(lunge.obj1.number));
        lunge.obj2.len = ones(1,sum(lunge.obj2.number));

        %JL0831 save the global variables ind1_count and ind2_count, probable lunging, to the output of the aggress_court.mat.
        ind1_count = find(ind1_count); 
        ind2_count = find(ind2_count);
        
        % SAVE DATA INTO MATLAB BINARY FILES
        save([path  addname name  '_aggress_court.mat'],...
            'lunge','jump','charge','chase','court','tussl','wing',...
            'cop', 'min_b_len','max_gap','min_b_len_circ','max_gap_circ',...
            'b_len_threat','max_gap_threat','dt','nfrms','ifile', ...
            'ind1_count', 'ind2_count');
        save([path  addname name  '_analysis.mat'],'Y1','Y2',...
            'dist_c','ddist_c','dist_h','dist_t','dist_h1t2','dist_h2t1',...
            'food_dist','food_time','movdist','movdirdiff','proxi_time',...
            'proxi_ind','proxi_vel','proxi_mvdirdiff','nmov','dt','nfrms','ifile');
    else
        % READ IF DATA IS PRESENT
        fid2 = fopen([path  addname name  '_aggress_court.mat'],'r');
        load([path  addname name  '_aggress_court.mat']);
        load([path  addname name  '_analysis.mat']);
        fclose(fid); fclose(fid2);
    end


end