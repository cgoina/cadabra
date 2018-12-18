%% Analyze_distances
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
%% Analyze distances of flies

function FID = analyze_distances(obj1,obj2,lunges,copu,dists,proxi,nmovs,dts,gen_ident,gen_ident1,params,chamber,FID)

if params.plots.stat.dist,
    % FLY DISTANCE FROM FOOD AREA
    % utest = utestpairs(gen_ident1);
    ngens = length(obj1);
    x = 1:ngens; y.mea = zeros(ngens,2); y.sem = y.mea; y.fmin = y.mea;
    y.fmax = y.mea; y.min = y.mea; y.max = y.mea;
    X = []; G = []; X2 = []; G2 = [];
    % Prepare data
    data = win_los_data('food_dist','',obj1,obj2,lunges,copu.pre,dts,params);
    for igen=1:ngens,
        % Statistics for each genotype
        % For barplots
        y.mea(igen,1) = mean(data{igen}.win);
        y.sem(igen,1) = std(data{igen}.win)/sqrt(nmovs(igen));
        y.fmin(igen,1) = - quantile(data{igen}.win,0.25) + y.mea(igen,1);
        y.fmax(igen,1) = quantile(data{igen}.win,0.75) - y.mea(igen,1);
        y.min(igen,1) = min(data{igen}.win); y.max(igen,1) = max(data{igen}.win);

        if params.oneobj, data{igen}.los(1,:) = 0; end
        y.mea(igen,2) = mean(data{igen}.los);
        y.sem(igen,2) = std(data{igen}.los)/sqrt(nmovs(igen));
        y.fmin(igen,2) = - quantile(data{igen}.los,0.25) + y.mea(igen,2);
        y.fmax(igen,2) = quantile(data{igen}.los,0.75) - y.mea(igen,2);
        y.min(igen,2) = min(data{igen}.los); y.max(igen,2) = max(data{igen}.los);

        % For boxplots
        X2 = [X2 data{igen}.win data{igen}.los];
        G2 = [G2 ones(1,length(data{igen}.win))*(2*igen-1) ones(1,length(data{igen}.los))*2*igen];
    end
    % Plot distance from food area
    if numel(find(params.plots.stat.flydist_vec == 1)),
        maxy = min(chamber.width,chamber.height);
        ytit = 'food distance [mm]'; titl = 'Distance from Food Area';
        utest = utestpairs(gen_ident1);
        FID = plot_bar_msem(y,0,X2,G2,maxy,ytit,titl,FID,gen_ident,params);
    end

    % RELATIVE FLY DISTANCE
    if numel(find(params.plots.stat.flydist_vec == 2)),
        X = []; G = [];
        x = 1:ngens; y.mea = zeros(ngens,1); y.sem = y.mea; y.fmin = y.mea; y.fmax = y.mea; y.min = y.mea; y.max = y.mea;
        for igen=1:ngens,
            % In case of one fly relate distance to center of arena
            if params.oneobj,
                r = sqrt(obj1{igen}.x.^2 + obj1{igen}.y.^2);
            else
                r = dists{igen}.distc;
            end
            % Statistics for each genotype
            % For barplots
            y.mea(igen) = mean(r);
            y.sem(igen) = std(r) / sqrt(nmovs(igen));
            y.fmin(igen) = - quantile(r,0.25) + y.mea(igen);
            y.fmax(igen) = quantile(r,0.75) - y.mea(igen);
            y.min(igen) = min(r); y.max(igen) = max(r);
            % For boxplots
            X = [X r]; G = [G ones(1,numel(r))*igen];
        end
        % Occurrence or time spent?
        if params.tim,
            maxy = params.range.stat.tim.flydist;
        else
            maxy = params.range.stat.occ.flydist;
        end
        % Plot data
        ytit = 'rel. fly distance [mm]'; titl = 'Distance between Flies';
        FID = plot_bar_msem(y,0,X,G,maxy,ytit,titl,FID,gen_ident,params);
    end

    % FLY PROXIMITY PERIOD
    if numel(find(params.plots.stat.flydist_vec == 3)),
        % TIME FLIES SPENT IN DIFFERENT RANGES FROM FOOD
        % FOOD LOCATED IN ARENA CENTER
        x = 1:ngens; y1.mea = zeros(ngens,1); y1.fmin = zeros(ngens,1); y1.fmax = zeros(ngens,1); y2 = y1; y3 = y1;
        for igen=1:ngens,
            % Statistics (mean, sem, min, max)
            y1.mea(igen) = mean(proxi.times{igen}.near);
            y1.fmin(igen) = - quantile(proxi.times{igen}.near,0.25) + y1.mea(igen);
            y1.fmax(igen) = quantile(proxi.times{igen}.near,0.75) - y1.mea(igen);
            y1.min(igen) = min(proxi.times{igen}.near); y1.max(igen) = max(proxi.times{igen}.near);
            y2.mea(igen) = mean(proxi.times{igen}.mid);
            y2.fmin(igen) = - quantile(proxi.times{igen}.mid,0.25) + y2.mea(igen);
            y2.fmax(igen) = quantile(proxi.times{igen}.mid,0.75) - y2.mea(igen);
            y2.min(igen) = min(proxi.times{igen}.mid); y2.max(igen) = max(proxi.times{igen}.mid);
            y3.mea(igen) = mean(proxi.times{igen}.far);
            y3.fmin(igen) = - quantile(proxi.times{igen}.far,0.25) + y3.mea(igen);
            y3.fmax(igen) = quantile(proxi.times{igen}.far,0.75) - y3.mea(igen);
            y3.min(igen) = min(proxi.times{igen}.far); y3.max(igen) = max(proxi.times{igen}.far);
        end
        % PLOT DATA AS NEAR-MID-FAR RANGE STACK FOR EACH GENOTYPE
        if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
        figure(FID); clf;
        set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
        if ngens > 1,
            h = bar(x,[y1.mea(:,1) y2.mea(:,1) y3.mea(:,1)],params.barwidth,'stack');
            set_xtick_label(gen_ident,45,params.cross_to,params.axisfontsize);
        else
            dat = [y1.mea(:,1) y2.mea(:,1) y3.mea(:,1)];
            h = bar([dat;nan(size(dat))],params.barwidth,'stack'); xlim([.5 1.5]);
            set(gca,'XTickLabel',gen_ident);
        end
        set(gca,'XTick',1:ngens);
        colormap([params.flycol.win ; params.flycol.los ; params.flycol.ext]); hold on;
        t = 0:300:params.max_frames; ytick_lab = cell(1,length(t));
        for i=1:length(t), ytick_lab{i} = num2str(t(i)/60); end
        set(gca,'YTick',t); set(gca,'YTickLabel',ytick_lab); set(gca,'FontSize',params.axisfontsize);
        ylabel('fly proximity time [min] (mean)');
        axis([0 ngens+1 0 params.max_frames]);
        legend('< 5 mm','5 - 10 mm','> 10 mm'); legend('boxoff');
        titl = 'Fly Proximity Period';
        set(FID,'Name',titl); mtit(titl,'FontSize',params.axisfontsize+2,'yoff',.04);
        if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
    end
    
    % TIME SPENT ON FOOD
    if numel(find(params.plots.stat.flydist_vec == 4)),
        x = 1:ngens;
        y1.mea = zeros(ngens,2);
        y1.fmin = zeros(ngens,2); y1.fmax = zeros(ngens,2);
        y1.min = zeros(ngens,2); y1.max = zeros(ngens,2);
        y2 = y1; y3 = y1;
        for igen=1:ngens,
            % Compute Statistics (mean, quartiles, min, max) for each genotype
            if params.winlos,
                ind1 = lunges{igen}.obj1.number >= lunges{igen}.obj2.number;
            else
                ind1 = ones(1,length(lunges{igen}.obj1.number));
            end
            ind2 = logical(1 - ind1);
            % Near (< 5 mm)
            tmp_w = [obj1{igen}.food_time.near(ind1) obj2{igen}.food_time.near(ind2)];
            tmp_l = [obj1{igen}.food_time.near(ind2) obj2{igen}.food_time.near(ind1)];
            tmp = [tmp_w' tmp_l'];
            y1.mea(igen,:) = mean(tmp,1);
            y1.fmin(igen,:) = - quantile(tmp,0.25,1) + y1.mea(igen,:);
            y1.fmax(igen,:) = quantile(tmp,0.75,1) - y1.mea(igen,:);
            y1.min(igen,:) = min(tmp); y1.max(igen,:) = max(tmp);
            % Mid-range (5 mm - 10 mm)
            tmp_w = [obj1{igen}.food_time.mid(ind1) obj2{igen}.food_time.mid(ind2)];
            tmp_l = [obj1{igen}.food_time.mid(ind2) obj2{igen}.food_time.mid(ind1)];
            tmp = [tmp_w' tmp_l'];
            y2.mea(igen,:) = mean(tmp,1);
            y2.fmin(igen,:) = - quantile(tmp,0.25,1) + y2.mea(igen,:);
            y2.fmax(igen,:) = quantile(tmp,0.75,1) - y2.mea(igen,:);
            y2.min(igen,:) = min(tmp); y3.max(igen,:) = max(tmp);
            % Far (? 10 mm)
            tmp_w = [obj1{igen}.food_time.far(ind1) obj2{igen}.food_time.far(ind2)];
            tmp_l = [obj1{igen}.food_time.far(ind2) obj2{igen}.food_time.far(ind1)];
            tmp = [tmp_w' tmp_l'];
            y3.mea(igen,:) = mean(tmp,1);
            y3.fmin(igen,:) = - quantile(tmp,0.25,1) + y3.mea(igen,:);
            y3.fmax(igen,:) = quantile(tmp,0.75,1) - y3.mea(igen,:);
            y3.min(igen,:) = min(tmp); y3.max(igen,:) = max(tmp);
        end
        
        % PLOT DATA AS STACK
        if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
        figure(FID); clf;
        set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
        % Two flies?
        if ~params.oneobj,
            d1 = -params.bardispl; d2 = -d1; w = params.twobarswidth;
            for i=1:ngens,
                patch([x(i)+d1-w x(i)+d1+w x(i)+d1+w x(i)+d1-w x(i)+d1-w],...
                    [0 0 y1.mea(i,1) y1.mea(i,1) 0],params.flycol.win); hold on;
                patch([x(i)+d1-w x(i)+d1+w x(i)+d1+w x(i)+d1-w x(i)+d1-w],...
                    [y1.mea(i,1) y1.mea(i,1) y1.mea(i,1)+y2.mea(i,1) y1.mea(i,1)+y2.mea(i,1) y1.mea(i,1)],params.flycol.los);
                patch([x(i)+d1-w x(i)+d1+w x(i)+d1+w x(i)+d1-w x(i)+d1-w],...
                    [y1.mea(i,1)+y2.mea(i,1) y1.mea(i,1)+y2.mea(i,1) y1.mea(i,1)+y2.mea(i,1)+y3.mea(i,1) ...
                    y1.mea(i,1)+y2.mea(i,1)+y3.mea(i,1) y1.mea(i,1)+y2.mea(i,1)],params.flycol.ext);
                patch([x(i)+d2-w x(i)+d2+w x(i)+d2+w x(i)+d2-w x(i)+d2-w],...
                    [0 0 y1.mea(i,2) y1.mea(i,2) 0],params.flycol.win);
                patch([x(i)+d2-w x(i)+d2+w x(i)+d2+w x(i)+d2-w x(i)+d2-w],...
                    [y1.mea(i,2) y1.mea(i,2) y1.mea(i,2)+y2.mea(i,2) y1.mea(i,2)+y2.mea(i,2) y1.mea(i,2)],params.flycol.los);
                patch([x(i)+d2-w x(i)+d2+w x(i)+d2+w x(i)+d2-w x(i)+d2-w],...
                    [y1.mea(i,2)+y2.mea(i,2) y1.mea(i,2)+y2.mea(i,2) y1.mea(i,2)+y2.mea(i,2)+y3.mea(i,2) ...
                    y1.mea(i,2)+y2.mea(i,2)+y3.mea(i,2) y1.mea(i,2)+y2.mea(i,2)],params.flycol.ext);
            end
            %         h = bar(x-.25,[y1.mea(:,1) y2.mea(:,1) y3.mea(:,1)],barw,'stack');
            %         colormap([params.flycol.win ; params.flycol.los ; params.flycol.ext]); hold on;
            %         h = bar(x-.25+barw,[y1.mea(:,2) y2.mea(:,2) y3.mea(:,2)],barw,'stack');
        else
            d1 = 0; w = params.barwidth;
            for i=1:ngens,
                patch([x(i)+d1-w x(i)+d1+w x(i)+d1+w x(i)+d1-w x(i)+d1-w],...
                    [0 0 y1.mea(i,1) y1.mea(i,1) 0],params.flycol.win); hold on;
                patch([x(i)+d1-w x(i)+d1+w x(i)+d1+w x(i)+d1-w x(i)+d1-w],...
                    [y1.mea(i,1) y1.mea(i,1) y1.mea(i,1)+y2.mea(i,1) y1.mea(i,1)+y2.mea(i,1) y1.mea(i,1)],params.flycol.los);
                patch([x(i)+d1-w x(i)+d1+w x(i)+d1+w x(i)+d1-w x(i)+d1-w],...
                    [y1.mea(i,1)+y2.mea(i,1) y1.mea(i,1)+y2.mea(i,1) y1.mea(i,1)+y2.mea(i,1)+y3.mea(i,1) ...
                    y1.mea(i,1)+y2.mea(i,1)+y3.mea(i,1) y1.mea(i,1)+y2.mea(i,1)],params.flycol.ext);
            end
            %         h = bar(x-.25,[y1.mea(:,1) y2.mea(:,1) y3.mea(:,1)],params.barwidth,'stack');
            %         colormap([params.flycol.win ; params.flycol.los ; params.flycol.ext]); hold on;
        end
        set(gca,'XTick',1:ngens); set_xtick_label(gen_ident,45,params.cross_to,params.axisfontsize);
        t = 0:300:params.max_frames; ytick_lab = cell(1,length(t));
        for i=1:length(t), ytick_lab{i} = num2str(t(i)/60); end
        set(gca,'YTick',t); set(gca,'YTickLabel',ytick_lab); set(gca,'FontSize',params.axisfontsize);
        ylabel('time on food [min] (median)');
        axis([0 ngens+1 0 params.max_frames]);
        h = legend('< 5 mm','5 - 10 mm','> 10 mm'); legend('boxoff');
        titl = 'Time Spent in Different Distances from Food';
        if ~params.oneobj,
            if params.winlos, titl = [titl ' L: winner, R: loser']; else titl = [titl ' L: fly 1, R: fly 2']; end;
        end
        set(FID,'Name',titl); mtit(titl,'FontSize',params.axisfontsize+2,'yoff',.04);
        if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
    end
end

% TRAVEL-DISTANCE (pre-copulation in male-female assays)
if params.plots.stat.traveldist,
    % Extract travel distance data
    data = win_los_data('movdists','max',obj1,obj2,lunges,copu.pre,dts,params);
    ngens = length(obj1);
    x = 1:ngens; y = []; y.mea = zeros(ngens,2); y.sem = y.mea;
    yt.mea = zeros(ngens,1); yt.sem = yt.mea; X = []; G = []; X2 = []; G2 = [];
    data1 = cell(ngens,1); data2 = data1;
    for igen=1:ngens,
        % Statistics for each genotype
        if params.oneobj, data{igen}.los(1,:) = 0; end
        y.mea(igen,:) = [mean(data{igen}.win) mean(data{igen}.los)]/1000;
        y.sem(igen,:) = [std(data{igen}.win)/sqrt(numel(data{igen}.win)) ...
            std(data{igen}.los)/sqrt(numel(data{igen}.los))]/1000;
        tmp = data{igen}.win + data{igen}.los;
        % For boxplots
        X = [X tmp/1000]; G = [G ones(1,length(tmp))*igen];
        if ~params.oneobj,
            X2 = [X2 data{igen}.win/1000 data{igen}.los/1000];
            G2 = [G2 ones(1,length(data{igen}.win))*(2*igen-1) ones(1,length(data{igen}.los))*2*igen];
        end
        % For barplots
        yt.mea(igen) = mean(tmp)/1000;
        yt.sem(igen) = std(tmp)/sqrt(numel(tmp))/1000;
        data1{igen} = data{igen}.win/1000; data2{igen} = data{igen}.los/1000;
    end

    % TEXT OUTPUT OF LUNGES AND TUSSLING PER METER AND PER MINUTE
    if params.bool_xls,
        % Travel Distance [m]
        textoutput([data1 data2],gen_ident,'distance [m]','Travel',params);
    end
    
    % In case of two flies take the average travel distance
    if ~params.oneobj,
        yt.mea = yt.mea / 2; yt.sem = yt.sem / 2; X = X / 2;
    end
    
    % PLOT DATA
    max_val = params.range.stat.occ.walk;
    titl = 'Total Distance Traveled per Fly';
    ytitl = 'traveled distance [m]';
    utest = utestpairs(gen_ident1);
    if ~params.oneobj,
        FID = plot_bar_msem(y,utest,X2,G2,max_val,ytitl,titl,FID,gen_ident,params); 
    end

    FID = plot_bar_msem(yt,utest,X,G,max_val,ytitl,titl,FID,gen_ident,params);
end
