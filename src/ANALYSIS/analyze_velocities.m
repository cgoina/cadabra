%% Analyze_velocities
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
%% Analyze fly velocities and compute various statistics, i.e. velocities
%% at different proximities between flies

function FID = analyze_velocities(obj1,obj2,lunges,proxi,dts,nframes,gen_ident,params,FID)

if params.tim,
    max_vel = params.range.stat.tim.vel;
else
    max_vel = params.range.stat.occ.vel;
end
ngens = length(lunges);
% FLY VELOCITIES UNDER DIFFERENT PROXIMITIES
if params.plots.stat.velocity,
    if ~params.oneobj,
        y1.mea = zeros(ngens,2); y1.min = y1.mea; y1.max = y1.mea;
        y1.sem = y1.mea;
        y2 = y1; y3 = y1;
        X1 = []; G1 = []; X2 = []; G2 = []; X3 = []; G3 = [];
        for igen=1:ngens,
            % Distinguish between winner and loser
            ind1 = lunges{igen}.obj1.number >= lunges{igen}.obj2.number;
            ind2 = logical(1 - ind1); ind1 = find(ind1); ind2 = find(ind2);

            % COMPUTE STATISTICS FOR FLY VELOCITIES AT RELATIVE DISTANCES 
            % Near (< 5 mm)
            ind = [0 cumsum(proxi.indi{igen}.near)]; indw = []; indl = [];
            for j=1:numel(ind1), indw = [indw ind(ind1(j))+1:ind(ind1(j)+1)]; end
            for j=1:numel(ind2), indl = [indl ind(ind2(j))+1:ind(ind2(j)+1)]; end
            tmp_w = [obj1{igen}.proxi_vel.near(indw) obj2{igen}.proxi_vel.near(indl)];
            tmp_l = [obj1{igen}.proxi_vel.near(indl) obj2{igen}.proxi_vel.near(indw)];
            tmp = [tmp_w' tmp_l'];
            y1.mea(igen,:) = mean(tmp,1);
            y1.sem(igen,:) = std(tmp,[],1)/sqrt(numel(lunges{igen}.mov));
            y1.min(igen,:) = min(tmp); y1.max(igen,:) = max(tmp);
            if params.oneobj,
                X1 = [X1 tmp]; G1 = [G1 ones(1,length(tmp))*igen];
            else
                X1 = [X1 tmp_w tmp_l];
                G1 = [G1 ones(1,length(tmp_w))*(2*igen-1) ones(1,length(tmp_l))*2*igen];
            end
            % Mid-range (5 mm - 10 mm)
            ind = [0 cumsum(proxi.indi{igen}.mid)]; indw = []; indl = [];
            for j=1:numel(ind1), indw = [indw ind(ind1(j))+1:ind(ind1(j)+1)]; end
            for j=1:numel(ind2), indl = [indl ind(ind2(j))+1:ind(ind2(j)+1)]; end
            tmp_w = [obj1{igen}.proxi_vel.mid(indw) obj2{igen}.proxi_vel.mid(indl)];
            tmp_l = [obj1{igen}.proxi_vel.mid(indl) obj2{igen}.proxi_vel.mid(indw)];
            tmp = [tmp_w' tmp_l'];
            y2.mea(igen,:) = mean(tmp,1);
            y2.sem(igen,:) = std(tmp,[],1)/sqrt(numel(lunges{igen}.mov));
            y2.min(igen,:) = min(tmp); y2.max(igen,:) = max(tmp);
            if params.oneobj,
                X2 = [X2 tmp]; G2 = [G2 ones(1,length(tmp))*igen];
            else
                X2 = [X2 tmp_w tmp_l];
                G2 = [G2 ones(1,length(tmp_w))*(2*igen-1) ones(1,length(tmp_l))*2*igen];
            end
            % Far (? 10 mm)
            ind = [0 cumsum(proxi.indi{igen}.far)]; indw = []; indl = [];
            for j=1:numel(ind1), indw = [indw ind(ind1(j))+1:ind(ind1(j)+1)]; end
            for j=1:numel(ind2), indl = [indl ind(ind2(j))+1:ind(ind2(j)+1)]; end
            tmp_w = [obj1{igen}.proxi_vel.far(indw) obj2{igen}.proxi_vel.far(indl)];
            tmp_l = [obj1{igen}.proxi_vel.far(indl) obj2{igen}.proxi_vel.far(indw)];
            tmp = [tmp_w' tmp_l'];
            y3.mea(igen,:) = mean(tmp,1);
            y3.sem(igen,:) = std(tmp,[],1)/sqrt(numel(lunges{igen}.mov));
            y3.min(igen,:) = min(tmp); y3.max(igen,:) = max(tmp);
            if params.oneobj,
                X3 = [X3 tmp]; G3 = [G3 ones(1,length(tmp))*igen];
            else
                X3 = [X3 tmp_w tmp_l];
                G3 = [G3 ones(1,length(tmp_w))*(2*igen-1) ones(1,length(tmp_l))*2*igen];
            end
        end
        
        % PLOTS
        utest = utestpairs(gen_ident);
        if numel(find(params.plots.stat.vel_vec == 1)),
            ytitl = 'fly velocity [mm/s]';
            titl = 'Fly Proximity < 5 mm';
            FID = plot_bar_msem(y1,utest,X1,G1,max_vel,ytitl,titl,FID,gen_ident,params);

            titl = 'Fly Proximity 5 - 10 mm';
            FID = plot_bar_msem(y2,utest,X2,G2,max_vel,ytitl,titl,FID,gen_ident,params);

            titl = 'Fly Proximity > 10 mm';
            FID = plot_bar_msem(y3,utest,X3,G3,max_vel,ytitl,titl,FID,gen_ident,params);
        end
    end

    % FLY VELOCITIES UNDER DIFFERENT PROXIMITIES
    % PLOT OF ALL 3 RANGES
    if ~params.oneobj && numel(find(params.plots.stat.vel_vec == 1)),
        x = 1:ngens; y.mea = zeros(ngens,3); y.sem = y.mea;
        y.min = y.mea; y.max = y.mea;
        X = []; G = []; icnt = 0;
        for igen=1:ngens,
            % Near (< 5 mm)
            icnt = icnt + 1;
            data = [obj1{igen}.proxi_vel.near obj2{igen}.proxi_vel.near];
            y.mea(igen,1) = mean(data);
            y.sem(igen,1) = std(data)/sqrt(numel(lunges{igen}.mov));
            y.min(igen,1) = min(data); y.max(igen,1) = max(data);
            X = [X data]; G = [G ones(1,numel(data))*icnt];
            % Mid-range (5 mm - 10 mm)
            icnt = icnt + 1;
            data = [obj1{igen}.proxi_vel.mid obj2{igen}.proxi_vel.mid];
            y.mea(igen,2) = mean(data);
            y.sem(igen,2) = std(data)/sqrt(numel(lunges{igen}.mov));
            y.min(igen,2) = min(data); y.max(igen,2) = max(data);
            X = [X data]; G = [G ones(1,numel(data))*icnt];
            % Far (> 10 mm)
            icnt = icnt + 1;
            data = [obj1{igen}.proxi_vel.far obj2{igen}.proxi_vel.far];
            y.mea(igen,3) = mean(data);
            y.sem(igen,3) = std(data)/sqrt(numel(lunges{igen}.mov));
            y.min(igen,3) = min(data); y.max(igen,3) = max(data);
            X = [X data]; G = [G ones(1,numel(data))*icnt];
        end
        
        % PLOT
        utest = utestpairs(gen_ident);
        ytitl = 'fly velocity [mm/s]';
        titl = 'Fly Velocity under Different Proximities';
        FID = plot_bar_msem(y,utest,X,G,max_vel,ytitl,titl,FID,gen_ident,params);
    end
    
    % DURATON OF VELOCITIES
    if numel(find(params.plots.stat.vel_vec == 2)),
        if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
        figure(FID); clf;
        set(FID,'PaperOrientation','portrait','PaperPositionMode','manual','PaperPosition',[0 0 8.5 11]);
        minv = 0; maxv = 40; dv = 2; iv = minv:dv:maxv; col = jet(ngens); vsem = cell(ngens,1);
        h1 = subplot(2,1,1);
        for igen=1:ngens,
            % Determine time difference between two frames
            dt = sum(dts{igen}.*nframes{igen})/sum(nframes{igen});
            % Adjust for the case where every 2nd frame was removed
            % from obj1{}.v to save memory
            % Removed - we cannot assume velocity instances for frames
            % we did not consider
%             if dt < 1/15, fac = 2; else fac = 1; end
            % DETERMINE START/END OF EACH MOVIE
            ind = [0 find(obj1{igen}.t(2:end)-obj1{igen}.t(1:end-1)<0) length(obj1{igen}.t)];
            h = zeros(numel(ind),length(iv)); v1 = obj1{igen}.v; v2 = obj2{igen}.v; %sqrt(obj1{igen}.vx.^2 + obj1{igen}.vy.^2);
            % COMPUTE ONE VELOCITY HISTOGRAM PER MOVIE and multiply occurences 
            % with dt to get the total time spent in each velocity bin
            for j=1:numel(ind)-1,
                h(j,:) = hist([v1(ind(j)+1:ind(j+1)) ...
                               v2(ind(j)+1:ind(j+1))],iv) * dt; %* fac;
                %         h(j,:) = hist(v1(ind(j)+1:ind(j+1)),iv) * dt;
            end
            % STATISTICS
            vmed = mean(h); vsem{igen} = std(h) / sqrt(length(h));
            vmed(end) = 0; vsem{igen}(end) = 0;
            % PLOT on a semi-logarithmic scale for the time spent
            semilogy(iv,vmed,'Color',col(igen,:),'LineWidth',4); hold on;
            %     errorbar(iv,v.med,v.fmin,v.fmax,'ok','LineWidth',2); hold on;
        end
        set(h1,'FontSize',params.axisfontsize);
        xlim([minv maxv]); xlabel('velocity [mm/s]'); ylabel('MEAN duration [s]');
        axis([minv maxv 0.5 1000]); %legend(gen_ident);
        h = subplot(2,1,2);
        for igen=1:ngens,
            semilogy(iv,vsem{igen},'Color',col(igen,:),'LineWidth',4); hold on;
        end
        set(h,'FontSize',params.axisfontsize);
        xlim([minv maxv]); xlabel('velocity [mm/s]'); ylabel('SEM duration [s]');
        axis([minv maxv 0.5 1000]);
        titl = 'Total Duration of Velocities'; h = legend(gen_ident); legend(h,'boxoff');
        set(FID,'Name',titl); mtit(titl,'FontSize',params.axisfontsize+2,'yoff',.0);
        if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
    end
    
    if ~params.oneobj,
        if numel(find(params.plots.stat.vel_vec == 3)),
            % FLY VELOCITY - APPROACH
            x = 1:ngens; y1.mea = zeros(ngens,2);
            y1.min = y1.mea; y1.max = y1.mea; y1.sem = y1.mea; y2 = y1; y3 = y1;
            X1 = []; G1 = []; X2 = []; G2 = []; X3 = []; G3 = [];
            for igen=1:ngens,
                % Distinguish between fly 1/2 or winner/loser
                if ~params.winlos,
                    ind1 = ones(1,length(lunges{igen}.obj1.number));
                else
                    ind1 = lunges{igen}.obj1.number >= lunges{igen}.obj2.number;
                end
                ind2 = logical(1 - ind1); ind1 = find(ind1); ind2 = find(ind2);

                % FLY VELOCITIES UNDER DIFFERENT PROXIMITIES AND FOR 
                % HEAD DIRECTION TOWARDS THE OPPONENT (< 45 DEGREES)
                % Near (< 5 mm)
                ind = [0 cumsum(proxi.indi{igen}.near)]; indw = []; indl = [];
                for j=1:numel(ind1), indw = [indw ind(ind1(j))+1:ind(ind1(j)+1)]; end
                for j=1:numel(ind2), indl = [indl ind(ind2(j))+1:ind(ind2(j)+1)]; end
                data1 = []; data2 = [];
                if numel(indw),
                    data1 = [data1 obj1{igen}.proxi_vel.near(abs(obj1{igen}.proxi_mvdirdiff.near(indw)) < 45)];
                    data2 = [data2 obj2{igen}.proxi_vel.near(abs(obj2{igen}.proxi_mvdirdiff.near(indw)) < 45)];
                end
                if numel(indl),
                    data1 = [data1 obj2{igen}.proxi_vel.near(abs(obj2{igen}.proxi_mvdirdiff.near(indl)) < 45)];
                    data2 = [data2 obj1{igen}.proxi_vel.near(abs(obj1{igen}.proxi_mvdirdiff.near(indl)) < 45)];
                end
                % Prepare data for boxplots
                if params.oneobj,
                    X1 = [X1 data1]; G1 = [G1 ones(1,length(data1))*igen];
                else
                    X1 = [X1 data1 data2];
                    G1 = [G1 ones(1,length(data1))*(2*igen-1) ones(1,length(data2))*2*igen];
                end
                % Statistics for barplots
                y1.mea(igen,1) = mean(data1);
                y1.sem(igen,1) = std(data1)/sqrt(numel(lunges{igen}.mov));
                y1.min(igen,1) = min(data1);
                y1.max(igen,1) = max(data1);
                y1.mea(igen,2) = mean(data2);
                y1.sem(igen,2) = std(data2)/sqrt(numel(lunges{igen}.mov));
                y1.min(igen,2) = min(data2);
                y1.max(igen,2) = max(data2);

                % Mid-range (5 mm - 10 mm)
                ind = [0 cumsum(proxi.indi{igen}.mid)]; indw = []; indl = [];
                for j=1:numel(ind1), indw = [indw ind(ind1(j))+1:ind(ind1(j)+1)]; end
                for j=1:numel(ind2), indl = [indl ind(ind2(j))+1:ind(ind2(j)+1)]; end
                data1 = []; data2 = [];
                if numel(indw),
                    data1 = [data1 obj1{igen}.proxi_vel.mid(abs(obj1{igen}.proxi_mvdirdiff.mid(indw)) < 45)];
                    data2 = [data2 obj2{igen}.proxi_vel.mid(abs(obj2{igen}.proxi_mvdirdiff.mid(indw)) < 45)];
                end
                if numel(indl),
                    data1 = [data1 obj2{igen}.proxi_vel.mid(abs(obj2{igen}.proxi_mvdirdiff.mid(indl)) < 45)];
                    data2 = [data2 obj1{igen}.proxi_vel.mid(abs(obj1{igen}.proxi_mvdirdiff.mid(indl)) < 45)];
                end
                % Prepare data for boxplots
                if params.oneobj,
                    X2 = [X2 data1]; G2 = [G2 ones(1,length(data1))*igen];
                else
                    X2 = [X2 data1 data2];
                    G2 = [G2 ones(1,length(data1))*(2*igen-1) ones(1,length(data2))*2*igen];
                end
                % Statistics for barplots
                y2.mea(igen,1) = mean(data1);
                y2.sem(igen,1) = std(data1)/sqrt(numel(lunges{igen}.mov));
                y2.min(igen,1) = min(data1);
                y2.max(igen,1) = max(data1);
                y2.mea(igen,2) = mean(data2);
                y2.sem(igen,2) = std(data2)/sqrt(numel(lunges{igen}.mov));
                y2.min(igen,2) = min(data2);
                y2.max(igen,2) = max(data2);

                % Far (> 10 mm)
                ind = [0 cumsum(proxi.indi{igen}.far)]; indw = []; indl = [];
                for j=1:numel(ind1), indw = [indw ind(ind1(j))+1:ind(ind1(j)+1)]; end
                for j=1:numel(ind2), indl = [indl ind(ind2(j))+1:ind(ind2(j)+1)]; end
                data1 = []; data2 = [];
                if numel(indw),
                    data1 = [data1 obj1{igen}.proxi_vel.far(abs(obj1{igen}.proxi_mvdirdiff.far(indw)) < 45)];
                    data2 = [data2 obj2{igen}.proxi_vel.far(abs(obj2{igen}.proxi_mvdirdiff.far(indw)) < 45)];
                end
                if numel(indl),
                    data1 = [data1 obj2{igen}.proxi_vel.far(abs(obj2{igen}.proxi_mvdirdiff.far(indl)) < 45)];
                    data2 = [data2 obj1{igen}.proxi_vel.far(abs(obj1{igen}.proxi_mvdirdiff.far(indl)) < 45)];
                end
                % Prepare data for boxplots
                if params.oneobj,
                    X3 = [X3 data1]; G3 = [G3 ones(1,length(data1))*igen];
                else
                    X3 = [X3 data1 data2];
                    G3 = [G3 ones(1,length(data1))*(2*igen-1) ones(1,length(data2))*2*igen];
                end
                % Statistics for barplots
                y3.mea(igen,1) = mean(data1);
                y3.sem(igen,1) = std(data1)/sqrt(numel(lunges{igen}.mov));
                y3.min(igen,1) = min(data1);
                y3.max(igen,1) = max(data1);
                y3.mea(igen,2) = mean(data2);
                y3.sem(igen,2) = std(data2)/sqrt(numel(lunges{igen}.mov));
                y3.min(igen,2) = min(data2);
                y3.max(igen,2) = max(data2);
            end

            % PLOTS
            utest = utestpairs(gen_ident);
            ytitl = 'MEAN fly velocity [mm/s]';
            titl = 'Approach, Range < 5 mm';
            FID = plot_bar_msem(y1,utest,X1,G1,max_vel,ytitl,titl,FID,gen_ident,params);

            titl = 'Approach, Range 5 - 10 mm';
            FID = plot_bar_msem(y2,utest,X1,G1,max_vel,ytitl,titl,FID,gen_ident,params);

            titl = 'Approach, Range > 10 mm';
            FID = plot_bar_msem(y3,utest,X1,G1,max_vel,ytitl,titl,FID,gen_ident,params);
        end

        % FLY VELOCITY - ESCAPE
        if numel(find(params.plots.stat.vel_vec == 4)),
            x = 1:ngens; y1.mea = zeros(ngens,2);
            y1.fmin = y1.mea; y1.fmax = y1.mea; y2 = y1; y3 = y1;
            X1 = []; G1 = []; X2 = []; G2 = []; X3 = []; G3 = [];
            for igen=1:ngens,
                if ~params.winlos,
                    ind1 = ones(1,length(lunges{igen}.obj1.number));
                else
                    ind1 = lunges{igen}.obj1.number >= lunges{igen}.obj2.number;
                end
                ind2 = logical(1 - ind1); ind1 = find(ind1); ind2 = find(ind2);

                % FLY VELOCITIES UNDER DIFFERENT PROXIMITIES AND FOR 
                % HEAD DIRECTION AWAY FROM THE OPPONENT (> 135 DEGREES)
                % Near (< 5 mm)
                ind = [0 cumsum(proxi.indi{igen}.near)]; indw = []; indl = [];
                for j=1:numel(ind1), indw = [indw ind(ind1(j))+1:ind(ind1(j)+1)]; end
                for j=1:numel(ind2), indl = [indl ind(ind2(j))+1:ind(ind2(j)+1)]; end
                data1 = []; data2 = [];
                if numel(indw),
                    data1 = [data1 obj1{igen}.proxi_vel.near(abs(obj1{igen}.proxi_mvdirdiff.near(indw)) > 135)];
                    data2 = [data2 obj2{igen}.proxi_vel.near(abs(obj2{igen}.proxi_mvdirdiff.near(indw)) > 135)];
                end
                if numel(indl),
                    data1 = [data1 obj2{igen}.proxi_vel.near(abs(obj2{igen}.proxi_mvdirdiff.near(indl)) > 135)];
                    data2 = [data2 obj1{igen}.proxi_vel.near(abs(obj1{igen}.proxi_mvdirdiff.near(indl)) > 135)];
                end
                % Prepare data for boxplots
                if params.oneobj,
                    X1 = [X1 data1]; G1 = [G1 ones(1,length(data1))*igen];
                else
                    X1 = [X1 data1 data2];
                    G1 = [G1 ones(1,length(data1))*(2*igen-1) ones(1,length(data2))*2*igen];
                end
                % Statistics for barplots
                y1.mea(igen,1) = mean(data1);
                y1.sem(igen,1) = std(data1)/sqrt(numel(lunges{igen}.mov));
                y1.min(igen,1) = min(data1);
                y1.max(igen,1) = max(data1);
                y1.mea(igen,2) = mean(data2);
                y1.sem(igen,2) = std(data2)/sqrt(numel(lunges{igen}.mov));
                y1.min(igen,2) = min(data2);
                y1.max(igen,2) = max(data2);

                % Mid-range (5 mm - 10 mm)
                ind = [0 cumsum(proxi.indi{igen}.mid)]; indw = []; indl = [];
                for j=1:numel(ind1), indw = [indw ind(ind1(j))+1:ind(ind1(j)+1)]; end
                for j=1:numel(ind2), indl = [indl ind(ind2(j))+1:ind(ind2(j)+1)]; end
                data1 = []; data2 = [];
                if numel(indw),
                    data1 = [data1 obj1{igen}.proxi_vel.mid(abs(obj1{igen}.proxi_mvdirdiff.mid(indw)) > 135)];
                    data2 = [data2 obj2{igen}.proxi_vel.mid(abs(obj2{igen}.proxi_mvdirdiff.mid(indw)) > 135)];
                end
                if numel(indl),
                    data1 = [data1 obj2{igen}.proxi_vel.mid(abs(obj2{igen}.proxi_mvdirdiff.mid(indl)) > 135)];
                    data2 = [data2 obj1{igen}.proxi_vel.mid(abs(obj1{igen}.proxi_mvdirdiff.mid(indl)) > 135)];
                end
                % Prepare data for boxplots
                if params.oneobj,
                    X2 = [X2 data1]; G2 = [G2 ones(1,length(data1))*igen];
                else
                    X2 = [X2 data1 data2];
                    G2 = [G2 ones(1,length(data1))*(2*igen-1) ones(1,length(data2))*2*igen];
                end
                % Statistics for barplots                
                y2.mea(igen,1) = mean(data1);
                y2.sem(igen,1) = std(data1)/sqrt(numel(lunges{igen}.mov));
                y2.min(igen,1) = min(data1);
                y2.max(igen,1) = max(data1);
                y2.mea(igen,2) = mean(data2);
                y2.sem(igen,2) = std(data2)/sqrt(numel(lunges{igen}.mov));
                y2.min(igen,2) = min(data2);
                y2.max(igen,2) = max(data2);

                % Far (> 10 mm)
                ind = [0 cumsum(proxi.indi{igen}.far)]; indw = []; indl = [];
                for j=1:numel(ind1), indw = [indw ind(ind1(j))+1:ind(ind1(j)+1)]; end
                for j=1:numel(ind2), indl = [indl ind(ind2(j))+1:ind(ind2(j)+1)]; end
                data1 = []; data2 = [];
                if numel(indw),
                    data1 = [data1 obj1{igen}.proxi_vel.far(abs(obj1{igen}.proxi_mvdirdiff.far(indw)) > 135)];
                    data2 = [data2 obj2{igen}.proxi_vel.far(abs(obj2{igen}.proxi_mvdirdiff.far(indw)) > 135)];
                end
                if numel(indl),
                    data1 = [data1 obj2{igen}.proxi_vel.far(abs(obj2{igen}.proxi_mvdirdiff.far(indl)) > 135)];
                    data2 = [data2 obj1{igen}.proxi_vel.far(abs(obj1{igen}.proxi_mvdirdiff.far(indl)) > 135)];
                end
                % Prepare data for boxplots
                if params.oneobj,
                    X3 = [X3 data1]; G3 = [G3 ones(1,length(data1))*igen];
                else
                    X3 = [X3 data1 data2];
                    G3 = [G3 ones(1,length(data1))*(2*igen-1) ones(1,length(data2))*2*igen];
                end
                % Statistics for barplots                
                y3.mea(igen,1) = mean(data1);
                y3.sem(igen,1) = std(data1)/sqrt(numel(lunges{igen}.mov));
                y3.min(igen,1) = min(data1);
                y3.max(igen,1) = max(data1);
                y3.mea(igen,2) = mean(data2);
                y3.sem(igen,2) = std(data2)/sqrt(numel(lunges{igen}.mov));
                y3.min(igen,2) = min(data2);
                y3.max(igen,2) = max(data2);
            end

            % PLOTS
            utest = utestpairs(gen_ident);
            ytitl = 'MEAN fly velocity [mm/s]';
            titl = 'Escape, Range < 5 mm';
            FID = plot_bar_msem(y1,utest,X1,G1,max_vel,ytitl,titl,FID,gen_ident,params);

            titl = 'Escape, Range 5 - 10 mm';
            FID = plot_bar_msem(y2,utest,X1,G1,max_vel,ytitl,titl,FID,gen_ident,params);

            titl = 'Escape, Range > 10 mm';
            FID = plot_bar_msem(y3,utest,X1,G1,max_vel,ytitl,titl,FID,gen_ident,params);
        end
    end
end

