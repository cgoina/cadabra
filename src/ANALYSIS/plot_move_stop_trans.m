%% Plot_move_stop_trans
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
%% Analyze and plot move/stop actions of flies

% MOVE / STOP ACTIONS
function FID = plot_move_stop_trans(obj1,obj2,lunges,gen_ident,params,dts,FID)

if params.plots.stat.movestop,
ngens = length(obj1);

% min/max velocities; velocities < thres_min are considered as stops
thres_min = 0.5; thres_max = 100;
gens = 1:ngens; % which genotypes to plot?
max_gap = 2; min_len = 5; % frames


minv = 0; maxv = 60; dv = 2; iv = minv:dv:maxv; col = jet(ngens); 
dv = 60; minv = dv/2; maxv = params.max_frames-dv/2; iv = minv:dv:maxv;
db = .5; minb = db/2; maxb = 30-db/2; ib = minb:db:maxb;
col = jet(ngens); 
mov.win.loc_med = cell(ngens,1); mov.win.loc_sem = cell(ngens,1); 
mov.win.bou_med = cell(ngens,1); mov.win.bou_sem = cell(ngens,1);
mov.win.t = cell(ngens,1); mov.win.l = cell(ngens,1); 
mov.win.lunf = cell(ngens,1); mov.win.n = cell(ngens,1);
mov.los = mov.win; stp = mov; nmovs = zeros(1,numel(gens));

for igen=1:numel(gens),
    ind = [0 find(obj1{gens(igen)}.t(2:end)-obj1{gens(igen)}.t(1:end-1)<0)]; 
    ind = [ind length(obj1{gens(igen)}.t)];
    nmovs(igen) = numel(ind)-1;
    
    mov0.bouts = []; mov0.win.h = zeros(numel(ind)-1,length(iv));
    mov0.win.h2d = zeros(numel(ind)-1,length(iv),length(ib));
    mov0.los = mov0.win; stp0 = mov0;
    for j=1:nmovs(igen),
        % TIME & VELOCITIES FOR EACH MOVIE
        if dts{gens(igen)}(j) < 1/15, fac = 2; else fac =1; end
        tim = obj1{gens(igen)}.t(ind(j)+1:ind(j+1));
        v1 = obj1{gens(igen)}.v(ind(j)+1:ind(j+1)); 
        v2 = obj2{gens(igen)}.v(ind(j)+1:ind(j+1));
        if (lunges{gens(igen)}.obj1.number(j) < lunges{gens(igen)}.obj2.number(j)) && params.winlos,
            tmp = v2; v2 = v1; v1 = tmp;
        end
        % EXTRACT MOVE / STOP INDICES FOR EACH FLY
        indmov1 = real((v1 >= thres_min) & (v1 < thres_max)); 
        indmov2 = real((v2 >= thres_min) & (v2 < thres_max));
        indstp1 = real((v1 < thres_min));
        indstp2 = real((v2 < thres_min));

        % MOVE / STOP BOUT DETECTION
        % Detect stop-bouts first
%         [stp0.bouts,indstp1,indstp2] = detect_bouts(tim,indstp1,indstp2,min_len*fac,max_gap/fac,stp0.bouts);        
%         indmov1 = indmov1*0; indmov1(indstp1) = 1; indmov1 = 1-indmov1;
%         indmov2 = indmov2*0; indmov2(indstp2) = 1; indmov2 = 1-indmov2;        
%         [mov0.bouts,indmov1,indmov2] = detect_bouts(tim,indmov1,indmov2,min_len*fac,max_gap/fac,mov0.bouts);
        % Detect move-bouts first
        [mov0.bouts,indmov1,indmov2] = detect_bouts(tim,indmov1,indmov2,min_len*fac,max_gap/fac,mov0.bouts);
        indstp1 = indstp1*0; indstp1(indmov1) = 1; indstp1 = 1-indstp1;
        indstp2 = indstp2*0; indstp2(indmov2) = 1; indstp2 = 1-indstp2;
        [stp0.bouts,indstp1,indstp2] = detect_bouts(tim,indstp1,indstp2,min_len*fac,max_gap/fac,stp0.bouts);        
        
        % MOVE / STOP STATISTICS
        % 2d histograms over bout time and length, histograms over bout time
        stmw = sum(mov0.bouts.obj1.number) - mov0.bouts.obj1.number(end) + 1; enmw = sum(mov0.bouts.obj1.number);
        stml = sum(mov0.bouts.obj2.number) - mov0.bouts.obj2.number(end) + 1; enml = sum(mov0.bouts.obj2.number);
        stsw = sum(stp0.bouts.obj1.number) - stp0.bouts.obj1.number(end) + 1; ensw = sum(stp0.bouts.obj1.number);
        stsl = sum(stp0.bouts.obj2.number) - stp0.bouts.obj2.number(end) + 1; ensl = sum(stp0.bouts.obj2.number);
        if numel(mov0.bouts.obj1.len),
            mov0.win.h2d(j,:,:) = hist3([mov0.bouts.obj1.t(stmw:enmw)' , mov0.bouts.obj1.len(stmw:enmw)' * dts{gens(igen)}(j).*fac],{iv ib});
            mov0.win.h(j,:) = hist(tim(indmov1), iv) * dts{gens(igen)}(j).*fac;
        end
        if numel(mov0.bouts.obj2.len),
            mov0.los.h2d(j,:,:) = hist3([mov0.bouts.obj2.t(stml:enml)' , mov0.bouts.obj2.len(stml:enml)' * dts{gens(igen)}(j).*fac],{iv ib});
            mov0.los.h(j,:) = hist(tim(indmov2), iv) * dts{gens(igen)}(j).*fac;
        end
        if numel(stp0.bouts.obj1.len),
            stp0.win.h2d(j,:,:) = hist3([stp0.bouts.obj1.t(stsw:ensw)' , stp0.bouts.obj1.len(stsw:ensw)' * dts{gens(igen)}(j).*fac],{iv ib});
            stp0.win.h(j,:) = hist(tim(indstp1), iv) * dts{gens(igen)}(j).*fac;
        end
        if numel(stp0.bouts.obj2.len),
            stp0.los.h2d(j,:,:) = hist3([stp0.bouts.obj2.t(stsl:ensl)' , stp0.bouts.obj2.len(stsl:ensl)' * dts{gens(igen)}(j).*fac],{iv ib});
            stp0.los.h(j,:) = hist(tim(indstp2), iv) * dts{gens(igen)}(j).*fac;
        end
    end
    
    % COLLECT MOVE / STOP BOUT INFORMATION FOR WINNER & LOSER
    % bout lengths
    indl = [0 cumsum(mov0.bouts.obj1.number)]; lenv = [];
    for j=1:numel(indl)-1, lenv = [lenv mov0.bouts.obj1.len(indl(j)+1:indl(j+1)) * dts{gens(igen)}(j).*fac]; end
    mov.win.l{gens(igen)} = lenv;
    % bout times
    mov.win.t{gens(igen)} = mov0.bouts.obj1.t;
    % unfiltered data
    mov.win.lunf{gens(igen)} = mov0.bouts.obj1.len_unfilt .* dts{gens(igen)}.*fac;
    % number of bouts 
    mov.win.n{gens(igen)} = mov0.bouts.obj1.number; 
    % mean histogram per genotype (bout number per time bin)
    mov.win.loc_med{gens(igen)} = mean(mov0.win.h,1); 
    % mean 2d histograms per genotype (bout time and length)
    mov.win.bou_med{gens(igen)} = reshape(mean(mov0.win.h2d,1),length(iv),length(ib));

    indl = [0 cumsum(stp0.bouts.obj1.number)]; lenv = [];
    for j=1:numel(indl)-1, lenv = [lenv stp0.bouts.obj1.len(indl(j)+1:indl(j+1)) * dts{gens(igen)}(j).*fac]; end
    stp.win.l{gens(igen)} = lenv;
    stp.win.t{gens(igen)} = stp0.bouts.obj1.t; 
    stp.win.lunf{gens(igen)} = stp0.bouts.obj1.len_unfilt .* dts{gens(igen)}.*fac;
    stp.win.n{gens(igen)} = stp0.bouts.obj1.number;
    stp.win.loc_med{gens(igen)} = mean(stp0.win.h,1); 
    stp.win.bou_med{gens(igen)} = reshape(mean(stp0.win.h2d,1),length(iv),length(ib));

    indl = [0 cumsum(mov0.bouts.obj2.number)]; lenv = [];
    for j=1:numel(indl)-1, lenv = [lenv mov0.bouts.obj2.len(indl(j)+1:indl(j+1)) * dts{gens(igen)}(j).*fac]; end
    mov.los.l{gens(igen)} = lenv;
    mov.los.t{gens(igen)} = mov0.bouts.obj2.t; 
    mov.los.lunf{gens(igen)} = mov0.bouts.obj2.len_unfilt .* dts{gens(igen)}.*fac;
    mov.los.n{gens(igen)} = mov0.bouts.obj2.number;
    mov.los.loc_med{gens(igen)} = mean(mov0.los.h,1); 
    mov.los.bou_med{gens(igen)} = reshape(mean(mov0.los.h2d,1),length(iv),length(ib));

    indl = [0 cumsum(stp0.bouts.obj2.number)]; lenv = [];
    for j=1:numel(indl)-1, lenv = [lenv stp0.bouts.obj2.len(indl(j)+1:indl(j+1)) * dts{gens(igen)}(j).*fac]; end
    stp.los.l{gens(igen)} = lenv;
    stp.los.t{gens(igen)} = stp0.bouts.obj2.t; 
    stp.los.lunf{gens(igen)} = stp0.bouts.obj2.len_unfilt .* dts{gens(igen)}.*fac;
    stp.los.n{gens(igen)} = stp0.bouts.obj2.number;
    stp.los.loc_med{gens(igen)} = mean(stp0.los.h,1); 
    stp.los.bou_med{gens(igen)} = reshape(mean(stp0.los.h2d,1),length(iv),length(ib));

    if numel(ind)>2,
        mov.win.loc_sem{gens(igen)} = std(mov0.win.h,1) / sqrt(size(mov0.win.h,1));
        mov.win.bou_sem{gens(igen)} = reshape(std(mov0.win.h2d,1),length(iv),length(ib)) / sqrt(size(mov0.win.h2d,1));
        stp.win.loc_sem{gens(igen)} = std(stp0.win.h,1) / sqrt(size(stp0.win.h,1));
        stp.win.bou_sem{gens(igen)} = reshape(std(stp0.win.h2d,1),length(iv),length(ib)) / sqrt(size(stp0.win.h2d,1));
        mov.los.loc_sem{gens(igen)} = std(mov0.los.h,1) / sqrt(size(mov0.los.h,1));
        mov.los.bou_sem{gens(igen)} = reshape(std(mov0.los.h2d,1),length(iv),length(ib)) / sqrt(size(mov0.los.h2d,1));
        stp.los.loc_sem{gens(igen)} = std(stp0.los.h,1) / sqrt(size(stp0.los.h,1));
        stp.los.bou_sem{gens(igen)} = reshape(std(stp0.los.h2d,1),length(iv),length(ib)) / sqrt(size(stp0.los.h2d,1));
    else
        mov.win.loc_sem{gens(igen)} = mov0.win.h * 0;
        mov.win.bou_sem{gens(igen)} = zeros(length(iv),length(ib));
        stp.win.loc_sem{gens(igen)} = stp0.win.h * 0;
        stp.win.bou_sem{gens(igen)} = zeros(length(iv),length(ib));
        mov.los.loc_sem{gens(igen)} = mov0.los.h * 0;
        mov.los.bou_sem{gens(igen)} = zeros(length(iv),length(ib));
        stp.los.loc_sem{gens(igen)} = stp0.los.h * 0;
        stp.los.bou_sem{gens(igen)} = zeros(length(iv),length(ib));
    end
end

% MEAN TIME SPENT IN MOVING / STOPING PER GENOTYPE
win_vec = zeros(ngens,2); los_vec = win_vec;
for igen=1:numel(gens),
    win_vec(gens(igen),:) = [mean(mov.win.lunf{gens(igen)}) ...
                             mean(stp.win.lunf{gens(igen)})];
    los_vec(gens(igen),:) = [mean(mov.los.lunf{gens(igen)}) ...
                             mean(stp.los.lunf{gens(igen)})];
end

% MEAN BOUT LENGTH AND FREQUENCY PER GENOTYPE AND TIME BIN
div1 = sum(mov.win.bou_med{gens(igen)},2); ind1 = find(div1 == 0);
if numel(ind1), div1(ind1) = 1; end
div2 = sum(stp.win.bou_med{gens(igen)},2); ind2 = find(div2 == 0);
if numel(ind2), div2(ind2) = 1; end
% Bout length winner
y1l_w = mov.win.bou_med{gens(igen)}*ib'./div1; if numel(ind1), y1l_w(ind1) = 0; end
y2l_w = stp.win.bou_med{gens(igen)}*ib'./div2; if numel(ind2), y2l_w(ind2) = 0; end
% Bout frequency winner
y1f_w = sum(mov.win.bou_med{gens(igen)},2);
y2f_w = sum(stp.win.bou_med{gens(igen)},2);
y1l_l = mov.los.bou_med{gens(igen)}*ib'./sum(mov.los.bou_med{gens(igen)},2);
y2l_l = stp.los.bou_med{gens(igen)}*ib'./sum(stp.los.bou_med{gens(igen)},2);
% Bout frequency loser
y1f_l = sum(mov.los.bou_med{gens(igen)},2);
y2f_l = sum(stp.los.bou_med{gens(igen)},2);
if params.autoscale, 
    max_bl = ceil(max([y1l_w + y2l_w ; y1l_l + y2l_l]));
    max_bf = ceil(max([y1f_w + y2f_w ; y1f_l + y2f_l])); 
    max_bl = max_bl+ceil(.1*abs(max_bl)); % extend the display range by 10%
    max_bf = max_bf+ceil(.1*abs(max_bf)); % extend the display range by 10%
else
    max_bl = 10; max_bf = 20;
end


% PLOTS
% Move ; stop color
barcolors = [.61 .92 .39 ; .97 .46 .34];
% Time axis length and increment - just for ploting the time axis
minv = 0; maxv = params.max_frames; dt = 300; % dt = 5 min.
t = minv:dt:maxv; ytick_lab = cell(1,length(t));
for i=1:length(t), ytick_lab{i} = num2str(t(i)/60); end

% PLOT TOTAL TIME SPENT MOVING (GREEN) / STOPING (RED) AS STACKED BARS
if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
figure(FID); clf;
set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
if ~params.oneobj,
    if ngens > 1,
        h = bar((1:ngens)-.38,win_vec,.25,'stack'); colormap(barcolors); hold on;
        h = bar((1:ngens)-.08,los_vec,.25,'stack');
    else
        h = bar([.8 1.2],[win_vec ; los_vec],.5,'stack'); 
        colormap(barcolors); hold on;
    end
else
    if ngens > 1,
        h = bar((1:ngens),win_vec,.25,'stack'); colormap(barcolors); hold on;
    else
        h = bar([win_vec ; nan(size(win_vec))],.25,'stack'); 
        colormap(barcolors); hold on;
    end
end
set(gca,'XTick',1:ngens);
if ngens > 1,
    set_xtick_label(gen_ident,45,params.cross_to,params.axisfontsize);
else
    set(gca,'XTickLabel',gen_ident);
end
set(gca,'YTick',t); set(gca,'YTickLabel', ytick_lab); set(gca,'FontSize',params.axisfontsize);
ylim([0 params.max_frames]); ylabel('MEAN period per fly pair [min]');
axis([0 ngens+1 0 params.max_frames]); h = legend('move','stop'); legend('boxoff');
titl = ['MOVE/STOP periods: v > ' num2str(thres_min) ' mm/s < ' num2str(thres_max) ...
        ' mm/s, min. length = ' num2str(min_len*dts{1}(1),'%6.1f') ' s'];
if ~params.oneobj, 
    if params.winlos, titl = [titl ', L: Winner, R: Loser']; else titl = [titl ', L: Fly 1, R: Fly 2']; end
end
set(FID,'Name',titl); mtit(titl,'FontSize',params.axisfontsize+2,'yoff',.0);
if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end


% PLOT BOUT ANALYSIS FOR WINNER OR FLY 1
appnd = '-append';
if length(gens)>=7, max_plts = 7; else max_plts = length(gens); end
for igen=1:length(gens),
    if mod(igen,max_plts) == 1 || length(gens) == 1,
        if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
        figure(FID); clf;
        set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
        icnt = 1;
    end
    
    % Plot time spent moving (green) / stoping (red) per one-minute time bin and genotype
    % dv = 60 seconds
    h = subplot(max_plts,3,3*(icnt-1)+1);
    errorbar(iv,stp.win.loc_med{gens(igen)}/dv*100,stp.win.loc_sem{gens(igen)}/dv*100,'o','Color',barcolors(2,:),'LineWidth',2); hold on;
    plot(iv,stp.win.loc_med{gens(igen)}/dv*100,'ob','MarkerEdgeColor',barcolors(2,:),'MarkerFaceColor','y','MarkerSize',6,'LineWidth',3,'Color',col(gens(igen),:)); hold on;
    errorbar(iv,mov.win.loc_med{gens(igen)}/dv*100,mov.win.loc_sem{gens(igen)}/dv*100,'o','Color',barcolors(1,:),'LineWidth',2); hold on;
    plot(iv,mov.win.loc_med{gens(igen)}/dv*100,'ob','MarkerEdgeColor',barcolors(1,:),'MarkerFaceColor','y','MarkerSize',6,'LineWidth',3,'Color',col(gens(igen),:)); hold on;
    set(h,'FontSize',params.axisfontsize);
    xlim([minv maxv]); 
    if ~mod(igen,max_plts) || igen == length(gens), xlabel('time [min]'); end
    if icnt == round(max_plts/2), ylabel('locomotion, % time (MEAN)'); end
    axis([minv maxv 0 100]); set(gca,'XTick',t); set(gca,'XTickLabel',ytick_lab);
    set(gca,'FontSize',params.axisfontsize);
    
    % Plot mean move / stop bout length per one-minute time bin and genotype
    h = subplot(max_plts,3,3*(icnt-1)+2);
    h1 = bar(iv,[y1l_w y2l_w],'stack'); colormap(barcolors); hold on;
    set(h,'FontSize',params.axisfontsize);
    xlim([minv maxv]); 
    if ~mod(igen,max_plts) || igen == length(gens), xlabel('time [min]'); end
    if icnt == round(max_plts/2), ylabel('MEAN bout length [s]'); end
    axis([minv maxv 0 max_bl]); title(gen_ident{gens(igen)});
    set(gca,'XTick',t); set(gca,'XTickLabel',ytick_lab);
    set(gca,'FontSize',params.axisfontsize);
    
    % Plot mean move / stop bout frequency per one-minute time bin and genotype
    h = subplot(max_plts,3,3*(icnt-1)+3);
    h1 = bar(iv,[y1f_w y2f_w],'stack'); colormap(barcolors); hold on;
    set(h,'FontSize',params.axisfontsize);
    xlim([minv maxv]); 
    if ~mod(igen,max_plts) || igen == length(gens), xlabel('time [min]'); end
    if icnt == round(max_plts/2), ylabel(['frequency of bouts / ' num2str(dv) ' s']); end
    axis([minv maxv 0 max_bf]); set(gca,'XTick',t); set(gca,'XTickLabel',ytick_lab);
    set(gca,'FontSize',params.axisfontsize);
    icnt = icnt + 1;
    if ~mod(igen,max_plts) || igen == length(gens),
        if ~params.oneobj, 
            if params.winlos, 
                titl = 'WINNER '; 
            else
                titl = 'FLY 1 ';
            end
        else
            titl = ''; 
        end
        titl = [titl 'MOVE/STOP bouts: v > ' num2str(thres_min) ' mm/s < ' num2str(thres_max) ' mm/s, min. length = ' num2str(min_len*dts{1}(1),'%6.1f') ' s'];
        set(FID,'Name',titl); mtit(titl,'FontSize',params.axisfontsize+2,'yoff',.03);
        if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
    end
end

% PLOT BOUT ANALYSIS FOR LOSER OR FLY 2
if ~params.oneobj,
    for igen=1:length(gens),
        if mod(igen,max_plts) == 1 || length(gens) == 1,
            if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
            figure(FID); clf;
            set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
            icnt = 1;
        end
        
        % Plot time spent moving (green) / stoping (red) per one-minute time bin and genotype
        % dv = 60 seconds
        h = subplot(max_plts,3,3*(icnt-1)+1);
        errorbar(iv,stp.los.loc_med{gens(igen)}/dv*100,stp.los.loc_sem{gens(igen)}/dv*100,'o','Color',barcolors(2,:),'LineWidth',2); hold on;
        plot(iv,stp.los.loc_med{gens(igen)}/dv*100,'ob','MarkerEdgeColor',barcolors(2,:),'MarkerFaceColor','y','MarkerSize',6,'LineWidth',3,'Color',col(gens(igen),:)); hold on;
        errorbar(iv,mov.los.loc_med{gens(igen)}/dv*100,mov.los.loc_sem{gens(igen)}/dv*100,'o','Color',barcolors(1,:),'LineWidth',2); hold on;
        plot(iv,mov.los.loc_med{gens(igen)}/dv*100,'ob','MarkerEdgeColor',barcolors(1,:),'MarkerFaceColor','y','MarkerSize',6,'LineWidth',3,'Color',col(gens(igen),:)); hold on;
        set(h,'FontSize',params.axisfontsize);
        xlim([minv maxv]);
        if ~mod(igen,max_plts) || igen == length(gens), xlabel('time [min]'); end
        if icnt == round(max_plts/2), ylabel('locomotion, % time (MEAN)'); end
        axis([minv maxv 0 100]); set(gca,'XTick',t); set(gca,'XTickLabel',ytick_lab);
        set(gca,'FontSize',params.axisfontsize);
        
        % Plot mean move / stop bout length per one-minute time bin and genotype
        h = subplot(max_plts,3,3*(icnt-1)+2);
        h1 = bar(iv,[y1l_l y2l_l],'stack'); colormap(barcolors); hold on;
        set(h,'FontSize',params.axisfontsize);
        xlim([minv maxv]);
        if ~mod(igen,max_plts) || igen == length(gens), xlabel('time [min]'); end
        if icnt == round(max_plts/2), ylabel('MEAN bout length [s]'); end
        axis([minv maxv 0 max_bl]); title(gen_ident{gens(igen)});
        set(gca,'XTick',t); set(gca,'XTickLabel',ytick_lab);
        set(gca,'FontSize',params.axisfontsize);

        % Plot mean move / stop bout frequency per one-minute time bin and genotype
        h = subplot(max_plts,3,3*(icnt-1)+3);
        h1 = bar(iv,[y1f_l y2f_l],'stack'); colormap(barcolors); hold on;
        set(h,'FontSize',params.axisfontsize);
        xlim([minv maxv]);
        if ~mod(igen,max_plts) || igen == length(gens), xlabel('time [min]'); end
        if icnt == round(max_plts/2), ylabel(['frequency of bouts / ' num2str(dv) ' s']); end
        axis([minv maxv 0 max_bf]); set(gca,'XTick',t); set(gca,'XTickLabel',ytick_lab);
        set(gca,'FontSize',params.axisfontsize);
        icnt = icnt + 1;
        if ~mod(igen,max_plts) || igen == length(gens),
            if ~params.oneobj,
                if params.winlos,
                    titl = 'LOSER ';
                else
                    titl = 'FLY 2 ';
                end
            else
                titl = '';
            end
            titl = [titl 'MOVE/STOP bouts: v > ' num2str(thres_min) ' mm/s < ' num2str(thres_max) ' mm/s, min. length = ' num2str(min_len*dts{1}(1),'%6.1f') ' s'];
            set(FID,'Name',titl); mtit(titl,'FontSize',params.axisfontsize+2,'yoff',.03);
            if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
        end
    end
end

end