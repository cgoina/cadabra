%% Compute time series information on aggressive,
%% chasing, and courtship actions

function FID = beh_occ_ts(wings,lunges,tussls,chases,courts,gen_ident,ngens,nmovs,dts,params,FID)

dtb = 30; % [s]

icnt = 1; nplt = 2;
k2 = floor(sqrt(nplt)); k1 = ceil(nplt/k2);
ind = strfind(params.PSFileN,params.slash); path = params.PSFileN(1:ind(end));
[s,mess,messid] = mkdir([path 'tables' params.slash]); %#ok<NASGU>
for igen=1:ngens,
    if ~mod(icnt-1,nplt),
        if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
        icnt = 1;
        figure(FID); clf;
        set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
    end
    
    fname = [path 'tables' params.slash 'behavior-timeseries_' cell2mat(gen_ident(igen)) '.txt'];
    FIDT = fopen(fname,'w');
    fprintf(FIDT,['Time [s] spent in aggression, chasing, courtship per ' ...
            num2str(dtb,'%4g') 's bin - ' cell2mat(gen_ident(igen)) '\n\n']);
    fprintf(FIDT,'%11s','Time point');
    for imov=1:nmovs(igen),
        fprintf(FIDT,' %6s %6s %6s ','aggr','chase','court');
    end
            
    dt = median(dts{igen});
    binsz = ceil(dtb / dt);    
    nfrms = ceil(params.max_frames/min(dts{igen}));
    aggr_ts = zeros(nfrms,nmovs(igen));
    chas_ts = aggr_ts; cour_ts = aggr_ts;
    cl = [0 cumsum(lunges{igen}.number)];
    ct = [0 cumsum(tussls{igen}.number)];
    cwthr1 = [0 cumsum(wings.threat{igen}.obj1.number)]; cwthr2 = [0 cumsum(wings.threat{igen}.obj2.number)];
    cch1 = [0 cumsum(chases{igen}.obj1.number)]; cch2 = [0 cumsum(chases{igen}.obj2.number)];
    cwextl1 = [0 cumsum(wings.ext.l{igen}.obj1.number)]; cwextl2 = [0 cumsum(wings.ext.l{igen}.obj2.number)];
    cwextr1 = [0 cumsum(wings.ext.r{igen}.obj1.number)]; cwextr2 = [0 cumsum(wings.ext.r{igen}.obj2.number)];
    ccrt1 = [0 cumsum(courts{igen}.obj1.number)]; ccrt2 = [0 cumsum(courts{igen}.obj2.number)];
    
    for imov=1:nmovs(igen),
        l = [0 wings.threat{igen}.obj1.len];
        st = sum(l(1:cwthr1(imov)+1))+1; en = sum(l(1:cwthr1(imov+1)+1));
        aggr_ts(wings.threat{igen}.obj1.lab(st:en)',imov) = aggr_ts(wings.threat{igen}.obj1.lab(st:en)',imov) + 1;
        l = [0 wings.threat{igen}.obj2.len];
        st = sum(l(1:cwthr2(imov)+1))+1; en = sum(l(1:cwthr2(imov+1)+1));
        aggr_ts(wings.threat{igen}.obj2.lab(st:en)',imov) = aggr_ts(wings.threat{igen}.obj2.lab(st:en)',imov) + 1;
        st = cl(imov)+1; en = cl(imov+1);
        aggr_ts(lunges{igen}.lab(st:en)',imov) = aggr_ts(lunges{igen}.lab(st:en)',imov) + 1;
        l = [0 tussls{igen}.len];
        st = sum(l(1:ct(imov)+1))+1; en = sum(l(1:ct(imov+1)+1));
        aggr_ts(tussls{igen}.ind(st:en)',imov) = aggr_ts(tussls{igen}.ind(st:en)',imov) + 1;
        
        l = [0 chases{igen}.obj1.len];
        st = sum(l(1:cch1(imov)+1))+1; en = sum(l(1:cch1(imov+1)+1));
        chas_ts(chases{igen}.obj1.lab(st:en)',imov) = chas_ts(chases{igen}.obj1.lab(st:en)',imov) + 1;
        l = [0 chases{igen}.obj2.len];
        st = sum(l(1:cch2(imov)+1))+1; en = sum(l(1:cch2(imov+1)+1));
        chas_ts(chases{igen}.obj2.lab(st:en)',imov) = chas_ts(chases{igen}.obj2.lab(st:en)',imov) + 1;
        
        l = [0 wings.ext.l{igen}.obj1.len];
        st = sum(l(1:cwextl1(imov)+1))+1; en = sum(l(1:cwextl1(imov+1)+1));
        cour_ts(wings.ext.l{igen}.obj1.lab(st:en)',imov) = cour_ts(wings.ext.l{igen}.obj1.lab(st:en)',imov) + 1;
        l = [0 wings.ext.l{igen}.obj2.len];
        st = sum(l(1:cwextl2(imov)+1))+1; en = sum(l(1:cwextl2(imov+1)+1));
        cour_ts(wings.ext.l{igen}.obj2.lab(st:en)',imov) = cour_ts(wings.ext.l{igen}.obj2.lab(st:en)',imov) + 1;
        l = [0 wings.ext.r{igen}.obj1.len];
        st = sum(l(1:cwextr1(imov)+1))+1; en = sum(l(1:cwextr1(imov+1)+1));
        cour_ts(wings.ext.r{igen}.obj1.lab(st:en)',imov) = cour_ts(wings.ext.r{igen}.obj1.lab(st:en)',imov) + 1;
        l = [0 wings.ext.r{igen}.obj2.len];
        st = sum(l(1:cwextr2(imov)+1))+1; en = sum(l(1:cwextr2(imov+1)+1));
        cour_ts(wings.ext.r{igen}.obj2.lab(st:en)',imov) = cour_ts(wings.ext.r{igen}.obj2.lab(st:en)',imov) + 1;
        l = [0 courts{igen}.obj1.len];
        st = sum(l(1:ccrt1(imov)+1))+1; en = sum(l(1:ccrt1(imov+1)+1));
        cour_ts(courts{igen}.obj1.lab(st:en)',imov) = cour_ts(courts{igen}.obj1.lab(st:en)',imov) + 1;
        l = [0 courts{igen}.obj2.len];
        st = sum(l(1:ccrt2(imov)+1))+1; en = sum(l(1:ccrt2(imov+1)+1));
        cour_ts(courts{igen}.obj2.lab(st:en)',imov) = cour_ts(courts{igen}.obj2.lab(st:en)',imov) + 1;
    end
    
    % SUM ACTIONS INTO dtb BINS
    nbins = ceil(size(aggr_ts,1)/binsz);
    aggr = zeros(nbins,imov); chas = aggr; cour = aggr;
    for iii=1:nbins,
        en = min(iii*binsz,length(aggr_ts));
        aggr(iii,:) = sum(aggr_ts((iii-1)*binsz+1:en,:),1);
        chas(iii,:) = sum(chas_ts((iii-1)*binsz+1:en,:),1);
        cour(iii,:) = sum(cour_ts((iii-1)*binsz+1:en,:),1);
        
        fprintf(FIDT,'%11g',iii*dtb);
        for imov=1:nmovs(igen),
            fprintf(FIDT,' %6.2f %6.2f %6.2f ',...
                aggr(iii,imov)*dt,chas(iii,imov)*dt,cour(iii,imov)*dt);
        end
        fprintf(FIDT,'\n');
    end
    
    % MEAN AND S.E.M. OF OCCURRENCES
    aggress.mean = mean(aggr,2); aggress.sem = std(aggr,[],2)/sqrt(nmovs(igen));
    chasing.mean = mean(chas,2); chasing.sem = std(chas,[],2)/sqrt(nmovs(igen));
    courting.mean = mean(cour,2); courting.sem = std(cour,[],2)/sqrt(nmovs(igen));
    
    % PERCENTAGES
    sum_all = aggress.mean + chasing.mean + courting.mean;
    sum_all(~sum_all) = 1;
    aggress.percent = aggress.mean./sum_all*100;
    chasing.percent = chasing.mean./sum_all*100;
    courting.percent = courting.mean./sum_all*100;
%     aggress.percent = aggress.mean./binsz*100;
%     chasing.percent = chasing.mean./binsz*100;
%     courting.percent = courting.mean./binsz*100;
    
    % PLOTS
    h = subplot('Position',[mod(icnt-1,k2)/(k2*1.1)+.05, ...
        0.9-fix((icnt-1)/k2)/(k1*1.1)-1/(k1*1.5), ...
        1/(k2*1.5), 1/(k1*1.5)]); hold on;
    bar([aggress.percent chasing.percent courting.percent],'stacked','Edgecolor','none');
    barcolors = [.8 .2 .2 ; .2 .8 .2 ; .2 .2 .8]; colormap(barcolors);
%     barcolors = [params.flycol.win ; params.flycol.los ; params.flycol.ext]; colormap(barcolors);
%     barcolors = [0 0 0 ; .8 .8 .8 ; .3 .3 .3]; colormap(barcolors);
    set(h,'FontSize',params.axisfontsize);
    dtick = 120; % [s]
    t = 0:dtick/dtb:params.max_frames/dtb;
    set(gca,'XTick',t); xtick_lab = cell(1,length(t));
    for i=1:length(t), xtick_lab{i} = num2str(t(i)*dtb/60); end
    set(gca,'XTickLabel',xtick_lab);
    if ~mod(icnt-1,k2), ylabel('Percentage'); end
    if (icnt > (k1-1)*k2) || (icnt == 1 && igen == ngens), 
        xlabel('time [min]'); 
    end
    axis([0 params.max_frames/dtb 0 100]); axis tight;
    tit = gen_ident{igen}; in = strfind(tit,'_'); tit(in) = '-';
    title(tit);
    if ~mod(icnt,nplt) || igen == ngens,
        h = legend('Aggression','Chasing','Courtship');
        titl = ['Percentage of actions over all actions per ' ...
            num2str(dtb/60,'%4.2f') ' minute bin'];
        set(FID,'Name',titl); mtit(titl,'FontSize',params.axisfontsize+2,'yoff',.03);
        if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
    end
    icnt = icnt + 1;
end
