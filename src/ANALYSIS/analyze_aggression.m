%% Analyze_aggression
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
%% Analyze detected aggressive actions (lunging, tussling, wing threats) 
%% and compute various statistics

function FID = analyze_aggression(obj1,obj2,lunges,wings,chases,courts,tussls,copu,nmovs,dts,nframes,gen_ident,gen_ident1,params,chamber,FID)

if ~params.courtship,
    ngens = length(lunges);
    
    %% DENDROGRAMS from various position histograms (lunges, tussling, wing threat,
    %% chasing, circling)
    if (ngens >= 3) && params.plots.stat.dendro,
        minx = -chamber.width/2; maxx = chamber.width/2; miny = -chamber.height/2; maxy = chamber.height/2; dv = 2;
        xy = cell(ngens,1); lg = xy; tus = xy; wthr = xy; chs = xy; crt = xy;
        xv = minx:dv:maxx; yv = miny:dv:maxy;
        for igen=1:ngens,
            xy{igen} = hist3([[obj1{igen}.x obj2{igen}.x]',[obj1{igen}.y obj2{igen}.y]'],{xv yv})/nmovs(igen) + 1;
            if numel(lunges{igen}.obj1.x1) || numel(lunges{igen}.obj2.x1),
                lg{igen} = hist3([[lunges{igen}.obj1.x1 lunges{igen}.obj2.x1]', ...
                    [lunges{igen}.obj1.y1 lunges{igen}.obj2.y1]'],{xv yv})/nmovs(igen) + 1;
            else
                lg{igen} = ones(length(xv),length(yv));
            end
            if numel(tussls{igen}.x),
                tus{igen} = hist3([tussls{igen}.x', tussls{igen}.y'],{xv yv})/nmovs(igen) + 1;
            else
                tus{igen} = ones(length(xv),length(yv));
            end
            if numel(wings.threat{igen}.obj1.x1) || numel(wings.threat{igen}.obj2.x1),
                wthr{igen}  = hist3([[wings.threat{igen}.obj1.x1 wings.threat{igen}.obj2.x1]', ...
                    [wings.threat{igen}.obj1.y1 wings.threat{igen}.obj2.y1]'],{xv yv})/nmovs(igen) + 1;
            else
                wthr{igen} = ones(length(xv),length(yv));
            end
            if numel(chases{igen}.obj1.x1) || numel(chases{igen}.obj2.x1),
                chs{igen} = hist3([[chases{igen}.obj1.x1 chases{igen}.obj2.x1]', ...
                    [chases{igen}.obj1.y1 chases{igen}.obj2.y1]'],{xv yv})/nmovs(igen) + 1;
            else
                chs{igen} = ones(length(xv),length(yv));
            end
            if numel(courts{igen}.obj1.x1) || numel(courts{igen}.obj2.x1),
                crt{igen} = hist3([[courts{igen}.obj1.x1 courts{igen}.obj2.x1]', ...
                    [courts{igen}.obj1.y1 courts{igen}.obj2.y1]'],{xv yv})/nmovs(igen) + 1;
            else
                crt{igen} = ones(length(xv),length(yv));
            end
        end

        % Put all information into one matrix
        len = length(xv)*length(yv);
        X = zeros(ngens,len*5);
        for igen=1:ngens,
            %         X(igen,:) = [reshape(xy{igen},1,len) reshape(lg{igen},1,len) ...
            %             reshape(chs{igen},1,len)];
            X(igen,:) = [reshape(lg{igen},1,len) reshape(tus{igen},1,len) reshape(wthr{igen},1,len) ...
                reshape(chs{igen},1,len) reshape(crt{igen},1,len)];
        end
        
        % Try to normalize all features, so that they all have an equal weight
        % Z = X;
        % med = median(X); st = std(X); med = repmat(med,[ngens 1]); st = repmat(st,[ngens 1]);
        % X = (X - med) ./ st; X(isnan(X)) = 0;
        X(:,1:len) = (X(:,1:len) - mean(median(X(:,1:len)))) / mean(std(X(:,1:len)));
        X(:,len+1:2*len) = (X(:,len+1:2*len) - mean(median(X(:,len+1:2*len)))) / mean(std(X(:,len+1:2*len)));
        X(:,2*len+1:3*len) = (X(:,2*len+1:3*len) - mean(median(X(:,2*len+1:3*len)))) / mean(std(X(:,2*len+1:3*len)));
        X(:,3*len+1:4*len) = (X(:,3*len+1:4*len) - mean(median(X(:,3*len+1:4*len)))) / mean(std(X(:,3*len+1:4*len)));
        X(:,4*len+1:5*len) = (X(:,4*len+1:5*len) - mean(median(X(:,4*len+1:5*len)))) / mean(std(X(:,4*len+1:5*len)));
        % Sphere data as one way to normalize them
        % [U,S,V] = svd(X);
        % figure(12); set(12,'PaperOrientation','portrait','PaperPositionMode','manual','PaperPosition',[0 0 8.5 11]); set(12,'Name','log Eigenvalues'); title('log Eigenvalues','Fontsize',16);
        % loglog(diag(S)); axis equal;
        % X = U * S(:,1:ngens) * V(1:ngens,1:ngens);

        % Measure pairwise distances, i.e., the euclidean
        C = pdist(X,'euclidean');
        % figure(1); colormap('hot'); imagesc(squareform(C));
        % axis square; xlabel('sequence no.'); ylabel('sequence no.');
        % title('correlation-distance matrix'); colorbar;
        % Compute linkage and plot dendogram
        L = linkage(C,'complete');
        if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
        figure(FID); clf; [H,T,P] = dendrogram(L,'colorthreshold','default'); hold on;
        gen_id_perm = gen_ident(P); set_xtick_label(gen_id_perm,45,params.cross_to);
        titl = 'Dendrogram from Heatmaps';
        set(FID,'Name',titl); mtit(titl,'FontSize',params.axisfontsize+2,'yoff',.04);
        if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
    end

    %% TUSSLING HEATMAPS and RASTER PLOTS
    mint = 0; maxt = (round(nframes{1}(1)*dts{1}(1)/60)+1)*60; dt1 = 60;
    minx = -chamber.width/2; maxx = chamber.width/2;
    miny = -chamber.height/2; maxy = chamber.height/2; dv = (maxx-minx)/30;
    xt = mint:dt1:maxt; xv = minx:dv:maxx; yv = miny:dv:maxy;
    mea_img = zeros(numel(yv),numel(xv),ngens); sem_img = mea_img;
    distr = cell(1,ngens); mov = distr;
    for igen=1:ngens,
        if params.tim, dt = sum(nframes{igen}.*dts{igen})/sum(nframes{igen}); else dt = 1; end
        nmov = length(tussls{igen}.number); img = zeros(numel(xv),numel(yv),nmov);
        cumind = [0 cumsum(tussls{igen}.number)]; cumindl = [1 tussls{igen}.len];
        for j=1:nmov,
            if sum(tussls{igen}.number(j)) > 0,
                % Compute a 2d histogram of fly positions while tussling for each movie
                st = sum(cumindl(1:cumind(j)+1)); en = sum(cumindl(2:cumind(j+1)+1));
                img(:,:,j) = hist3([tussls{igen}.x(st:en)' , tussls{igen}.y(st:en)'], ...
                    {minx:dv:maxx miny:dv:maxy})*dt;
                % Prepare time stamps (distr) and movie/fly pair numbers (mov) for raster
                % plots
                st = cumind(j)+1; en = cumind(j+1);
                distr{igen} = [distr{igen} tussls{igen}.t(st:en)/60];
                mov{igen} = [mov{igen} ones(1,numel(st:en))*j];
            end
        end
        % Compute mean and sem (not plotted) 2d histogram of fly positions while tussling
        mea_img(:,:,igen) = mean(img,3)'; 
        sem_img(:,:,igen) = std(img,[],3)'/sqrt(nmov);
    end

    % Tussling plots
    if params.plots.heatmaps.tussling,
        % Auto scaling?
        if params.autoscale,
            max_val = ceil(max(max(max(mea_img))));
            max_val = max_val+ceil(.1*abs(max_val)); % extend the display range by 10%
        else
            if params.tim,
                max_val = params.range.heatmaps.tim.tussl;
            else
                max_val = params.range.heatmaps.occ.tussl;
            end
        end
        % TUSSLING HEATMAPS (one per genotype)
        if numel(find(params.plots.heatmaps.tussl_vec == 1)),
            if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
            figure(FID); clf;
            set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
            for igen=1:ngens,
                h = subplot('Position',[mod(igen-1,params.k(2))/(params.k(2)*1.1)+.05, 0.9-fix((igen-1)/params.k(2))/(params.k(1)*1.1)-1/(params.k(1)*1.5), ...
                    1/(params.k(2)*1.5), 1/(params.k(1)*1.5)]); hold on;
                xlim([-chamber.width/2 chamber.width/2]); ylim([-chamber.height/2 chamber.height/2]);
                set(h,'FontSize',params.axisfontsize);
                title(cell2mat(gen_ident(igen)),'FontSize',params.axisfontsize); axis image;
                imagesc(xv,yv,mea_img(:,:,igen),[0 max_val]); drawnow;
                colormap('default'); colormap([0 0 0 ; jet]);
            end
            pos = get(gca,'Position');
            ch = colorbar('location','EastOutside');
            if params.tim,
                set(get(ch,'YLabel'),'String','cumulative time [s]');
            else
                set(get(ch,'YLabel'),'String','occurence');
            end
            set(ch,'FontSize',params.axisfontsize);
            set(ch,'Position',[0.92 pos(2) 0.02 pos(4)]);
            titl = 'Histograms of Fly Tussling Positions (MEAN)';
            set(FID,'Name',titl); mtit(titl,'FontSize',params.axisfontsize+2,'yoff',.04);
            if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
        end
        
        % TUSSLING RASTER PLOTS (one per genotype)
        if numel(find(params.plots.heatmaps.tussl_vec == 2)),
            COL =  .8*[[.87 .76 .54]*.6 ; .87 .76 .54]; COL2 = .8*[[.71 .82 .49]*.6 ; .71 .82 .49];
            if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
            figure(FID); clf;
            set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
            stpx = 5; xaxis = 0:stpx:params.max_frames/60;
            xtick_lab = cell(1,length(xaxis));
            for i=1:length(xaxis), xtick_lab{i} = num2str((i-1)*stpx); end
            miny = 0; maxy = 1; yaxis = miny:.1:maxy;
            ytick_lab = cell(1,length(yaxis)); hmax = 0;
            for igen=1:ngens,
                nmov0 = length(lunges{igen}.obj1.number); nmov = nmov0 * 1.5;
                h = hist(distr{igen},0:1:30); h = h/nmov0;
                if hmax < max(h), hmax = max(h); end
            end
            fac = 1/ceil(hmax/.3);
            for i=1:length(yaxis),
                if i < 5,
                    ytick_lab{i} = num2str(yaxis(i)/fac);
                else
                    ytick_lab{i} = ' ';
                end
            end
            for igen=1:ngens,
                h = subplot('Position',[mod(igen-1,params.k(2))/(params.k(2)*1.1)+.05, 0.9-fix((igen-1)/params.k(2))/(params.k(1)*1.1)-1/(params.k(1)*1.5), ...
                    1/(params.k(2)*1.5), 1/(params.k(1)*1.5)]);
                set(h,'FontSize',params.axisfontsize); hold on;
                nmov0 = length(lunges{igen}.obj1.number); nmov = nmov0 * 1.5;
                dw = 1/nmov; db = .8*dw/2;
                % Plot the raster lines
                for i=1:numel(distr{igen}),
                    if mod(i,2),
                        if mod(mov{igen}(i),2), c = COL(1,:); else c = COL2(1,:); end
                    else
                        if mod(mov{igen}(i),2), c = COL(2,:); else c = COL2(2,:); end
                    end
                    plot([distr{igen}(i) distr{igen}(i)],dw/2+[(mov{igen}(i)+nmov-nmov0)/nmov-dw-db (mov{igen}(i)+nmov-nmov0)/nmov-dw+db],'Color',c); hold on;
                end
                % Plot the histogram
                h = hist(distr{igen},0:1:30); h = h/nmov0*fac;
                dbar = .2;
                for ii=0:max(xaxis),
                    patch([ii-dbar ii+dbar ii+dbar ii-dbar ii-dbar],...
                        [0 0 h(ii+1) h(ii+1) 0],'k');
                end
                set(gca,'YTick',yaxis); set(gca,'YTickLabel',ytick_lab); hold on;
                set(gca,'XTick',xaxis); set(gca,'XTickLabel',xtick_lab);
                colormap('default'); colormap([0 0 0 ; jet]);
                if ~mod(igen-1,params.k(2)), ylabel('occurrence  |  fly pair -->'); end
                if (igen > (params.k(1)-1)*params.k(2)), xlabel('time [min]'); end
                title(cell2mat(gen_ident(igen)),'FontSize',params.axisfontsize);
                axis([xaxis(1) xaxis(end) yaxis(1) yaxis(end)]);                
            end
            titl = 'Fly Tussling Distribution';
            set(FID,'Name',titl);
            mtit(titl,'FontSize',params.axisfontsize+2,'yoff',.04);
            if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
        end
    end

    % LUNGING HEATMAPS and RASTER PLOTS
    if params.plots.heatmaps.lunging,
        dxy.x = 0; dxy.y = dxy.x;
        % Plot the positions of the opponent as well?
        plot_oppon = 0;
        % Cumulative time or occurence
        if params.tim,
            max_val = params.range.heatmaps.tim.lunge;
        else
            max_val = params.range.heatmaps.occ.lunge;
        end

        % LUNGING HEATMAPS
        if numel(find(params.plots.heatmaps.lunge_vec == 1)),
            titl = 'Histograms of Fly Lunge Positions';
            FID = plot_heatmaps(lunges,plot_oppon,titl,max_val,0,0,chamber,gen_ident,dxy,dts,nframes,FID,params);
        end

        % LUNGING RASTER PLOTS
        if numel(find(params.plots.heatmaps.lunge_vec == 2)),
            titl = 'Fly Lunge Distribution';
            FID = plot_2ddistr(lunges,titl,max_val,gen_ident,dxy,dts,nframes,FID,params);
        end
    end

    % WING THREAT HEATMAPS and RASTER PLOTS
    if params.plots.heatmaps.wingthreat,
        dxy.x = 0; dxy.y = dxy.x;
        % Plot the positions of the opponent as well?
        plot_oppon = 0;
        % Cumulative time or occurence
        if params.tim,
            max_val = params.range.heatmaps.tim.wingthreat;
        else
            max_val = params.range.heatmaps.occ.wingthreat;
        end

        % WING THREAT HEATMAPS
        if numel(find(params.plots.heatmaps.wingthreat_vec == 1)),
            titl = 'Histograms of Wing Threat Positions';
            FID = plot_heatmaps(wings.threat,plot_oppon,titl,max_val,0,0,chamber,gen_ident,dxy,dts,nframes,FID,params);
        end

        % WING THREAT RASTER PLOTS
        if numel(find(params.plots.heatmaps.wingthreat_vec == 2)),
            titl = 'Wing Threat Distribution';
            FID = plot_2ddistr(wings.threat,titl,max_val,gen_ident,dxy,dts,nframes,FID,params);
        end
    end
    
    
    % LUNGING / TUSSLING BEHAVIOR ANALYSIS
    % MIN. LATENCIES TO THE FIRST LUNGE / TUSSLING
	% Compute min. latencies
    latency = zeros(1,ngens); latency_tussls = latency;
    for igen=1:ngens,
        data = tussls{igen}.t;
        if numel(data),
            latency_tussls(igen) = min(data);
        end
        if numel(lunges{igen}.t),
            latency(igen) = min(lunges{igen}.t);
        end
    end
    minv = 0;
	% Autoscale histograms?
    if params.autoscale, 
        maxvl = ceil(max(latency/60))*60; dtl = 0; 
        maxvt = ceil(max(latency_tussls/60))*60; dtt = 0;
        maxvl = maxvl+ceil(.1*abs(maxvl)); % extend the display range by 10%
        maxvt = maxvt+ceil(.1*abs(maxvt));
    else
        maxvl = params.max_frames; maxvt = params.max_frames; 
        dtl = 300; dtt = 300; % = 5 min.
    end
	
	% Plot histograms of min. latencies to the first lunge
    if params.plots.stat.lunging && ~params.plots.movieclips && numel(find(params.plots.stat.lunge_vec == 1)),
        maxv = maxvl; dt = dtl; x = 1:ngens;
        if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
        figure(FID); clf;
        set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
		% Plot bars
        patch_plot(x',latency',params);
        axis([0 ngens+1 minv maxv]);
        set(gca,'FontSize',params.axisfontsize);
        if dt,
            t = minv:dt:maxv; ytick_lab = cell(1,length(t));
            for i=1:length(t), ytick_lab{i} = num2str(t(i)/60); end
            set(gca,'YTick',minv:dt:maxv); set(gca,'YTickLabel',ytick_lab);
            ylabel('latency [min]');
        else
            ylabel('latency [s]');
        end
        set(gca,'XTick',x); set_xtick_label(gen_ident,45,params.cross_to);
        titl = 'Minimal Latency to the 1st Lunge';
        set(FID,'Name',titl); mtit(titl,'FontSize',params.axisfontsize+2,'yoff',.0);
        if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
    end

	% Plot histograms of min. latencies to the first tussling
    if params.plots.stat.tussl && ~params.plots.movieclips && numel(find(params.plots.stat.tussl_vec == 1)),
        maxv = maxvt; dt = dtt; x = 1:ngens;
        if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
        figure(FID); clf;
        set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
		% Plot bars
        patch_plot(x',latency_tussls',params);
        axis([0 ngens+1 minv maxv]);
        set(gca,'FontSize',params.axisfontsize);
        if dt,
            t = minv:dt:maxv; ytick_lab = cell(1,length(t));
            for i=1:length(t), ytick_lab{i} = num2str(t(i)/60); end
            set(gca,'YTick',minv:dt:maxv); set(gca,'YTickLabel',ytick_lab);
            ylabel('latency [min]');
        else
            ylabel('latency [s]');
        end
        set(gca,'XTick',x); set_xtick_label(gen_ident,45,params.cross_to);
        titl = 'Minimal Latency to the 1st Tussling';
        set(FID,'Name',titl); mtit(titl,'FontSize',params.axisfontsize+2,'yoff',.0);
        if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
    end

	% ANALSYS OF LATENCIES TO THE FIRST LUNGE AND TUSSLING
	% CUMULATIVE NUMBER OF LUNGES / TUSSLING OVER TIME
    utest = utestpairs(gen_ident1);
    minv = 0; maxv = params.max_frames; dt = 300; % dt = 5 min.
    t = minv:dt:maxv; xtick_lab = cell(1,length(t));
    for i=1:length(t), xtick_lab{i} = num2str(t(i)/60); end

    mint = 0; maxt = params.max_frames; dv = 5; it = mint:dv:maxt; 
    maxy = 0; maxy_t = 0;
    x = mint:dv:maxt; mid = round((maxt-mint)/dv*(5/6));
    c_mea = cell(ngens,1); c_sem = c_mea; c_mea_tussl = c_mea; 
	c_sem_tussl = c_mea; c_diff = c_mea;
    clear latency; latency.mea = zeros(ngens,1); latency.sem = latency.mea; 
	latency.max = latency.mea; latency.min = latency.mea;
    latency_tussl= latency;
    Xl = []; Gl = []; Xt = []; Gt = [];
    for igen=1:ngens,
        % Compute histograms & cumulative distributions of lunges (data & cdata) 
        % and tussling (data_tussl & cdata_tussl)
        imov = [0 cumsum(lunges{igen}.number)]; 
		imov_tussl = [0 cumsum(tussls{igen}.number)];
        data = zeros(length(x),nmovs(igen)); data_tussl = data;
        for i=1:length(imov)-1,
            if lunges{igen}.number(i) ~=0,
                data(:,i) = hist(lunges{igen}.t(imov(i)+1:imov(i+1)),mint:dv:maxt);
            end
        end
        for i=1:length(imov_tussl)-1,
            if tussls{igen}.number(i) ~=0,
                data_tussl(:,i) = hist(tussls{igen}.t(imov_tussl(i)+1:imov_tussl(i+1)),...
				mint:dv:maxt);
            end
        end
        cdata = cumsum(data,1); cdata_tussl = cumsum(data_tussl,1);

        % Compute mean and sem histogram & cumulative distribution 
        % of lunges (mea_h, sem_h & c_mea, c_sem)
        mea_h = mean(data,2); sem_h = std(data,[],2) / sqrt(length(data));
        c_mea{igen} = cumsum(mea_h); c_sem{igen} = cumsum(mea_h + sem_h) - c_mea{igen};
        if max(c_mea{igen}) > maxy, maxy = ceil(max(c_mea{igen})); end        
        % Determine minimum latency to the first lunge for each fly pair
        diff = cdata - 0.5; in = [];
        for i=1:length(imov)-1, 
            ind = find(diff(:,i) > 0); 
            if numel(ind), in = [in ind(1)]; end; 
        end
        % Compute mean, sem, min, max latency to the first lunge
        if ~isempty(in),
            latency.mea(igen) = mean(it(in)); latency.sem(igen) = std(it(in)) / sqrt(length(in));
            latency.min(igen) = min(it(in)); latency.max(igen) = max(it(in));
            Xl = [Xl it(in)]; Gl = [Gl ones(1,numel(in))*igen];
        else
            Xl = [Xt 0]; Gl = [Gl ones(1,numel(1))*igen];
        end

        % Compute mean and sem histogram & cumulative distribution 
        % of tussling (mea_h, sem_h & c_mea_tussl, c_sem_tussl)
        mea_h = mean(data_tussl,2); sem_h = std(data_tussl,[],2) / sqrt(length(data_tussl));
        c_mea_tussl{igen} = cumsum(mea_h); c_sem_tussl{igen} = cumsum(mea_h + sem_h) - c_mea_tussl{igen};
        if max(c_mea_tussl{igen}) > maxy_t, maxy_t = ceil(max(c_mea_tussl{igen})); end
        % Determine minimum latency to the first tussling for each fly pair
        diff = cdata_tussl - 0.5; in = [];
        for i=1:length(imov_tussl)-1,
            ind = find(diff(:,i) > 0);
            if numel(ind), in = [in ind(1)]; end
        end
        % Compute mean, sem, min, max latency to the first tussling
        if ~isempty(in),
            latency_tussl.mea(igen) = mean(it(in)); latency_tussl.sem(igen) = std(it(in)) / sqrt(length(in));
            latency_tussl.min(igen) = min(it(in)); latency_tussl.max(igen) = max(it(in));
            Xt = [Xt it(in)]; Gt = [Gt ones(1,numel(in))*igen];
        else
            Xt = [Xt 0]; Gt = [Gt ones(1,numel(1))*igen];
        end
        % Difference between cumulative distribution of lunging and tussling
        % Not plotted for far
        c_diff{igen} = (c_mea{igen} - c_mea_tussl{igen});
    end
    [s,gen] = sort(latency.mea);
    % No autoscale y-axis?
    if ~params.autoscale, 
        maxy = params.range.stat.occ.lunge*2; 
        maxy_t = params.range.stat.occ.tussl*2; 
    end

    % PLOTS
    if params.plots.stat.lunging && ~params.plots.movieclips,
        if numel(find(params.plots.stat.lunge_vec == 2)),
            % Plot mean cumulative number of lunges over time
            if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
            figure(FID); clf;
            set(FID,'PaperOrientation','portrait','PaperPositionMode','manual','PaperPosition',[0 0 8.5 11]);
            col = jet(ngens);
            h = subplot(3,1,1);
            for igen=1:ngens,
                semilogy(x,c_mea{gen(igen)},'Color',col(igen,:),'LineWidth',2); hold on;
                text(x(mid),c_mea{igen}(mid),['\leftarrow ' gen_ident{igen}],'Fontsize',6,'Color','k');
            end
            set(h,'FontSize',params.axisfontsize);
            set(gca,'XTick',minv:dt:maxv); set(gca,'XTickLabel',xtick_lab);
            legend(gen_ident(gen));
            grid on; ylim([0 maxy]); xlabel('time [min]');
            ylabel('MEAN cumulative number of lunges');
            title('MEAN total number of lunges over time per fly pair');
            axis([minv maxv .5 maxy]);
            
            % Plot sem cumulative number of lunges over time
            h = subplot(3,1,2);
            for igen=1:ngens,
                semilogy(x,c_sem{gen(igen)},'Color',col(igen,:),'LineWidth',2); hold on;
                %     text(x(mid),c_sem{igen}(mid),['\leftarrow ' gen_ident{igen}],'Fontsize',6,'Color','k');
            end
            set(h,'FontSize',params.axisfontsize);
            set(gca,'XTick',minv:dt:maxv); set(gca,'XTickLabel',xtick_lab);
            grid on; ylim([0 maxy]); xlabel('time [min]');
            ylabel('SEM of cumulative number of lunges');
            title('SEM of total number of lunges over time per fly pair');
            axis([minv maxv .5 maxy]);

            % Plot cumulative number of lunges over time per fly pair
            h = subplot(3,1,3);
            col = gray(3);
            data = zeros(length(x),sum(nmovs)); icnt = 0;
            for igen=1:ngens,
                imov = [0 cumsum(lunges{igen}.number)];
                for i=1:length(imov)-1,
                    icnt = icnt + 1;
                    if lunges{igen}.number(i) ~=0,
                        data(:,icnt) = hist(lunges{igen}.t(imov(i)+1:imov(i+1)),[mint:dv:maxt]);
                    end
                end
            end
            mea_h = mean(data'); sem_h = std(data') / sqrt(length(data));
            c_mea = cumsum(mea_h); c_sem = cumsum(mea_h + sem_h);
            semilogy(x,c_mea,'Color','k','LineWidth',5); hold on;
            semilogy(x,c_sem,'.','Color',col(2,:),'LineWidth',2); hold on;
            ind = find(c_mea >= .5);
            if numel(ind),
                text(x(ind(1)),7,'\downarrow first lunge');
                text(x(ind(1))+30,3,['t = ' num2str(x(ind(1))) ' s']);
            end
            set(h,'FontSize',params.axisfontsize);
            set(gca,'XTick',minv:dt:maxv); set(gca,'XTickLabel',xtick_lab);
            grid on; ylim([0 maxy]); xlabel('time [min]'); ylabel('cumulative number of lunges');
            title('cummulative number of lunges over time per fly pair'); axis([minv maxv .5 maxy]);
            titl = 'Lunge Behavior over Time';
            set(FID,'Name',titl); mtit(titl,'FontSize',params.axisfontsize+2,'yoff',.04);
            if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
        end

        if numel(find(params.plots.stat.lunge_vec == 3)),
            % Histogram of latency to the first lunge per genotype
            if params.autoscale,
                maxv = ceil(max(latency.max/60))*60; dt = 0;%maxv/5;
                maxv = maxv+ceil(.1*abs(maxv)); % extend the display range by 10%
            else
                maxv = maxt; dt = 300; % = 5 min.
            end
            if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
            figure(FID); clf;
            set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
            xval = 1:ngens;
            if ~params.boxplot,
                % Bars with mean, sem, min, max
                patch_plot(xval,latency.mea,params);
                d1 = 0; d2 = -d1; w = .2;
                for i=1:length(xval),
                    plot([xval(i)+d1 xval(i)+d1],[latency.mea(i)-latency.sem(i) latency.mea(i)+latency.sem(i)],'k','Linewidth',1);
                    plot([xval(i)+d1-w xval(i)+d1+w],[latency.mea(i)-latency.sem(i) latency.mea(i)-latency.sem(i)],'k','Linewidth',1);
                    plot([xval(i)+d1-w xval(i)+d1+w],[latency.mea(i)+latency.sem(i) latency.mea(i)+latency.sem(i)],'k','Linewidth',1);
                end
                plot(xval,latency.min,'+','MarkerEdgeColor','k','MarkerSize',4,'LineWidth',1); hold on;
                plot(xval,latency.max,'+','MarkerEdgeColor','k','MarkerSize',4,'LineWidth',1); hold on;
            else
                % Box plots
                boxplot(Xl,Gl,'whisker',params.whisker,'plotstyle','traditional','notch','marker'); hold on;
            end
            if isstruct(utest) && numel(Xl), plot_utests(utest,Xl,Gl,latency,maxt,params); end
            axis([0 ngens+1 mint maxv]);
            set(gca,'FontSize',params.axisfontsize);
            if dt,
                t = minv:dt:maxv; ytick_lab = cell(1,length(t));
                for i=1:length(t), ytick_lab{i} = num2str(t(i)/60); end
                set(gca,'YTick',minv:dt:maxv); set(gca,'YTickLabel',ytick_lab);
                ylabel('latency [min]');
            else
                ylabel('latency [s]');
            end
            set(gca,'XTick',1:ngens); set_xtick_label(gen_ident,45,params.cross_to,params.axisfontsize);
            titl = 'Latency to the 1st Lunge'; mtit(titl,'FontSize',params.axisfontsize+2,'yoff',.0);
            set(FID,'Name',titl);
            if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
        end
    end

    if params.plots.stat.tussl && ~params.plots.movieclips,
        % Histogram of latency to the first tussling per genotype
        if numel(find(params.plots.stat.tussl_vec == 2)),
            if params.autoscale,
                maxv = ceil(max(latency_tussl.max/60))*60; dt = 0;%maxv/5;
                maxv = maxv+ceil(.1*abs(maxv)); % extend the display range by 10%
            else
                maxv = maxt; dt = 300; % = 5 min.
            end
            if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
            figure(FID); clf;
            set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
            % subplot(2,1,2);
            xval = 1:ngens;
            if ~params.boxplot,
                % Bars with mean, sem, min, max
                patch_plot(xval,latency_tussl.mea,params);
                d1 = 0; d2 = -d1; w = .2; xval = 1:ngens;
                for i=1:length(xval),
                    plot([xval(i)+d1 xval(i)+d1],[latency_tussl.mea(i)-latency_tussl.sem(i) latency_tussl.mea(i)+latency_tussl.sem(i)],'k','Linewidth',1);
                    plot([xval(i)+d1-w xval(i)+d1+w],[latency_tussl.mea(i)-latency_tussl.sem(i) latency_tussl.mea(i)-latency_tussl.sem(i)],'k','Linewidth',1);
                    plot([xval(i)+d1-w xval(i)+d1+w],[latency_tussl.mea(i)+latency_tussl.sem(i) latency_tussl.mea(i)+latency_tussl.sem(i)],'k','Linewidth',1);
                end
                plot(xval,latency_tussl.min,'+','MarkerEdgeColor','k','MarkerSize',4,'LineWidth',1); hold on;
                plot(xval,latency_tussl.max,'+','MarkerEdgeColor','k','MarkerSize',4,'LineWidth',1); hold on;
            else
                % Box plots
                boxplot(Xt,Gt,'whisker',params.whisker,'plotstyle','traditional','notch','marker'); hold on;
            end
            if isstruct(utest) && numel(Xt), plot_utests(utest,Xt,Gt,latency_tussl,maxt,params); end
            axis([0 ngens+1 mint maxv]);
            set(gca,'FontSize',params.axisfontsize);
            set(gca,'XTick',1:ngens); set_xtick_label(gen_ident,45,params.cross_to,params.axisfontsize);
            if dt,
                t = minv:dt:maxv; ytick_lab = cell(1,length(t));
                for i=1:length(t), ytick_lab{i} = num2str(t(i)/60); end
                set(gca,'YTick',minv:dt:maxv); set(gca,'YTickLabel',ytick_lab);
                ylabel('latency [min]');
            else
                ylabel('latency [s]');
            end
            titl = 'Latency to the 1st Tussling'; mtit(titl,'FontSize',params.axisfontsize+2,'yoff',.0);
            set(FID,'Name',titl);
            if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
        end

        if numel(find(params.plots.stat.tussl_vec == 3)),
            % Plot mean cumulative number of tussling over time per fly pair
            if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
            figure(FID); clf;
            set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
            col = jet(ngens);
            for igen=1:ngens,
                semilogy(x,c_mea_tussl{gen(igen)},'Color',col(igen,:),'LineWidth',2); hold on;
                text(x(mid),c_mea_tussl{igen}(mid),['\leftarrow ' gen_ident{igen}],'Fontsize',6,'Color','k');
            end
            set(gca,'FontSize',params.axisfontsize);
            set(gca,'XTick',minv:dt:maxv); set(gca,'XTickLabel',xtick_lab);
            grid on; ylim([0 maxy_t]); xlabel('time [min]'); ylabel('MEAN cumulative number of tussls');
            legend(gen_ident(gen));
            axis([minv maxv 0 maxy_t]);
            titl = 'MEAN cummulative number of tussls over time per fly pair';
            mtit(titl,'FontSize',params.axisfontsize+2,'yoff',.1);
            set(FID,'Name',titl);
            if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
        end
    end

    % PERCENTAGE OF FLY PAIRS PERFORMING AT LEAST ONE LUNGE / TUSSLING
    utest = utestpairs(gen_ident1);
    min_num_events = 1;
    lungepairs.mea = zeros(ngens,1); tusslpairs.mea = lungepairs.mea;
    lungepairs.sem = zeros(ngens,1); tusslpairs.sem = lungepairs.sem;
    Xl = []; Gl = []; Xt = []; Gt = [];
    for igen=1:ngens,
        datal = lunges{igen}.number >= min_num_events;
        datat = tussls{igen}.number >= min_num_events;
        lungepairs.mea(igen) = mean(datal) * 100;
        lungepairs.sem(igen) = 0;
        tusslpairs.mea(igen) = mean(datat) * 100;
        tusslpairs.sem(igen) = 0;
        Xl = [Xl datal*100.]; Gl = [Gl ones(1,length(datal))*igen];
        Xt = [Xt datat*100.]; Gt = [Gt ones(1,length(datat))*igen];
    end

    if params.plots.stat.lunging && ~params.plots.movieclips && numel(find(params.plots.stat.lunge_vec == 4)),
        maxy = 100; ytit = 'fly pairs [%]'; titl = ['Percentage of Fly Pairs >=' num2str(min_num_events) ' Lunge'];
        FID = plot_bar_msem(lungepairs,0,Xl,Gl,maxy,ytit,titl,FID,gen_ident,params);
    end

    if params.plots.stat.tussl && ~params.plots.movieclips && numel(find(params.plots.stat.tussl_vec == 4)),
        maxy = 100; ytit = 'fly pairs [%]'; titl = ['Percentage of Fly Pairs >=' num2str(min_num_events) ' Tussl'];
        FID = plot_bar_msem(tusslpairs,0,Xl,Gl,maxy,ytit,titl,FID,gen_ident,params);
    end

    % PERCENTAGE OF FLY PAIRS WITH N NUMBER OF LUNGES
    if params.plots.stat.lunging && ~params.plots.movieclips && numel(find(params.plots.stat.lunge_vec == 5)),
        % Plot cumlative distribution of percentage of fly pairs over number of lunges
        % Fly pairs that performed ? the displayed number of lunges
        if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
        figure(FID); clf;
        set(FID,'PaperOrientation','portrait','PaperPositionMode','manual','PaperPosition',[0 0 8.5 11]);
        subplot(3,1,1);
        mint = 0; maxt = 20; dv = 1; maxy = 100;
        x = mint:dv:maxt; data = zeros(length(x),ngens); col = jet(ngens);
        for igen=1:ngens,
            data(:,igen) = hist(lunges{igen}.obj1.number + lunges{igen}.obj2.number,[mint:dv:maxt])/nmovs(igen);
            tmp = cumsum(data(:,igen))*100;
            plot(x,tmp,'Color',col(igen,:),'LineWidth',4); hold on;
        end
        mid = 3; %round((maxt-mint)/dv*(1/6));
        % for igen=1:ngens,
        %     tmp = cumsum(data(:,igen))*100;
        %     text(x(mid),tmp(mid),['\leftarrow ' gen_ident{igen}],'Fontsize',12,'Color','k'); hold on;
        % end
        grid on; ylim([0 maxy]); xlabel('# lunges'); ylabel('fly pairs [%]');
        title('% of fly pairs with <= number of lunges'); axis([mint maxt 0 maxy]);
        % legend(gen_ident);

        % Plot mean & sem cumlative distribution of percentage of fly pairs over number of lunges
        % Fly pairs that performed ? the displayed number of lunges
        subplot(3,1,2);
        mea_h = mean(data'); sem_h = std(data') / sqrt(length(data));
        c_mea = cumsum(mea_h)*100; c_sem = cumsum(mea_h + sem_h)*100;
        plot(x,c_mea,'Color','k','LineWidth',5); hold on;
        plot(x,c_sem,'.','Color',col(1,:),'LineWidth',2); hold on;
        grid on; ylim([0 maxy]); xlabel('# lunges'); ylabel('fly pairs [%]');
        title('% of fly pairs with <= number of lunges'); axis([mint maxt 0 maxy]); %legend(gen_ident);

        % Plot cumlative distribution of percentage of fly pairs over number of lunges
        % Fly pairs that performed ? the displayed number of lunges
        subplot(3,1,3);
        for igen=1:ngens,
            tmp = 100-cumsum(data(:,igen))*100;
            plot(x,tmp,'Color',col(igen,:),'LineWidth',4); hold on;
        end
        mid = 3; 
        grid on; ylim([0 maxy]); xlabel('# lunges'); ylabel('fly pairs [%]');
        title('% of fly pairs with >= number of lunges'); axis([mint maxt 0 maxy]);
        legend(gen_ident);
        titl = 'Lunge Behavior per Fly Pair';
        set(FID,'Name',titl); mtit(titl,'FontSize',params.axisfontsize,'yoff',.04);
        if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
    end

    %% LUNGE FREQUENCY OVER TIME
    if params.plots.stat.lunging && ~params.plots.movieclips && numel(find(params.plots.stat.lunge_vec == 6)),
        % Compute a histogram of lunge distribution over time for one-minute bins
        % for each fly pair and genotype, sum over both flies
        dv = 60; minv = dv/2; maxv = params.max_frames-dv/2; iv = minv:dv:maxv;
        h2dm = zeros(ngens,length(iv)); h2dsem = h2dm; txtdata = cell(ngens,1);
        for igen=1:ngens,
            nmov = length(lunges{igen}.mov); h2d0 = zeros(nmov,length(iv));
            cumind1 = [0 cumsum(lunges{igen}.obj1.number)]; cumind2 = [0 cumsum(lunges{igen}.obj2.number)];
            for j=1:nmov,
                if lunges{igen}.obj1.number(j)>0,
                    st = cumind1(j)+1; en = cumind1(j+1);
                    h2d0(j,:) = hist(lunges{igen}.obj1.t(st:en)',iv);
                end
                if lunges{igen}.obj2.number(j)>0,
                    st = cumind2(j)+1; en = cumind2(j+1);
                    h2d0(j,:) = h2d0(j,:) + hist(lunges{igen}.obj2.t(st:en)',iv);
                end
            end
            % Compute mean and sem lunge distribution over time per genotype
            h2dm(igen,:) = mean(h2d0,1); h2dsem(igen,:) = std(h2d0,[],1)/sqrt(nmov);
            txtdata{igen} = h2dm(igen,:);
        end

        % Text output
        textoutput(txtdata,gen_ident,['mean frequency [per ' num2str(dv) ' s]'],'Lunges',params,'abs. time [s]',iv+30);

        % Plot mean lunge distribution over time per genotype
        if params.autoscale, 
            maxy = ceil(max(max(h2dm))); 
            maxy = maxy+ceil(.1*abs(maxy)); % extend the display range by 10%
        else
            maxy = 10; 
        end
        iy = 1:maxy;
        minv = 0; maxv = params.max_frames; dt = 300; % dt = 5 min.
        t = minv:dt:maxv; xtick_lab = cell(1,length(t));
        for i=1:length(t), xtick_lab{i} = num2str(t(i)/60); end
        if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
        figure(FID); clf;
        set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
        for igen=1:ngens,
            h = subplot('Position',[mod(igen-1,params.k(2))/(params.k(2)*1.1)+.05, 0.9-fix((igen-1)/params.k(2))/(params.k(1)*1.1)-1/(params.k(1)*1.5), ...
                1/(params.k(2)*1.5), 1/(params.k(1)*1.5)]); hold on;
            title(cell2mat(gen_ident(igen)),'FontSize',params.axisfontsize);
%             patch_plot(iv',h2dm(igen,:)',params);        
            bar(iv,h2dm(igen,:),params.barwidth); colormap(params.flycol.win);
            axis([minv maxv 0 maxy]); set(gca,'XTick',minv:dt:maxv);
            set(gca,'XTick',minv:dt:maxv); set(gca,'XTickLabel',xtick_lab);
            set(gca,'YTick',iy);
            xlim([minv maxv]);
            if igen>2*params.k(2), xlabel('time [min]'); end
            if ~mod(igen-1,params.k(2)), ylabel(['frequency of lunges / ' num2str(dv) ' s']); end
        end
        titl = 'MEAN lunge frequency over time';
        set(FID,'Name',titl); mtit(titl,'FontSize',params.axisfontsize+2,'yoff',.04);
        if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
    end

    %% TUSSLING BOUT FREQUENCY OVER TIME
    if params.plots.stat.tussl && ~params.plots.movieclips,
        % Compute a histogram of tussling distribution over time and duration 
        % for one-minute bins for each fly pair and genotype, sum over both flies
        dv = 60; minv = dv/2; maxv = params.max_frames-dv/2; iv = minv:dv:maxv;
        db = .1; minb = db/2; maxb = 10-db/2; ib = minb:db:maxb;
        h2dm = zeros(ngens,length(iv),length(ib)); h2dsem = h2dm;
        txtdata = cell(ngens,1); txtdata_f = cell(ngens,1);
        for igen=1:ngens,
            nmov = length(tussls{igen}.mov); h2d0 = zeros(nmov,length(iv),length(ib));
            cumind = [0 cumsum(tussls{igen}.number)];
            for j=1:nmov,
                if numel(tussls{igen}.number(j)),
                    st = cumind(j)+1; en = cumind(j+1);
                    h2d0(j,:,:) = hist3([tussls{igen}.t(st:en)', ...
                        tussls{igen}.len(st:en)'*dts{igen}(j)],{iv ib});
                end
            end
            % Compute mean and sem tussling distribution over time and duration per genotype
            h2dm(igen,:,:) = mean(h2d0,1); h2dsem(igen,:,:) = std(h2d0,[],1)/sqrt(nmov);
        end

        % Compute mean and sem tussling distribution over time per genotype
        % tussling bout duration (mea), bout frequency (fmea)
        y.mea = zeros(ngens,numel(iv)); y.sem = y.mea; y.min = y.mea; y.max = y.mea;
        y.fmea = y.mea;
        for igen=1:ngens,
            data = reshape(h2dm(igen,:,:),numel(iv),numel(ib));
            sdata = sum(data,2); sdata(sdata == 0) = 0.001;
            y.mea(igen,:) = data * ib' ./ sdata;
            y.fmea(igen,:) = sum(data,2);
            txtdata{igen} = y.mea(igen,:); 
            txtdata_f{igen} = y.fmea(igen,:);
        end

        % Plots
        minv = 0; maxv = params.max_frames; dt = 120; % dt = 5 min.
        t = minv:dt:maxv; xtick_lab = cell(1,length(t));
        for i=1:length(t), xtick_lab{i} = num2str(t(i)/60); end
        % Plot mean tussling bout duration distribution over time per genotype
        if numel(find(params.plots.stat.tussl_vec == 5)),
            % Text output
            textoutput(txtdata,gen_ident,'mean bout period [s]','Tussling',params,'abs. time [s]',iv+30);

            if params.autoscale, 
                maxy = ceil(max(max(y.mea))); 
                maxy = maxy+ceil(.1*abs(maxy)); % extend the display range by 10%
            else
                maxy = 10; 
            end
            if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
            figure(FID); clf;
            set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
            for igen=1:ngens,
                h = subplot('Position',[mod(igen-1,params.k(2))/(params.k(2)*1.1)+.05, 0.9-fix((igen-1)/params.k(2))/(params.k(1)*1.1)-1/(params.k(1)*1.5), ...
                    1/(params.k(2)*1.5), 1/(params.k(1)*1.5)]); hold on;
                title(cell2mat(gen_ident(igen)),'FontSize',params.axisfontsize);
                %             patch_plot(iv',y.mea(igen,:)',params);
                bar(iv,y.mea(igen,:),params.barwidth); colormap(params.flycol.win);
                %     plot(iv,y.min(igen,:),'+','MarkerEdgeColor','r','MarkerSize',4,'LineWidth',2); hold on;
                %     plot(iv,y.max(igen,:),'+','MarkerEdgeColor','r','MarkerSize',4,'LineWidth',2); hold on;
                axis([minv maxv 0 maxy]); set(gca,'XTick',minv:dt:maxv);
                set(gca,'XTick',minv:dt:maxv); set(gca,'XTickLabel',xtick_lab);
                xlim([minv maxv]);
                if igen>2*params.k(2), xlabel('time [min]'); end
                if ~mod(igen-1,params.k(2)), ylabel('bout [s]'); end
            end
            titl = 'MEAN tussling bout period over time';
            set(FID,'Name',titl); mtit(titl,'FontSize',params.axisfontsize+2,'yoff',.04);
            if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
        end
        
        % Plot mean tussling bout frequency distribution over time per genotype
        if numel(find(params.plots.stat.tussl_vec == 6)),
             % Text output
            textoutput(txtdata_f,gen_ident,['mean bout frequency [per ' num2str(dv) ' s]'],'Tussling',params,'abs. time [s]',iv+30);

            if params.autoscale, 
                maxy = ceil(max(max(y.fmea))); 
                maxy = maxy+ceil(.1*abs(maxy)); % extend the display range by 10%
            else
                maxy = 2; 
            end
            iy = 1:maxy;
            if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
            figure(FID); clf;
            set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
            for igen=1:ngens,
                h = subplot('Position',[mod(igen-1,params.k(2))/(params.k(2)*1.1)+.05, 0.9-fix((igen-1)/params.k(2))/(params.k(1)*1.1)-1/(params.k(1)*1.5), ...
                    1/(params.k(2)*1.5), 1/(params.k(1)*1.5)]); hold on;
                title(cell2mat(gen_ident(igen)),'FontSize',params.axisfontsize);
                %             patch_plot(iv',y.fmea(igen,:)',params);
                bar(iv,y.fmea(igen,:),params.barwidth); colormap(params.flycol.win);
                axis([minv maxv 0 maxy]); set(gca,'XTick',minv:dt:maxv);
                set(gca,'XTick',minv:dt:maxv); set(gca,'XTickLabel',xtick_lab);
                set(gca,'YTick',iy);
                xlim([minv maxv]);
                if igen>2*params.k(2), xlabel('time [min]'); end
                if ~mod(igen-1,params.k(2)), ylabel(['frequency of bouts / ' num2str(dv) ' s']); end
            end
            titl = 'MEAN tussling bout frequency over time';
            set(FID,'Name',titl); mtit(titl,'FontSize',params.axisfontsize+2,'yoff',.04);
            if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
        end
    end

    %% LUNGING HISTOGRAMS
    if params.plots.stat.lunging && ~params.plots.movieclips && numel(find(params.plots.stat.lunge_vec == 7)),
        % Plot lunges per genotype and for fly 1 & 2 or loser & winner
        bool_sumobj = 0;
        titl = 'Lunges'; ylab = 'lunges';
        if params.tim,
            maxy = params.range.stat.tim.lunge;
        else
            maxy = params.range.stat.occ.lunge;
        end
        FID = plot_loser_winner_stat(lunges,FID,maxy,gen_ident,params,ylab,dts,titl,bool_sumobj,lunges,0);

        % Plot lunges per genotype both flies together (bool_sumobj = 1)
        utest = utestpairs(gen_ident1);
        bool_sumobj = 1;
        titl = 'Lunges'; ylab = 'lunges';
        maxy = 2*maxy;
        FID = plot_loser_winner_stat(lunges,FID,maxy,gen_ident,params,ylab,dts,titl,bool_sumobj,lunges,utest);
    end

    %% WING THREAT HISTOGRAMS
    if params.plots.stat.wingthreat && ~params.plots.movieclips,
        % Plot wing threats per genotype and for fly 1 & 2 or loser & winner
        bool_sumobj = 0;
        titl = 'Wing Threats'; ylab = 'wing threats';
        if params.tim,
            maxy = params.range.stat.tim.wingthreat;
        else
            maxy = params.range.stat.occ.wingthreat;
        end
        FID = plot_loser_winner_stat(wings.threat,FID,maxy,gen_ident,params,ylab,dts,titl,bool_sumobj,lunges,0);

        % Plot wing threats per genotype both flies together (bool_sumobj = 1)
        utest = utestpairs(gen_ident1);
        bool_sumobj = 1;
        titl = 'Wing Threats'; ylab = 'wing threats';
        maxy = 2*maxy;
        FID = plot_loser_winner_stat(wings.threat,FID,maxy,gen_ident,params,ylab,dts,titl,bool_sumobj,lunges,utest);
    end


    %% TUSSLING NUMBER & BOUT PERIODS
    if params.plots.stat.tussl && ~params.plots.movieclips,
        y.mea = zeros(ngens,1); y.sem = y.mea; y.min = y.mea; y.max = y.mea;
        data = cell(ngens,1); sdata = data;
        X = []; G = [];
        if numel(find(params.plots.stat.tussl_vec == 7)),
            % PLOTS
            utest = utestpairs(gen_ident1);
            if ~params.tim,
                % Plot number of tussling per genotype
                maxy = params.range.stat.occ.tussl;
                for igen=1:ngens,
                    if numel(tussls{igen}.number),
                        X = [X tussls{igen}.number]; G = [G ones(1,length(tussls{igen}.number))*igen];
                        y.mea(igen) = mean(tussls{igen}.number);
                        y.sem(igen) = std(tussls{igen}.number) / sqrt(numel(tussls{igen}.number));
                        y.min(igen) = min(tussls{igen}.number); y.max(igen) = max(tussls{igen}.number);
                    end
                    data{igen} = [data{igen} tussls{igen}.number]; sdata{igen} = [sdata{igen} tussls{igen}.len*dt/60];
                end
                if params.bool_xls, textoutput(data,gen_ident,'number','Tussling',params); end
                ytit = 'number of tussling'; titl = 'Number of Tussling';
                FID = plot_bar_msem(y,utest,X,G,maxy,ytit,titl,FID,gen_ident,params);
            else
                % Plot tussling duration per genotype
                maxy = params.range.stat.tim.tussl;
                for igen=1:ngens,
                    dt = sum(dts{igen}.*nframes{igen})/sum(nframes{igen});
                    data = tussls{igen}.len*dt;
                    if numel(data),
                        X = [X data]; G = [G ones(1,length(data))*igen];
                        y.mea(igen) = mean(data);
                        y.sem(igen) = std(data) / sqrt(numel(data));
                        y.min(igen) = min(data); y.max(igen) = max(data);
                    end
                end
                if params.bool_xls, textoutput(data,gen_ident,'time spent [min]','Tussling',params); end
                ytit = 'tussling bout period [s]'; titl = 'Tussling Bout Periods';
                FID = plot_bar_msem(y,utest,X,G,maxy,ytit,titl,FID,gen_ident,params);
            end
        end
    end
    

    %%% LUNGES / TUSSLS PER METER
    if (params.plots.stat.lunging || params.plots.stat.tussl) && ~params.plots.movieclips,
        % Compute travel distance and travel time for normalization
        data = win_los_data('movdists','max',obj1,obj2,lunges,copu.pre,dts,params);
        datat= win_los_data('t','numel',obj1,obj2,lunges,copu.pre,dts,params);
        
        x = 1:ngens; yt.mea = zeros(ngens,1); yt.sem = yt.mea;
        lunges_per_m.mea = yt.mea; tussls_per_m.mea = yt.mea;
        lunges_per_m.sem = yt.mea; tussls_per_m.sem = yt.mea;
        lunges_per_min.mea = yt.mea; tussls_per_min.mea = yt.mea;
        lunges_per_min.sem = yt.mea; tussls_per_min.sem = yt.mea;
        tmpl = cell(ngens,1); tmpt = tmpl; tmplt = tmpl; tmptt = tmpl; 
        Xl = []; Xt = []; Gl = []; Gt = []; Xlt = []; Xtt = []; Glt = []; Gtt = [];
        for igen=1:ngens,
            if params.oneobj, data{igen}.los = 0; end
            % Mean travel distance and mean travel time
            tmp = (data{igen}.win + data{igen}.los) / 1000/2;
            tim = (datat{igen}.win + datat{igen}.los) .*dts{igen} / 60/2;
            % Lunges & tussling per meter (tmpl,tmpt) and per minute (tmplt,tmptt)
            tmpl{igen} = lunges{igen}.number ./ tmp;
            tmpt{igen} = tussls{igen}.number ./ tmp;
            tmplt{igen} = lunges{igen}.number ./ tim;
            tmptt{igen} = tussls{igen}.number ./ tim;

            % Prepare data for box plots
            Xl = [Xl tmpl{igen}]; Gl = [Gl ones(1,length(tmpl{igen}))*igen];
            Xt = [Xt tmpt{igen}]; Gt = [Gt ones(1,length(tmpt{igen}))*igen];
            Xlt = [Xlt tmplt{igen}]; Glt = [Glt ones(1,length(tmplt{igen}))*igen];
            Xtt = [Xtt tmptt{igen}]; Gtt = [Gtt ones(1,length(tmptt{igen}))*igen];

            % Compute mean, sem
            lunges_per_m.mea(igen) = mean(tmpl{igen});
            lunges_per_m.sem(igen) = std(tmpl{igen})/sqrt(numel(tmpl{igen}));
            lunges_per_min.mea(igen) = mean(tmplt{igen});
            lunges_per_min.sem(igen) = std(tmplt{igen})/sqrt(numel(tmplt{igen}));
            if numel(tussls{igen}.t),
                tussls_per_m.mea(igen) = mean(tmpt{igen});
                tussls_per_m.sem(igen) = std(tmpt{igen})/sqrt(numel(tmpt{igen}));
                tussls_per_min.mea(igen) = mean(tmptt{igen});
                tussls_per_min.sem(igen) = std(tmptt{igen})/sqrt(numel(tmptt{igen}));
            end
        end
        
        % TEXT OUTPUT OF LUNGES AND TUSSLING PER METER AND PER MINUTE
        if params.bool_xls,
            % Lunges per meter
            if params.plots.stat.lunging && numel(find(params.plots.stat.lunge_vec == 8)),
                textoutput(tmpl,gen_ident,'number [per meter]','Lunges',params);
            end
            % Tussling per meter
            if params.plots.stat.tussl && numel(find(params.plots.stat.tussl_vec == 8)),
                textoutput(tmpt,gen_ident,'number [per meter]','Tussling',params);
            end
            % Lunges per minute
            if  params.plots.stat.lunging && numel(find(params.plots.stat.lunge_vec == 9)),
                textoutput(tmplt,gen_ident,'number [per minute]','Lunges',params);
            end
            % Tussling per minute
            if params.plots.stat.tussl && numel(find(params.plots.stat.tussl_vec == 9)),
                textoutput(tmptt,gen_ident,'number [per minute]','Tussling',params);
            end
        end
        
        % Determine genotype pairs to be compared
        utest = utestpairs(gen_ident1);

        % Plots (bars with mean/sem or box plots) including statistics
        if params.plots.stat.lunging && numel(find(params.plots.stat.lunge_vec == 8)),
            ytit = 'lunges [per meter]'; titl = 'Lunges per Meter';
            FID = plot_bar_msem(lunges_per_m,utest,Xl,Gl,params.range.stat.occ.lunge/2,ytit,titl,FID,gen_ident,params);
        end
        if params.plots.stat.tussl && numel(find(params.plots.stat.tussl_vec == 8)),
            ytit = 'tussling [per meter]'; titl = 'Tusslings per Meter';
            FID = plot_bar_msem(tussls_per_m,utest,Xt,Gt,5,ytit,titl,FID,gen_ident,params);
        end
        if  params.plots.stat.lunging && numel(find(params.plots.stat.lunge_vec == 9)),
            ytit = 'lunges [per minute]'; titl = 'Lunges per Minute';
            FID = plot_bar_msem(lunges_per_min,utest,Xlt,Glt,10,ytit,titl,FID,gen_ident,params);
        end
        if params.plots.stat.tussl && numel(find(params.plots.stat.tussl_vec == 9)),
            ytit = 'tussling [per minute]'; titl = 'Tusslings per Minute';
            FID = plot_bar_msem(tussls_per_min,utest,Xtt,Gtt,2,ytit,titl,FID,gen_ident,params);
        end
    end

    %% SCATTERPLOTS LUNGES <-> TUSSLING <-> CHASING
    if params.plots.stat.scat,
        tussl.mea = []; lunge.mea = []; chase.mea = [];
        tussl.sem = []; lunge.sem = []; chase.sem = [];
        wing_thr.mea = []; wing_thr.sem = [];
        wing_ext.lr.mea = zeros(1,ngens); wing_ext.b.mea = zeros(1,ngens); wing_fli = wing_ext;
        % Collect mean & sem number of actions per genotype
        for igen=1:ngens,
            tussl.mea(igen) = mean(tussls{igen}.number);
            lunge.mea(igen) = mean(lunges{igen}.obj1.number+lunges{igen}.obj2.number);
            chase.mea(igen) = mean(chases{igen}.obj1.number+chases{igen}.obj2.number);
            tussl.sem(igen) = std(tussls{igen}.number) / sqrt(nmovs(igen));
            lunge.sem(igen) = std(lunges{igen}.obj1.number+lunges{igen}.obj2.number) / sqrt(nmovs(igen));
            chase.sem(igen) = std(chases{igen}.obj1.number+chases{igen}.obj2.number) / sqrt(nmovs(igen));
            wing_thr.mea(igen) = mean(wings.threat{igen}.obj1.number+wings.threat{igen}.obj2.number);
            wing_thr.sem(igen) = std(wings.threat{igen}.obj1.number+wings.threat{igen}.obj2.number) / sqrt(nmovs(igen));
            for i=1:nmovs(igen),
                if lunges{igen}.obj1.number(i) > lunges{igen}.obj2.number(i),
                    wing_ext.lr.mea(igen) = wing_ext.lr.mea(igen) + wings.ext.l{igen}.obj1.number(i) + wings.ext.r{igen}.obj1.number(i);
                    wing_ext.b.mea(igen)  = wing_ext.b.mea(igen) + wings.ext.b{igen}.obj1.number(i);
                    wing_fli.lr.mea(igen) = wing_fli.lr.mea(igen) + wings.fli.l{igen}.obj1.number(i) + wings.fli.r{igen}.obj1.number(i);
                    wing_fli.b.mea(igen)  = wing_fli.b.mea(igen) + wings.fli.b{igen}.obj1.number(i);
                else
                    wing_ext.lr.mea(igen) = wing_ext.lr.mea(igen) + wings.ext.l{igen}.obj2.number(i) + wings.ext.r{igen}.obj2.number(i);
                    wing_ext.b.mea(igen)  = wing_ext.b.mea(igen) + wings.ext.b{igen}.obj2.number(i);
                    wing_fli.lr.mea(igen) = wing_fli.lr.mea(igen) + wings.fli.l{igen}.obj2.number(i) + wings.fli.r{igen}.obj2.number(i);
                    wing_fli.b.mea(igen)  = wing_fli.b.mea(igen) + wings.fli.b{igen}.obj2.number(i);
                end
            end
            wing_ext.lr.mea(igen) = wing_ext.lr.mea(igen) / nmovs(igen);
            wing_ext.b.mea(igen) = wing_ext.b.mea(igen) / nmovs(igen);
            wing_fli.lr.mea(igen) = wing_fli.lr.mea(igen) / nmovs(igen);
            wing_fli.b.mea(igen) = wing_fli.b.mea(igen) / nmovs(igen);
            txt(igen) = gen_ident(igen);
        end
        % Scatterplots
        if numel(find(params.plots.stat.scat_vec == 1)),
            FID = scatterplot(lunge,tussl,txt,'lunges','tussling',' ',' ',FID,params);
        end
        if numel(find(params.plots.stat.scat_vec == 2)),
            FID = scatterplot(wing_thr,lunge,txt,'wing threats','lunges',' ',' ',FID,params);
        end
        if numel(find(params.plots.stat.scat_vec == 3)),
            FID = scatterplot(wing_thr,tussl,txt,'wing threats','tussling',' ',' ',FID,params);
        end
        if numel(find(params.plots.stat.scat_vec == 4)),
            FID = scatterplot(chase,lunge,txt,'chases','lunging',' ',' ',FID,params);
        end
        if numel(find(params.plots.stat.scat_vec == 5)),
            FID = scatterplot(chase,tussl,txt,'chases','tussling',' ',' ',FID,params);
        end
        if numel(find(params.plots.stat.scat_vec == 6)),
            FID = scatterplot(wing_ext.lr,lunge,txt,'one wing ext.','lunges',' ',' ',FID,params);
        end
    end
end

end