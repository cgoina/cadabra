%% Plot_wing_heatm_raster
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
%% Plot heatmaps and raster plots for wing extensions

function FID = plot_wing_heatm_raster(lunges,wings,dts,nframes,dxy,gen_ident,params,chamber,FID)
% WING EXTENSIONS
if params.tim,
    max_val = params.range.heatmaps.tim.wing;
else
    max_val = params.range.heatmaps.occ.wing;
end
plot_oppon = 0;
min_bout = wings.min_bout.long; max_gap = wings.max_gap.long;
ngens = numel(lunges);

% LEFT WING
if params.plots.heatmaps.ext_leftwing,
    tmp = wings.ext.l;
    % Heatmaps
    if numel(find(params.plots.heatmaps.wing_vec == 1)),
        titl = 'Histograms of Fly Positions - Left Wing Extended';
        [FID,mea_l,sem_l] = plot_heatmaps(tmp,plot_oppon,titl,max_val,min_bout,max_gap,chamber,gen_ident,dxy,dts,nframes,FID,params,1);
    end
    % Raster plots
    if numel(find(params.plots.heatmaps.wing_vec == 2)),
        titl = 'Distribution - Left Wing Extended';
        [FID,distr_l] = plot_2ddistr(tmp,titl,10,gen_ident,dxy,dts,nframes,FID,params);
    end
end
% RIGHT WING
if params.plots.heatmaps.ext_rightwing,
    tmp = wings.ext.r;
    % Heatmaps
    if numel(find(params.plots.heatmaps.wing_vec == 1)),
        titl = 'Histograms of Fly Positions - Right Wing Extended';
        [FID,mea_r,sem_r] = plot_heatmaps(tmp,plot_oppon,titl,max_val,min_bout,max_gap,chamber,gen_ident,dxy,dts,nframes,FID,params,1);
    end
    % Raster plots
    if numel(find(params.plots.heatmaps.wing_vec == 2)),
        titl = 'Distribution - Right Wing Extended';
        [FID,distr_r] = plot_2ddistr(tmp,titl,10,gen_ident,dxy,dts,nframes,FID,params);
    end
end

if params.plots.heatmaps.ext_diffleftrightwing && numel(find(params.plots.heatmaps.wing_vec == 1)),
    % PLOT HEATMAP DIFFERENCE BETWEEN EXTENSION OF LEFT-RIGHT WING
    minx = -chamber.width/2+dxy.x; maxx = chamber.width/2-dxy.x;
    miny = -chamber.height/2+dxy.y; maxy = chamber.height/2-dxy.y; dv = (maxx-minx)/30;
    xv = minx:dv:maxx; yv = miny:dv:maxy; xaxis = minx:10:maxx; yaxis = miny:10:maxy;
    xtick_lab = cell(1,length(xaxis));for i=1:length(xaxis), xtick_lab{i} = num2str(-chamber.width/2+dxy.x+(i-1)*5); end
    ytick_lab = cell(1,length(yaxis));for i=1:length(yaxis), ytick_lab{i} = num2str(-chamber.height/2+dxy.y+(i-1)*5); end
    % EXTRACT 2D HISTOGRAMS OF FLY POSITIONS WHILE EXTENDING LEFT/RIGHT WING
    [FID,mea_l,sem_l] = plot_heatmaps(wings.ext.l,0,'',max_val,min_bout,max_gap,chamber,gen_ident,dxy,dts,nframes,FID,params,0);
    [FID,mea_r,sem_r] = plot_heatmaps(wings.ext.r,0,'',max_val,min_bout,max_gap,chamber,gen_ident,dxy,dts,nframes,FID,params,0);

    % Autoscale?
    if params.autoscale,
        ymax = ceil(max(max(max(abs(mea_l-mea_r)))));
        ymax = ymax+ceil(.1*abs(ymax)); % extend the display range by 10%
    else
        ymax = max_val/2;
    end
    
    % Plot heatmaps
    if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
    figure(FID); clf;
    set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
    for igen=1:ngens,
        h = subplot('Position',[mod(igen-1,params.k(2))/(params.k(2)*1.1)+.05, 0.9-fix((igen-1)/params.k(2))/(params.k(1)*1.1)-1/(params.k(1)*1.5), ...
            1/(params.k(2)*1.5), 1/(params.k(1)*1.5)]); hold on;
        xlim([-chamber.width/2+dxy.x chamber.width/2-dxy.x]); ylim([-chamber.height/2+dxy.y chamber.height/2-dxy.y]);
        set(h,'FontSize',params.axisfontsize);
        title(cell2mat(gen_ident(igen)),'FontSize',params.axisfontsize); axis image;
        imagesc(xv,yv,mea_l(:,:,igen)-mea_r(:,:,igen),[-ymax ymax]); drawnow;
        %     set(gca,'XTick',xaxis); set(gca,'XTickLabel',xtick_lab);
        %     set(gca,'YTick',yaxis); set(gca,'YTickLabel',ytick_lab); hold on;
        if ~mod(igen-1,params.k(2)), ylabel('distance [mm]'); end
        if (igen > (params.k(1)-1)*params.k(2)), xlabel('distance [mm]'); end
        col = jet; for j=1:4, col(ceil(length(col)/2)+j-2,:) = [0 0 0]; end
        colormap(col);
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
    titl = 'Comparison LEFT-RIGHT Wing Extended';
    set(FID,'Name',titl);
    if min_bout >= 0,
        mtit([titl ' (MEAN), b > ' num2str(min_bout*1000,'%4.0f') ' ms, g < ' ...
            num2str(max_gap*1000,'%4.0f') ' ms'],'FontSize',params.axisfontsize+2,'yoff',.04);
    else
        mtit([titl ' (MEAN)'],'FontSize',params.axisfontsize+2,'yoff',.04);
    end
    if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
end

if params.plots.heatmaps.ext_onewing,
    % PLOT HEATMAPS AND RASTERPLOTS FOR SINGLE WING EXTENSION (left or right)
    minx = -chamber.width/2+dxy.x; maxx = chamber.width/2-dxy.x;
    miny = -chamber.height/2+dxy.y; maxy = chamber.height/2-dxy.y;
    dv = (maxx-minx)/30;
    xv = minx:dv:maxx; yv = miny:dv:maxy; xaxis = minx:5:maxx; yaxis = miny:5:maxy;
    xtick_lab = cell(1,length(xaxis));for i=1:length(xaxis), xtick_lab{i} = num2str(-chamber.width/2+dxy.x+(i-1)*5); end
    ytick_lab = cell(1,length(yaxis));for i=1:length(yaxis), ytick_lab{i} = num2str(-chamber.height/2+dxy.y+(i-1)*5); end
    % EXTRACT 2D HISTOGRAMS OF FLY POSITIONS AND RASTERPLOTS 
    % WHILE EXTENDING LEFT/RIGHT WING
    [FID,mea_l,sem_l] = plot_heatmaps(wings.ext.l,0,'',max_val,min_bout,max_gap,chamber,gen_ident,dxy,dts,nframes,FID,params,0);
    [FID,mea_r,sem_r] = plot_heatmaps(wings.ext.r,0,'',max_val,min_bout,max_gap,chamber,gen_ident,dxy,dts,nframes,FID,params,0);
    [FID,distr_l,mov_l] = plot_2ddistr(wings.ext.l,'',10,gen_ident,dxy,dts,nframes,FID,params,0);
    [FID,distr_r,mov_r] = plot_2ddistr(wings.ext.r,'',10,gen_ident,dxy,dts,nframes,FID,params,0);

    % Autoscale?
    if params.autoscale,
        ymax = ceil(max(max(max(abs(mea_l-mea_r)))));
        ymax = ymax+ceil(.1*abs(ymax)); % extend the display range by 10%
    else
        ymax = max_val;
    end

    if numel(find(params.plots.heatmaps.wing_vec == 1)),
        % PLOT HEATMAPS (SUM OF LEFT AND RIGHT WING)
        if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
        figure(FID); clf;
        set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
        for igen=1:ngens,
            h = subplot('Position',[mod(igen-1,params.k(2))/(params.k(2)*1.1)+.05, 0.9-fix((igen-1)/params.k(2))/(params.k(1)*1.1)-1/(params.k(1)*1.5), ...
                1/(params.k(2)*1.5), 1/(params.k(1)*1.5)]); hold on;
            xlim([-chamber.width/2+dxy.x chamber.width/2-dxy.x]); ylim([-chamber.height/2+dxy.y chamber.height/2-dxy.y]);
            set(h,'FontSize',params.axisfontsize);
            title(cell2mat(gen_ident(igen)),'FontSize',params.axisfontsize); axis image;
            imagesc(xv,yv,mea_l(:,:,igen)+mea_r(:,:,igen),[0 ymax]); drawnow;
            %     set(gca,'XTick',xaxis); set(gca,'XTickLabel',xtick_lab);
            %     set(gca,'YTick',yaxis); set(gca,'YTickLabel',ytick_lab); hold on;
            if ~mod(igen-1,params.k(2)), ylabel('distance [mm]'); end
            if (igen > (params.k(1)-1)*params.k(2)), xlabel('distance [mm]'); end
            col = jet; for j=1:2, col(j,:) = [0 0 0]; end
            colormap(col);
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
        titl = 'ONE Wing Extended';
        set(FID,'Name',titl);
        if min_bout >= 0,
            mtit([titl ' (MEAN), b > ' num2str(min_bout*1000,'%4.0f') ' ms, g < ' ...
                num2str(max_gap*1000,'%4.0f') ' ms'],'FontSize',params.axisfontsize+2,'yoff',.04);
        else
            mtit([titl ' (MEAN)'],'FontSize',params.axisfontsize+2,'yoff',.04);
        end
        if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
    end

    if numel(find(params.plots.heatmaps.wing_vec == 2)),
        % RASTERPLOTS - ONE WING EXTENDED
        stpx = 5; 
        xaxis = 0:stpx:params.max_frames/60;
        xtick_lab = cell(1,length(xaxis)); 
        for i=1:length(xaxis), xtick_lab{i} = num2str((i-1)*stpx); end
        COL =  .8*[[.87 .76 .54]*.6 ; .87 .76 .54]; % color of fly 1
        COL2 = .8*[[.71 .82 .49]*.6 ; .71 .82 .49]; % color of fly 2
        if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
        figure(FID); clf;
        set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
        miny = 0; maxy = 1; yaxis = miny:.1:maxy;
        ytick_lab = cell(1,length(yaxis));
        for i=1:length(yaxis),
            if i < 5,
                ytick_lab{i} = num2str(yaxis(i)*10);
            else
                ytick_lab{i} = ' ';
            end
        end
        %     if params.k(1) == 1, params.k(1) = 2; end
        for igen=1:ngens,
            h = subplot('Position',[mod(igen-1,params.k(2))/(params.k(2)*1.1)+.05, 0.9-fix((igen-1)/params.k(2))/(params.k(1)*1.1)-1/(params.k(1)*1.5), ...
                1/(params.k(2)*1.5), 1/(params.k(1)*1.5)]);
            set(h,'FontSize',params.axisfontsize); hold on;
            % RASTERPLOT OF WING EXTENSIONS
            nmov0 = length(lunges{igen}.obj1.number); nmov = nmov0 * 1.5;
            dw = 1/nmov; db = .8*dw/2;
            % LEFT WING
            for i=1:numel(distr_l{igen}),
                % Flies get different colors
                if mod(i,2),
                    if mod(mov_l{igen}(i),2), c = COL(1,:); else c = COL2(1,:); end
                else
                    if mod(mov_l{igen}(i),2), c = COL(2,:); else c = COL2(2,:); end
                end
                plot([distr_l{igen}(i) distr_l{igen}(i)],dw/2+[(mov_l{igen}(i)+nmov-nmov0)/nmov-dw-db (mov_l{igen}(i)+nmov-nmov0)/nmov-dw+db],'Color',c); hold on;
            end
            % RIGHT WING
            for i=1:numel(distr_r{igen}),
                % Flies get different colors
                if mod(i,2),
                    if mod(mov_r{igen}(i),2), c = COL(1,:); else c = COL2(1,:); end
                else
                    if mod(mov_r{igen}(i),2), c = COL(2,:); else c = COL2(2,:); end
                end
                plot([distr_r{igen}(i) distr_r{igen}(i)],dw/2+[(mov_r{igen}(i)+nmov-nmov0)/nmov-dw-db (mov_r{igen}(i)+nmov-nmov0)/nmov-dw+db],'Color',c); hold on;
            end            
            % HISTOGRAM - NUMBER OF WING EXTENSIONS PER TIME BIN
            h = hist([distr_l{igen} distr_r{igen}],0:1:30); h = h/nmov0/10;
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
        titl = 'ONE Wing Extended';
        set(FID,'Name',titl);
        mtit(titl,'FontSize',params.axisfontsize+2,'yoff',.04);
    if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
    end
end

% WING FLICK HEATMAPS AND RASTERPLOTS (future task)
% DATA PREPARED BUT NOT PLOTTED, BECAUSE NOT VERIFIED

