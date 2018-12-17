%% Plot_velocity_heatmaps
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
%% Plot velocity heatmaps

function FID = plot_velocity_heatmaps(obj1,obj2,gen_ident,params,chamber,path,FID)

% VELOCITY & DIRECTION HEATMAPS
if params.plots.heatmaps.velocity,
    % PREPARE DATA
    ngens = numel(obj1);
    minx = -chamber.width/2; maxx = chamber.width/2;
    miny = -chamber.height/2; maxy = chamber.height/2;
    dr = 2; minv = -20; maxv = 20; dv = 1;
    HistFileN = [path 'x_' num2str(minx) '-' num2str(maxx) ...
        '_y_' num2str(miny) '-' num2str(maxy) '_dr_' num2str(dr) ...
        '_v_' num2str(minv) '-' num2str(maxv) '_dv_' num2str(dv)];
    fid = fopen([HistFileN '.mat'],'r');
    if (fid < 0),
        xv = minx:dr:maxx; yv = miny:dr:maxy;
        hists.mea.x = zeros(length(xv)-1,length(yv)-1,ngens); hists.mea.y = hists.mea.x;
        hists.std.x = hists.mea.x; hists.std.y = hists.mea.x;
        for igen=1:ngens,
            tic;
            for ix=1:length(xv)-1,
                for iy=1:length(yv)-1,
                    ind1 = find(((obj1{igen}.x >= xv(ix)) & (obj1{igen}.x < xv(ix+1))) & ...
                        ((obj1{igen}.y >= yv(iy)) & (obj1{igen}.y < yv(iy+1))));
                    ind2 = find(((obj2{igen}.x >= xv(ix)) & (obj2{igen}.x < xv(ix+1))) & ...
                        ((obj2{igen}.y >= yv(iy)) & (obj2{igen}.y < yv(iy+1))));
                    if numel(ind1) || numel(ind2),
                        data.x = [obj1{igen}.vx(ind1) obj2{igen}.vx(ind2)];
                        data.y = [obj1{igen}.vy(ind1) obj2{igen}.vy(ind2)];
                        data.r = sqrt(data.x.^2+data.y.^2);
                        % Velocity heatmaps (mean, stdev) data for vx, vy
                        hists.mea.x(ix,iy,igen) = mean(data.x);
                        hists.mea.y(ix,iy,igen) = mean(data.y);
                        hists.std.x(ix,iy,igen) = std(data.x);
                        hists.std.y(ix,iy,igen) = std(data.y);
                    end
                end
            end
            ttt=toc; fprintf(1,'left %2.0f:%2.0f mm:ss\n',floor(ttt * (ngens-igen+1)/60),mod(ttt * (ngens-igen+1),60));
        end
        % Velocity heatmaps (mean, stdev) data for resulting v and direction phi
        hists.mea.rr = sqrt(hists.mea.x.^2+hists.mea.y.^2);
        hists.std.rr = sqrt(hists.std.x.^2+hists.std.y.^2);
        hists.mea.phi= atan2(hists.mea.y,hists.mea.x);
        hists.std.phi= atan2(hists.std.y,hists.std.x);
        save(HistFileN, 'hists');
    else
        load(HistFileN);
    end
    % Autoscale?
    if params.autoscale,
        max_val = ceil(max(max(max(hists.mea.rr))));
        max_val = max_val+ceil(.1*abs(max_val)); % extend the display range by 10%
    else
        % Time spent or occurrence?
        if params.tim,
            max_val = params.range.heatmaps.tim.pos;
        else
            max_val = params.range.heatmaps.occ.pos;
        end
    end
    % PLOT MEAN VELOCITIES
    if numel(find(params.plots.heatmaps.velocity_vec == 1)),
        if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
        figure(FID); clf;
        set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
        col = jet; for j=1:1, col(j,:) = [0 0 0]; end
        xv = minx:dr:maxx; yv = miny:dr:maxy; colormap(col);
        for igen=1:ngens,
            h = subplot('Position',[mod(igen-1,params.k(2))/(params.k(2)*1.1)+.05, ...
                0.9-fix((igen-1)/params.k(2))/(params.k(1)*1.1)-1/(params.k(1)*1.5), ...
                1/(params.k(2)*1.5), 1/(params.k(1)*1.5)]); hold on;
            set(h,'FontSize',params.axisfontsize);
            xlim([minx maxx]); ylim([miny maxy]);
            title(cell2mat(gen_ident(igen)),'FontSize',params.axisfontsize); axis image;
            imagesc(xv,yv,(hists.mea.rr(:,:,igen))',[0 max_val]); drawnow;
            if ~mod(igen-1,params.k(2)), ylabel('distance [mm]'); end
            if (igen > (params.k(1)-1)*params.k(2)), xlabel('distance [mm]'); end
        end
        pos = get(gca,'Position');
        ch = colorbar('location','EastOutside');
        set(get(ch,'YLabel'),'String','velocity [mm/s]');
        set(ch,'FontSize',params.axisfontsize);
        set(ch,'Position',[0.92 pos(2) 0.02 pos(4)]);
        titl = 'Mean Velocities';
        set(FID,'Name',titl); mtit(titl,'FontSize',params.axisfontsize+2,'yoff',.04);
        if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
    end

    if numel(find(params.plots.heatmaps.velocity_vec == 2)),
        % PLOT MEAN VELOCITY DIRECTIONS
        if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
        figure(FID); clf;
        set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
        col = jet; for j=1:2, col(ceil(length(col)/2)+j-1,:) = [0 0 0]; end
        xv = minx:dr:maxx; yv = miny:dr:maxy; colormap(col);
        for igen=1:ngens,
            h = subplot('Position',[mod(igen-1,params.k(2))/(params.k(2)*1.1)+.05, ...
                0.9-fix((igen-1)/params.k(2))/(params.k(1)*1.1)-1/(params.k(1)*1.5), ...
                1/(params.k(2)*1.5), 1/(params.k(1)*1.5)]); hold on;
            set(h,'FontSize',params.axisfontsize);
            xlim([minx maxx]); ylim([miny maxy]);
            title(cell2mat(gen_ident(igen)),'FontSize',params.axisfontsize); axis image;
            imagesc(xv,yv,rad2deg(hists.mea.phi(:,:,igen))',[-180 180]); drawnow;
            if ~mod(igen-1,params.k(2)), ylabel('distance [mm]'); end
            if (igen > (params.k(1)-1)*params.k(2)), xlabel('distance [mm]'); end
        end
        pos = get(gca,'Position');
        ch = colorbar('location','EastOutside');
        set(get(ch,'YLabel'),'String','velocity direction [deg]');
        set(ch,'FontSize',params.axisfontsize);
        set(ch,'Position',[0.92 pos(2) 0.02 pos(4)]);
        titl = 'Mean Directions';
        set(FID,'Name',titl); mtit(titl,'FontSize',params.axisfontsize+2,'yoff',.04);
        if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
    end
end

%% Velocity Histograms
% minv = -20; maxv = 20; dv = 1;
% xv = minv:dv:maxv; yv = minv:dv:maxv;
% figure(2); clf;
% set(2,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
% if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
% for igen=1:ngens,
%     subplot(params.k(1),params.k(2),igen); hold on; xlim([minv maxv]); ylim([minv maxv]); title(cell2mat(gen_ident(igen))); axis square;
%     nframes = numel(proxi.times{igen}.near) * params.max_frames/dts(igen);
%     img = hist3([[obj1{igen}.vy obj2{igen}.vy]',[obj1{igen}.vx obj2{igen}.vx]'],{minv:dv:maxv minv:dv:maxv})/nframes*100;
%     imagesc(xv,yv,img,[0 200/nframes*100]); drawnow;
%     if ~mod(igen-1,params.k(2)), ylabel('distance [mm]'); end
%     if (igen > (params.k(1)-1)*params.k(2)), xlabel('distance [mm]'); end
% end
% pos = get(gca,'Position');
% ch = colorbar('location','EastOutside'); 
% set(ch,'FontSize',params.axisfontsize);
% set(ch,'Position',[0.92 pos(2) 0.02 pos(4)]);
% titl = 'Histograms of Fly Velocities';
% set(2,'Name',titl); mtit(titl,'FontSize',14','yoff',.04);
% if params.pdf, print('-f2','-dpsc2',appnd,params.PSFileN); end
