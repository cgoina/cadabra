%% Analyze_positions
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
%% Analyze fly positions

function FID = analyze_positions(obj1,obj2,dts,nframes,gen_ident,params,chamber,FID)
% FLY POSITION HISTOGRAMS
if params.plots.heatmaps.position,
    % PREPARE DATA
    ngens = numel(obj1);
    minx = -chamber.width/2; maxx = chamber.width/2;
    miny = -chamber.height/2; maxy = chamber.height/2; dv = .5;
    xv = minx:dv:maxx; yv = miny:dv:maxy;
    mea_img = zeros(numel(yv),numel(xv),ngens); sem_img = mea_img;
    for igen=1:ngens,
        % Compute a 2d position histogram for each movie
        if params.tim,
            dt = sum(nframes{igen}.*dts{igen})/sum(nframes{igen});
        else
            dt = 1;
        end
        if dt < 1/15, fac = 2; else fac =1; end
        ind = [0 find(obj1{igen}.t(2:end)-obj1{igen}.t(1:end-1)<0) length(obj1{igen}.t)];
        img = zeros(numel(xv),numel(yv),numel(ind)-1);
        for j=1:numel(ind)-1,
            if ~params.oneobj,
                img(:,:,j) = hist3([[obj1{igen}.x(ind(j)+1:ind(j+1)) obj2{igen}.x(ind(j)+1:ind(j+1))]',...
                    [obj1{igen}.y(ind(j)+1:ind(j+1)) obj2{igen}.y(ind(j)+1:ind(j+1))]'],{minx:dv:maxx miny:dv:maxy})*dt;
            else
                img(:,:,j) = hist3([[obj1{igen}.x(ind(j)+1:ind(j+1))]',...
                    [obj1{igen}.y(ind(j)+1:ind(j+1))]'],{minx:dv:maxx miny:dv:maxy})*dt*fac;
            end
        end
        % Compute mean and sem histogram for each genotype
        mea_img(:,:,igen) = mean(img,3)'; sem_img(:,:,igen) = std(img,[],3)'/sqrt(numel(ind)-1);
    end
    % Autoscale?
    if params.autoscale,
        max_val = ceil(max(max(max(mea_img))));
        max_val = max_val+ceil(.1*abs(max_val)); % extend the display range by 10%
    else
        % Time spent or occurrence?
        if params.tim,
            max_val = params.range.heatmaps.tim.pos;
        else
            max_val = params.range.heatmaps.occ.pos;
        end
    end
    
    % PLOT FLY POSITION HISTOGRAMS - MEAN
    if numel(find(params.plots.heatmaps.position_vec == 1)),
        if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
        figure(FID); clf;
        page.w = 11; page.h = 8.5;
        set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 page.w page.h]);
        for igen=1:ngens,
            h = subplot('Position',[mod(igen-1,params.k(2))/(params.k(2)*1.1)+.05, ...
                0.9-fix((igen-1)/params.k(2))/(params.k(1)*1.1)-1/(params.k(1)*1.5), ...
                1/(params.k(2)*1.5), 1/(params.k(1)*1.5)]); hold on;
            %     subplot(params.k(1),params.k(2),igen); hold on;
            set(h,'FontSize',params.axisfontsize);
            xlim([-chamber.width/2 chamber.width/2]); ylim([-chamber.height/2 chamber.height/2]);
            title(cell2mat(gen_ident(igen)),'FontSize',params.axisfontsize); axis image;
            imagesc(xv,yv,mea_img(:,:,igen),[0 max_val]); drawnow;
            if ~mod(igen-1,params.k(2)), ylabel('distance [mm]'); end
            if (igen > (params.k(1)-1)*params.k(2)), xlabel('distance [mm]'); end
            colormap([0 0 0 ; jet]);
        end
        pos = get(gca,'Position');
        ch = colorbar('location','EastOutside');
        if params.tim,
            set(get(ch,'YLabel'),'String','cumulative time [s]');
        else
            set(get(ch,'YLabel'),'String','occurence');
        end
        set(ch,'FontSize',params.axisfontsize);
        set(ch,'Position',[0.94 pos(2) 0.01 pos(4)]);
        titl = 'Histograms of Fly Positions (MEAN)';
        set(FID,'Name',titl); mtit(titl,'FontSize',params.axisfontsize+2,'yoff',.02);
        if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
    end

    % FLY POSITION HISTOGRAMS - SEM
    if numel(find(params.plots.heatmaps.position_vec == 2)),
        if params.autoscale,
            max_val = ceil(max(max(max(sem_img))));
            max_val = max_val+ceil(.1*abs(max_val)); % extend the display range by 10%
            if ~max_val, max_val = 1; end
        end            
        if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
        figure(FID); clf;
        set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
        for igen=1:ngens,
            h = subplot('Position',[mod(igen-1,params.k(2))/(params.k(2)*1.1)+.05, ...
                0.9-fix((igen-1)/params.k(2))/(params.k(1)*1.1)-1/(params.k(1)*1.5), ...
                1/(params.k(2)*1.5), 1/(params.k(1)*1.5)]); hold on;
            set(h,'FontSize',params.axisfontsize);
            xlim([-chamber.width/2 chamber.width/2]); ylim([-chamber.height/2 chamber.height/2]);
            title(cell2mat(gen_ident(igen)),'FontSize',params.axisfontsize); axis image;
            imagesc(xv,yv,sem_img(:,:,igen),[0 max_val/2]); drawnow;
            if ~mod(igen-1,params.k(2)), ylabel('distance [mm]'); end
            if (igen > (params.k(1)-1)*params.k(2)), xlabel('distance [mm]'); end
            colormap([0 0 0 ; jet]);
        end
        pos = get(gca,'Position');
        ch = colorbar('location','EastOutside');
        if params.tim,
            set(get(ch,'YLabel'),'String','cumulative time [s]');
        else
            set(get(ch,'YLabel'),'String','occurence');
        end
        set(ch,'FontSize',params.axisfontsize);
        set(ch,'Position',[0.93 pos(2) 0.01 pos(4)]);
        titl = 'Histograms of Fly Positions (SEM)';
        set(FID,'Name',titl); mtit(titl,'FontSize',params.axisfontsize+2,'yoff',.04);
        if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
    end
end
