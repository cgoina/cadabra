%% Plot_heatmaps
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
%% Plot 2d-histograms, the so-called 'heatmaps'

% PLOT HEATMAPS
function [FID,mea_img,sem_img] = plot_heatmaps(field,plot_oppon,titl,ymax,min_bout,max_gap,chamber,gen_ident,dxy,dts,nframes,FID,params,bool_plot)

% PREPARE DATA
if nargin < 14, bool_plot = 1; end
if ~isfield(dxy,'x'), dxy.x = 0; end
if ~isfield(dxy,'y'), dxy.y = dxy.x; end

minx = -chamber.width/2+dxy.x; maxx = chamber.width/2-dxy.x; 
miny = -chamber.height/2+dxy.y; maxy = chamber.height/2-dxy.y; dv = (maxx-minx)/30;
xv = minx:dv:maxx; yv = miny:dv:maxy; ngens = length(field); 
xaxis = minx:10:maxx;
xtick_lab = cell(1,length(xaxis)); for i=1:length(xaxis), xtick_lab{i} = num2str(-chamber.width/2+dxy.x+(i-1)*10); end
yaxis = miny:10:maxy;
ytick_lab = cell(1,length(yaxis)); for i=1:length(yaxis), ytick_lab{i} = num2str(-chamber.height/2+dxy.y+(i-1)*10); end

mea_img = zeros(numel(yv),numel(xv),ngens); sem_img = mea_img;
for igen=1:ngens,
    if params.tim, fac = sum(nframes{igen}.*dts{igen})/sum(nframes{igen}); else fac = 1; end
    nmov = length(field{igen}.obj1.number); img = zeros(numel(xv),numel(yv),nmov);
    cumind1 = [0 cumsum(field{igen}.obj1.number)]; cumindl1 = [1 field{igen}.obj1.len];
    cumind2 = [0 cumsum(field{igen}.obj2.number)]; cumindl2 = [1 field{igen}.obj2.len];
    for j=1:nmov,
        % COMPUTE A 2D HISTOGRAM OF FLY POSITIONS WHILE 
        % PERFORMING AN ACTION FOR EACH MOVIE
        % Fly 1
        if sum(field{igen}.obj1.number(j)) > 0,
            st = sum(cumindl1(1:cumind1(j)+1)); en = sum(cumindl1(2:cumind1(j+1)+1));
            data = [field{igen}.obj1.x1(st:en)' , field{igen}.obj1.y1(st:en)'];
            % Add opponent?
            if plot_oppon,
                data = [data ; [field{igen}.obj1.x2(st:en)' , field{igen}.obj1.y2(st:en)']];
            end
            img(:,:,j) = hist3(data,{minx:dv:maxx miny:dv:maxy}).*fac;
        end
        % Fly 2
        if sum(field{igen}.obj2.number(j)) > 0,
            st = sum(cumindl2(1:cumind2(j)+1)); en = sum(cumindl2(2:cumind2(j+1)+1));
            data = [field{igen}.obj2.x1(st:en)', field{igen}.obj2.y1(st:en)'];
            if plot_oppon,
                data = [data ; [field{igen}.obj2.x2(st:en)', field{igen}.obj2.y2(st:en)']];
            end
            img(:,:,j) = img(:,:,j) + hist3(data,{minx:dv:maxx miny:dv:maxy})*fac;
        end
    end
    % COMPUTE A MEAN AND SEM (NOT PLOTTED) 2D HISTOGRAM
    % OF FLY POSITIONS PER GENOTYPE WHILE PERFORMING AN ACTION
    mea_img(:,:,igen) = mean(img,3)'; sem_img(:,:,igen) = std(img,[],3)'/sqrt(nmov);
end

% Autoscale?
if params.autoscale, 
    ymax = ceil(max(max(max(mea_img)))); 
    ymax = ymax+ceil(.1*abs(ymax)); % extend the display range by 10%
end

% PLOT HEATMAPS
if bool_plot,
    if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
    % MEAN
    figure(FID); clf;
    set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
    for igen=1:ngens,
        h = subplot('Position',[mod(igen-1,params.k(2))/(params.k(2)*1.1)+.05, 0.9-fix((igen-1)/params.k(2))/(params.k(1)*1.1)-1/(params.k(1)*1.5), ...
            1/(params.k(2)*1.5), 1/(params.k(1)*1.5)]); hold on;
        xlim([-chamber.width/2+dxy.x chamber.width/2-dxy.x]);
        ylim([-chamber.height/2+dxy.y chamber.height/2-dxy.y]);
        set(h,'FontSize',params.axisfontsize);
        title(cell2mat(gen_ident(igen)),'FontSize',params.axisfontsize); axis image;
        imagesc(xv,yv,mea_img(:,:,igen),[0 ymax]); drawnow;
        colormap('default'); colormap([0 0 0 ; jet]);
        if ~mod(igen-1,params.k(2)), ylabel('distance [mm]'); end
        if (igen > (params.k(1)-1)*params.k(2)), xlabel('distance [mm]'); end
    end
    pos = get(gca,'Position');
    ch = colorbar('location','EastOutside');
    set(ch,'FontSize',params.axisfontsize);
    if params.tim,
        set(get(ch,'YLabel'),'String','cumulative time [s]');
    else
        set(get(ch,'YLabel'),'String','occurence');
    end
    set(ch,'Position',[0.92 pos(2) 0.02 pos(4)]);
    set(FID,'Name',titl);
    if min_bout >= 0,
        mtit([titl ' (MEAN), b > ' num2str(min_bout*1000,'%4.0f') ' ms, g < ' ...
            num2str(max_gap*1000,'%4.0f') ' ms'],'FontSize',params.axisfontsize+2,'yoff',.04);
    else
        mtit([titl ' (MEAN)'],'FontSize',params.axisfontsize+2,'yoff',.04);
    end
    if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end

    % SEM
    % figure(FID+1); clf;
    % set(FID+1,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
    % for igen=1:ngens,
    %     h = subplot('Position',[mod(igen-1,params.k(2))/(params.k(2)*1.1)+.05, 0.9-fix((igen-1)/params.k(2))/(params.k(1)*1.1)-1/(params.k(1)*1.5), ...
    %                         1/(params.k(2)*1.5), 1/(params.k(1)*1.5)]); hold on;
    %     xlim([-chamber.width/2+dxy.x chamber.width/2-dxy.x]);
    %     ylim([-chamber.height/2+dxy.y chamber.height/2-dxy.y]);
    %     set(h,'FontSize',params.axisfontsize);
    %     title(cell2mat(gen_ident(igen)),'FontSize',params.axisfontsize); axis image;
    %     imagesc(xv,yv,sem_img(:,:,igen),[0 ymax]); drawnow;
    %     set(gca,'XTick',xaxis); set(gca,'XTickLabel',xtick_lab);
    %     set(gca,'YTick',yaxis); set(gca,'YTickLabel',ytick_lab); hold on;
    %     colormap([0 0 0 ; jet]);
    % end
    % pos = get(gca,'Position');
    % ch = colorbar('location','EastOutside');
    % set(ch,'Position',[0.92 pos(2) 0.02 pos(4)]);
    % set(FID+1,'Name',titl);
    % if min_bout >= 0,
    %     mtit([titl ' (SEM), b > ' num2str(min_bout*1000,'%4.0f') ' ms, g < ' ...
    %           num2str(max_gap*1000,'%4.0f') ' ms'],'FontSize',params.axisfontsize+2,'yoff',.04);
    % else
    %     mtit([titl ' (SEM)'],'FontSize',params.axisfontsize+2,'yoff',.04);
    % end
    % if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
    % if params.pdf, print(['-f' num2str(FID+1)],'-dpsc2',appnd,params.PSFileN); end
end

end