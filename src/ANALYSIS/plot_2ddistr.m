%% Plot_2ddistr
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
%% Raster plots for the given genotypes

% COMBINATION OF RASTERPLOTS AND HISTOGRAM
function [FID,distr,mov] = plot_2ddistr(field,titl,ymax,gen_ident,dxy,dts,nframes,FID,params,bool_plot)

if nargin < 10, bool_plot = 1; end
if ~isfield(dxy,'x'), dxy.x = 0; end
if ~isfield(dxy,'y'), dxy.y = dxy.x; end

ngens = length(field); stpx = 5; 
xaxis = 0:stpx:params.max_frames/60;
xtick_lab = cell(1,length(xaxis)); 
for i=1:length(xaxis), xtick_lab{i} = num2str((i-1)*stpx); end

% PREPARE DATA
distr = cell(1,ngens); mov = distr; maxmov = 0;
for igen=1:ngens,
    if params.tim, dt = sum(nframes{igen}.*dts{igen})/sum(nframes{igen}); else dt = 1; end
    nmov = length(field{igen}.obj1.number); if nmov>maxmov, maxmov = nmov; end
    cumind1 = [0 cumsum(field{igen}.obj1.number)]; cumindl1 = [1 field{igen}.obj1.len];
    cumind2 = [0 cumsum(field{igen}.obj2.number)]; cumindl2 = [1 field{igen}.obj2.len];
    for j=1:nmov,
%         if ~params.courtship,
            if sum(field{igen}.obj1.number(j)) > 0,
                % Prepare time stamps (distr) and movie/fly pair numbers (mov) 
                if sum(cumindl1)<=length(field{igen}.obj1.t),
                    st = sum(cumindl1(1:cumind1(j)+1)); en = sum(cumindl1(2:cumind1(j+1)+1));
                else
                    st = cumind1(j)+1; en = cumind1(j+1);
                end
                distr{igen} = [distr{igen} field{igen}.obj1.t(st:en)/60];
                mov{igen} = [mov{igen} ones(1,numel(st:en))*j];
            end
%         end
        if sum(field{igen}.obj2.number(j)) > 0,
            if sum(cumindl1)<=length(field{igen}.obj1.t),
                st = sum(cumindl2(1:cumind2(j)+1)); en = sum(cumindl2(2:cumind2(j+1)+1));
            else
                st = cumind2(j)+1; en = cumind2(j+1);
            end
            distr{igen} = [distr{igen} field{igen}.obj2.t(st:en)/60];
            mov{igen} = [mov{igen} ones(1,numel(st:en))*j];
        end
    end
end

% PLOTS
if bool_plot,
    COL =  .8*[[.87 .76 .54]*.6 ; .87 .76 .54]; COL2 = .8*[[.71 .82 .49]*.6 ; .71 .82 .49];
    if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end 
    figure(FID); clf;
    set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
    miny = 0; maxy = 1; yaxis = miny:.1:maxy;
    ytick_lab = cell(1,length(yaxis)); hmax = 0;
    for igen=1:ngens,
        nmov0 = length(field{igen}.obj1.number);
        h = hist(distr{igen},0:1:params.max_frames/60); h = h/nmov0;
        if hmax < max(h), hmax = max(h); end
    end
    fac = 1/ceil(hmax*3);
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
        nmov0 = length(field{igen}.obj1.number); nmov = nmov0 * 1.5;
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
        h = hist(distr{igen},0.5:1:params.max_frames/60); h = h/nmov0*fac;
        dbar = .2;
        for ii=1:max(xaxis),
            patch([ii-dbar ii+dbar ii+dbar ii-dbar ii-dbar]-.5,...
                  [0 0 h(ii) h(ii) 0],'k');
        end
        set(gca,'YTick',yaxis); set(gca,'YTickLabel',ytick_lab); hold on;
        set(gca,'XTick',xaxis); set(gca,'XTickLabel',xtick_lab);
        colormap('default'); colormap([0 0 0 ; jet]);
        if ~mod(igen-1,params.k(2)), ylabel('occurrence  |  fly pair -->'); end
        if (igen > (params.k(1)-1)*params.k(2)), xlabel('time [min]'); end
        title(cell2mat(gen_ident(igen)),'FontSize',params.axisfontsize);
        axis([xaxis(1) xaxis(end) yaxis(1) yaxis(end)]);
    end
    set(FID,'Name',titl);
    mtit(titl,'FontSize',params.axisfontsize+2,'yoff',.04);
    if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
end
%%