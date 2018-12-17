%% Plot_bar_msem
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
%% Plots statistics either as boxplot or barplot with mean and standard
%% error of the mean (SEM)

% BOXPLOT / BARPLOT ROUTINE
function FID = plot_bar_msem(data,utest,X,G,maxy,ytit,titl,FID,gen_ident,params)

ngens = size(data.mea,1);
x = 1:ngens;
maxx = [x(1)-.5 x(end)+.5];
% Autoscale?
if params.autoscale,
    maxy = [floor(min(X)) ceil(max(X))];
    maxy = [maxy(1)-ceil(.1*abs(maxy(1))) maxy(2)+ceil(.1*abs(maxy(2)))];
end
if (length(maxy) == 1),
    maxy = [0 maxy];
end
maxy = double(maxy);
if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
figure(FID); clf;
set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
if params.tim && ~numel(strfind(ytit,']')),
    maxt = maxy(2);
    yaxis = 0:60:maxt;
    ytick_lab = cell(1,length(yaxis)); for i=1:length(yaxis), ytick_lab{i} = num2str(yaxis(i)/60); end
    set(gca,'YTick',yaxis); set(gca,'YTickLabel',ytick_lab);
    ytit = [ytit ' [min]'];
end
if ~params.boxplot || ~numel(X),
    % BARPLOT
    % Plot bars (patches)
    patch_plot(x,data.mea,params);
    if (size(data.mea,2) > 1) && ~params.oneobj,
        d1 = -params.bardispl; d2 = -d1;
        if (size(data.mea,2) > 2), 
            w = params.semwidth*2/3;
        else
            w = params.semwidth;
        end
    else
        d1 = 0;
        w = params.semwidth;
    end
    % Add whiskers (standard error of the mean [sem])
    for i=1:ngens,
        plot([x(i)+d1 x(i)+d1],[data.mea(i,1)-data.sem(i,1) data.mea(i,1)+data.sem(i,1)],'k','Linewidth',1);
        plot([x(i)+d1-w x(i)+d1+w],[data.mea(i,1)-data.sem(i,1) data.mea(i,1)-data.sem(i,1)],'k','Linewidth',1);
        plot([x(i)+d1-w x(i)+d1+w],[data.mea(i,1)+data.sem(i,1) data.mea(i,1)+data.sem(i,1)],'k','Linewidth',1);
        if (size(data.mea,2) == 2) && ~params.oneobj,
            plot([x(i)+d2 x(i)+d2],[data.mea(i,2)-data.sem(i,2) data.mea(i,2)+data.sem(i,2)],'k');
            plot([x(i)+d2-w x(i)+d2+w],[data.mea(i,2)-data.sem(i,2) data.mea(i,2)-data.sem(i,2)],'k');
            plot([x(i)+d2-w x(i)+d2+w],[data.mea(i,2)+data.sem(i,2) data.mea(i,2)+data.sem(i,2)],'k');
        elseif (size(data.mea,2) == 3) && ~params.oneobj,
            plot([x(i) x(i)],[data.mea(i,2)-data.sem(i,2) data.mea(i,2)+data.sem(i,2)],'k');
            plot([x(i)-w x(i)+w],[data.mea(i,2)-data.sem(i,2) data.mea(i,2)-data.sem(i,2)],'k');
            plot([x(i)-w x(i)+w],[data.mea(i,2)+data.sem(i,2) data.mea(i,2)+data.sem(i,2)],'k');
            plot([x(i)+d2 x(i)+d2],[data.mea(i,3)-data.sem(i,3) data.mea(i,3)+data.sem(i,3)],'k');
            plot([x(i)+d2-w x(i)+d2+w],[data.mea(i,3)-data.sem(i,3) data.mea(i,3)-data.sem(i,3)],'k');
            plot([x(i)+d2-w x(i)+d2+w],[data.mea(i,3)+data.sem(i,3) data.mea(i,3)+data.sem(i,3)],'k');
        end
    end
    % Add minima/maxima in case the information is provided
    if isfield(data,'min'),
        if (size(data.mea,2) == 1) || params.oneobj,
            plot(x,data.min(:,1),'+','MarkerEdgeColor','k','MarkerSize',4,'LineWidth',1); hold on;
            plot(x,data.max(:,1),'+','MarkerEdgeColor','k','MarkerSize',4,'LineWidth',1); hold on;
        elseif (size(data.mea,2) == 2) && ~params.oneobj,
            plot(x+d1,data.min(:,1),'+','MarkerEdgeColor','k','MarkerSize',4,'LineWidth',1); hold on;
            plot(x+d1,data.max(:,1),'+','MarkerEdgeColor','k','MarkerSize',4,'LineWidth',1); hold on;
            plot(x+d2,data.min(:,2),'+','MarkerEdgeColor','k','MarkerSize',4,'LineWidth',1); hold on;
            plot(x+d2,data.max(:,2),'+','MarkerEdgeColor','k','MarkerSize',4,'LineWidth',1); hold on;
        elseif (size(data.mea,2) == 3) && ~params.oneobj,
            plot(x+d1,data.min(:,1),'+','MarkerEdgeColor','k','MarkerSize',4,'LineWidth',1); hold on;
            plot(x+d1,data.max(:,1),'+','MarkerEdgeColor','k','MarkerSize',4,'LineWidth',1); hold on;
            plot(x,data.min(:,2),'+','MarkerEdgeColor','k','MarkerSize',4,'LineWidth',1); hold on;
            plot(x,data.max(:,2),'+','MarkerEdgeColor','k','MarkerSize',4,'LineWidth',1); hold on;
            plot(x+d2,data.min(:,3),'+','MarkerEdgeColor','k','MarkerSize',4,'LineWidth',1); hold on;
            plot(x+d2,data.max(:,3),'+','MarkerEdgeColor','k','MarkerSize',4,'LineWidth',1); hold on;
        end
    end
    % Legend
    if ~params.oneobj && (size(data.mea,2) == 2),
        if params.winlos, h=legend('winner','loser'); else h=legend('fly 1','fly 2'); end
        legend(h,'boxoff');
    elseif ~params.oneobj && (size(data.mea,2) == 3),
        text(ngens,maxy(2),'< 5 mm, 5 - 10 mm, > 10 mm','HorizontalAlignment','Right','VerticalAlignment','Top');
    end
    if isstruct(utest) && numel(X), plot_utests(utest,X+min(X),G,data,max(maxy),params); end
    ylim([maxy(1) maxy(2)]); ylabel(ytit);
    axis([maxx(1) maxx(2) maxy(1) maxy(2)]);
    set(gca,'XTick',x,'FontSize',params.axisfontsize);
    set_xtick_label(gen_ident,45,params.cross_to,params.axisfontsize);
else
    % BOXPLOT
    % Plot boxes, whiskers, everything
    boxplot(X,G,'whisker',params.whisker,'plotstyle','traditional','notch','marker'); hold on;
    % Legend
    if ~params.oneobj && (size(data.mea,2) == 2),
        gen_ident2 = cell(1,2*ngens); icnt = 0;
        for i=2:2:2*ngens,
            icnt = icnt + 1;
            gen_ident2{i} = gen_ident{icnt};
        end
        ylim([maxy(1) maxy(2)]);
        axis([maxx(1) 2*ngens+1 maxy(1) maxy(2)]);
        set_xtick_label(gen_ident2,45,params.cross_to,params.axisfontsize);
        set(gca,'XTick',1.5:2:2*ngens-.5); 
        set(gca,'XTickLabel',' ');
        if params.winlos,
            text(2*ngens,maxy(2),'L: winner, R: loser','HorizontalAlignment','Right','VerticalAlignment','Top');
        else
            text(2*ngens,maxy(2),'L: fly 1, R: fly 2','HorizontalAlignment','Right','VerticalAlignment','Top');
        end
    elseif ~params.oneobj && (size(data.mea,2) == 3),
        gen_ident2 = cell(1,3*ngens); icnt = 0;
        for i=2:3:3*ngens,
            icnt = icnt + 1;
            gen_ident2{i} = gen_ident{icnt};
        end
        ylim([maxy(1) maxy(2)]);
        axis([maxx(1) 3*ngens+1 maxy(1) maxy(2)]);
        set_xtick_label(gen_ident2,45,params.cross_to,params.axisfontsize);
        set(gca,'XTick',1.5:3:3*ngens-.5);
        set(gca,'XTickLabel',' ');
        text(2*ngens,maxy(2),'< 5 mm, 5 - 10 mm, > 10 mm','HorizontalAlignment','Right','VerticalAlignment','Top');
    else
        if isstruct(utest) && numel(X), plot_utests(utest,X+min(X),G,data,max(maxy),params); end
        ylim([maxy(1) maxy(2)]);
        axis([maxx(1) maxx(2) maxy(1) maxy(2)]);
        set(gca,'XTick',x,'FontSize',params.axisfontsize);
        set_xtick_label(gen_ident,45,params.cross_to,params.axisfontsize);
        set(gca,'XTickLabel',' ');
    end
end
ylabel(ytit);
set(gca,'FontSize',params.axisfontsize);
if params.plots.stat.traveldist && params.courtship, titl = [titl ' (Pre-Copulation)']; end
set(FID,'Name',titl); mtit(titl,'FontSize',params.axisfontsize+2,'yoff',.04);
if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end

end
