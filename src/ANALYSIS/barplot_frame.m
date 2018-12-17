%% Barplot_frame
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
%% Plotting routine for 'analyze_copulation'

function FID = barplot_frame(frdata,dts,titl,maxt,gen_ident,mov_count,cross_to,FID,params)
% PLOT BARS WITH MEAN, SEM, MIN, MAX
ngens = length(frdata);
x = 1:ngens; y.mea = zeros(1,ngens); y.sem = y.mea; 
y.min = y.mea; y.max = y.mea; y.n = y.mea;
for igen=1:ngens,
    data = frdata{igen}.*dts{igen};
    data = data(data>0);
    if numel(data)>0,
        y.n(igen) = numel(data);
        y.mea(igen) = mean(data);
        y.sem(igen) = std(data) / sqrt(numel(data));
        y.min(igen) = min(data); y.max(igen) = max(data);
    end
end
if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
figure(FID); clf;
set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
bar(x,y.mea,.5); colormap([.87 .76 .54]); hold on;
d1 = 0; d2 = -d1; w = .2;
for i=1:length(x),
    plot([x(i)+d1 x(i)+d1],[y.mea(i)-y.sem(i) y.mea(i)+y.sem(i)],'k','Linewidth',2);
    plot([x(i)+d1-w x(i)+d1+w],[y.mea(i)-y.sem(i) y.mea(i)-y.sem(i)],'k','Linewidth',2);
    plot([x(i)+d1-w x(i)+d1+w],[y.mea(i)+y.sem(i) y.mea(i)+y.sem(i)],'k','Linewidth',2);
end
plot(x,y.min,'+','MarkerEdgeColor','r','MarkerSize',4,'LineWidth',2); hold on;
plot(x,y.max,'+','MarkerEdgeColor','r','MarkerSize',4,'LineWidth',2); hold on;
set(gca,'FontSize',params.axisfontsize);
set(gca,'XTick',1:ngens); set_xtick_label(gen_ident,90,cross_to,params.axisfontsize);

maxt = maxt*60; dt = 120; t = 0:dt:maxt; ytick_lab = cell(1,length(t));
for i=1:length(t), ytick_lab{i} = num2str(t(i)/60); end
set(gca,'YTick',0:dt:maxt); set(gca,'YTickLabel',ytick_lab);
ylim([0 maxt]); ylabel([titl ' [min]']);
axis([0 ngens+1 0 maxt]);
for igen=1:ngens,
    text(x(igen)-.2,maxt/30,[num2str(y.n(igen)) '/' num2str(mov_count(igen))],'FontSize',8);
end
set(FID,'Name',titl); mtit(titl,'FontSize',params.axisfontsize,'yoff',.04);
if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
end