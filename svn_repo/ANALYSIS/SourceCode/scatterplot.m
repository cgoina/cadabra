%% Scatterplot
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
%% Scatterplot for one action of provided genotypes showing mean and
%% standard error of the mean (SEM)

% SCATTERPLOT OF ONE ACTION AGAINST ANOTHER ACTION
function FID = scatterplot(x,y,txt,xtitl,ytitl,xunit,yunit,FID,params)

% x, y are two difference actions; txt is the genotype name
% x, y consist of ".mea" and ".sem", vectors of mean and sem occurrence
% frequency of an action per fly pair/movie
if max(x.mea)>10, maxx = ceil(double(max(x.mea))/10)*10; else maxx = ceil(max(x.mea)); end
if max(y.mea)>10, maxy = ceil(double(max(y.mea))/10)*10; else maxy = ceil(max(y.mea)); end
if (maxx == 0), maxx = 1; end
if (maxy == 0), maxy = 1; end

if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
figure(FID); clf;
set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
dlu = max(x.mea)/40; dtu = 0; %max(y)/40;
semcol = [0.8 0.8 0.8]; dwx = maxx/40; dwy = maxy/40;
% PLOT SEM BARS IN X AND Y DIRECTION
if isfield(x,'sem'),
    for i=1:length(x.mea),
        plot([x.mea(i)+.1-x.sem(i) x.mea(i)+.1-x.sem(i)],[y.mea(i)-dwy y.mea(i)+dwy],'Color',semcol); hold on;
        plot([x.mea(i)+.1+x.sem(i) x.mea(i)+.1+x.sem(i)],[y.mea(i)-dwy y.mea(i)+y.sem(i)*dwy],'Color',semcol);
        plot([x.mea(i)+.1-dwx x.mea(i)+.1+dwx],[y.mea(i)-y.sem(i) y.mea(i)-y.sem(i)],'Color',semcol);
        plot([x.mea(i)+.1-dwx x.mea(i)+.1+dwx],[y.mea(i)+y.sem(i) y.mea(i)+y.sem(i)],'Color',semcol);
        plot([x.mea(i)+.1 x.mea(i)+.1],[y.mea(i)-y.sem(i) y.mea(i)+y.sem(i)],'Color',semcol);
        plot([x.mea(i)+.1-x.sem(i) x.mea(i)+.1+x.sem(i)],[y.mea(i) y.mea(i)],'Color',semcol);
    end
end
% PLOT MEANS AS CIRCLES, ATTACH GENOTYPE NAME
plot([0 maxx],[0 maxy],'--','Color',[.71 .82 .49],'LineWidth',2); hold on;
plot(x.mea+.1,y.mea,'o','MarkerEdgeColor',[.87 .76 .54],'MarkerFaceColor','k','MarkerSize',9,'LineWidth',4);
text(x.mea+dlu,y.mea+dtu,txt,'Color','k','FontSize',params.axisfontsize);
% ADD A LEGEND WITH SOME STATISTICS (CORRELATION, 
% STANDARD DEVIATION IN X, Y DIRECTION)
if numel(y.mea) > 1,
    corr = corrcoef(y.mea,x.mea); corr = corr(1,2);
    text(maxx*.05,maxy*.95,['corr = ' num2str(corr,'%5.2f')],'Color','k','FontSize',params.axisfontsize);
    text(maxx*.05,maxy*.88,['\sigma_{lunges} = ' num2str(std(x.mea),'%5.2f')],'Color','k','FontSize',params.axisfontsize);
    text(maxx*.05,maxy*.81,['\sigma_{tussls} = ' num2str(std(y.mea),'%5.2f')],'Color','k','FontSize',params.axisfontsize);
end
set(gca,'FontSize',params.axisfontsize);
xlabel(['MEAN number of ' xtitl ' [' xunit ']']);
ylabel(['MEAN number of ' ytitl ' [' yunit ']']);
axis([0 maxx 0 maxy]);
titl = [xtitl ' <-> ' ytitl];
set(FID,'Name',titl); mtit(titl,'FontSize',params.axisfontsize+2,'yoff',.0);
if params.pdf, print(['-f' num2str(FID)],'-dpsc2','-append',params.PSFileN); end

end