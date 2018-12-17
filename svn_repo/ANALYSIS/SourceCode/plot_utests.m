%% Plot_utests
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
%% Perform Kruskal-Wallis ANOVA and Mann-Whitney U-test including 
%% Bonferroni correction on provided data and plot asterisks 
%% on top of a bar plot

% PERFORM AND PLOT STATISTICAL TESTS:
% K-W ANOVA AND PAIRWISE MW U-TESTS WITH BONF. CORRECTION
function plot_utests(utest,dataarr,labels,ys,maxy,params)

% KRUSKAL-WALLIS ANOVA
p = kruskalwallis(dataarr,labels,'off');

npairs = length(utest.ind_utest);
npairsk = length(utest.ind_over);

% TEST ALL POSSIBLE PAIRWISE COMBINATIONS, IN CASE NO OTHER
% INFORMATION WAS PROVIDED
if ~npairs && ~numel(utest.ind_kirop) && npairsk > 1,
    utest.ind_utest = zeros(length(utest.ind_over),2); icnt = 0;
    for i=1:length(utest.ind_over),
        for ii=i+1:length(utest.ind_over),
            icnt = icnt + 1;
            utest.ind_utest(icnt,:) = [i ii];
        end
    end
    npairs = length(utest.ind_utest);
end

% COMBINATORIAL NUMBER OF COMPARISIONS FOR BONFERRONI CORRECTION
if params.bonf,
    bonf = sum(1:npairs-1);
else
    bonf = 1;
end

if npairs,
    pval = zeros(1,npairs);
    for i=1:npairs,
        x1 = dataarr(labels == utest.ind_utest(i,1));
        y1 = dataarr(labels == utest.ind_utest(i,2));
        if numel(x1),
            % MANN-WHITNEY U-TEST
            pval(i) = ranksum(x1,y1);
        else
            pval(i) = 99;
        end
    end
    pmy = 0.9;
    for i=1:npairs,
        % IF (P < 0.05/BONF) PLOT ASTERISKS AND BRACKETS ON TOP OF GIVEN BARPLOT 
        if pval(i) < 0.05/bonf,
            pmy = pmy * 0.95;
            ind = utest.ind_utest(i,:);
            plot([ind(1) ind(1)],[pmy*maxy pmy*maxy-0.05*(maxy-(ys.mea(ind(1))+ys.sem(ind(1))))],'k');
            plot([ind(2) ind(2)],[pmy*maxy pmy*maxy-0.05*(maxy-(ys.mea(ind(2))+ys.sem(ind(2))))],'k');
            plot(ind,[maxy maxy]*pmy,'k');
            dx = ((ind(2) - ind(1))-.8)/2;
            patch([ind(1)+dx ind(2)-dx ind(2)-dx ind(1)+dx ind(1)+dx],...
                maxy*[(pmy-0.01) (pmy-0.01) (pmy+0.01) (pmy+0.01) (pmy-0.01)],[1 1 1],'EdgeColor','none');
            if (pval(i) <= 0.001/bonf), 
                stars = '\bullet\bullet\bullet'; 
            elseif (pval(i) > 0.001/bonf) && (pval(i) <= 0.01/bonf), 
                stars = '\bullet\bullet'; 
            else
                stars = '\bullet'; 
            end
            text(mean(ind),pmy*maxy,stars,'FontSize',10,'Color',[0 0 0],'HorizontalAlignment','center');
        end
    end
end

% SIMILAR TO BEFORE BUT WIRED FOR COMPARISONS WITH RESPECT TO 
% A CONTROL LINE ("kir")
if npairsk && numel(utest.ind_kirop),
    pvalk = zeros(1,npairsk);
    for i=1:npairsk,
        x1 = dataarr(labels == utest.ind_kirop);
        y1 = dataarr(labels == utest.ind_over(i));
        if numel(x1),
            pvalk(i) = ranksum(x1,y1); % Mann-Whitney U-Test
        else
            pvalk(i) = 99;
        end
    end
    pmy = 0.85;
    for i=1:npairsk,
        if pvalk(i) < 0.05/bonf,
            pmy = pmy * 0.95;
            ind = utest.ind_over(i); if (ind > utest.ind_kirop), ind = [ind utest.ind_kirop]; else ind = [utest.ind_kirop ind]; end
            plot([ind(1) ind(1)],[pmy*maxy pmy*maxy-0.05*(maxy-(ys.mea(ind(1))+ys.sem(ind(1))))],'k');
            plot([ind(2) ind(2)],[pmy*maxy pmy*maxy-0.05*(maxy-(ys.mea(ind(2))+ys.sem(ind(2))))],'k');
            plot(ind,[maxy maxy]*pmy,'k');
            dx = ((ind(2) - ind(1))-.8)/2;
            patch([ind(1)+dx ind(2)-dx ind(2)-dx ind(1)+dx ind(1)+dx],...
                maxy*[(pmy-0.01) (pmy-0.01) (pmy+0.01) (pmy+0.01) (pmy-0.01)],[1 1 1],'EdgeColor','none');
            if (pvalk(i) <= 0.001/bonf), 
                stars = '\bullet\bullet\bullet'; 
            elseif (pvalk(i) > 0.001/bonf) && (pvalk(i) <= 0.01/bonf), 
                stars = '\bullet\bullet'; 
            else
                stars = '\bullet';
            end
            text(mean(ind),pmy*maxy,stars,'FontSize',10,'Color',[0 0 0],'HorizontalAlignment','center');
        end
    end
end
text(max(labels),.95*maxy,['p\leq' num2str(p,'%4.3f')],'FontSize',10,'Color',[0 0 0],'HorizontalAlignment','right');
end