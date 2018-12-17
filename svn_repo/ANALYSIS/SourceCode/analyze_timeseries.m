%% Analyze_timeseries
% Copyright (C) 2008 Heiko Dankert, California Institute of Technology

% This file is part of QTRAK.
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
% * Heiko Dankert 03/2009
%
%% Plot time series (mean and s.e.m.) of various data

function FID = analyze_timeseries(obj1,obj2,nmovs,gen_ident,params,FID)

if params.plots.timeseries,
    if numel(find(params.plots.timeseries_vec == 1)),
        % PLOT VELOCITY TIME SERIES

        minv = 0; maxt = params.max_frames; dtp = 300; % dtp = 5 min.
        t = minv:dtp:maxt; xtick_lab = cell(1,length(t));
        for i=1:length(t), xtick_lab{i} = num2str(t(i)/60); end

        ngens = numel(gen_ident);
        COL1 = zeros(3,ngens); COL2 = COL1; col = lines(ngens)';
        for i=1:ngens,
            COL1(:,i) = (0.2+col(:,i)*0.6)*i/ngens; COL2(:,i) = (0.6+COL1(:,i)*0.3)*i/ngens;
        end
        t = cell(ngens,1); vm = t; vstd = t;
        maxv = 0;

        for igen=1:ngens,
            ind = [0 find((obj1{igen}.t(2:end)-obj1{igen}.t(1:end-1))<0) ...
                length(obj1{igen}.t)];
            lmov = ind(2:end)-ind(1:end-1);
            lmax = max(lmov); imax = find(lmov == lmax); imax = imax(1);
            t0 = 1:lmax;
            t{igen} = obj1{igen}.t(ind(imax)+1:ind(imax+1));
            dt = double(median(t{igen}(2:end)-t{igen}(1:end-1)));
            v0 = nan(lmax,nmovs(igen));
            for imov=1:nmovs(igen)
                %         bin_size = ceil(window_size/dt);
                %         v0(1:lmov(imov),imov) = ...
                %             filtfilt(ones(1,bin_size)/bin_size,1,...
                %             obj1{igen}.v(ind(imov)+1:ind(imov+1)) + ...
                %             obj2{igen}.v(ind(imov)+1:ind(imov+1)))';
                if params.oneobj,
                    v0(1:lmov(imov),imov) = ...
                        obj1{igen}.vs(ind(imov)+1:ind(imov+1))';
                else
                v0(1:lmov(imov),imov) = ...
                    obj1{igen}.vs(ind(imov)+1:ind(imov+1))' + ...
                    obj2{igen}.vs(ind(imov)+1:ind(imov+1))';
                end
            end
            vstd0 = nanstd(v0,[],2)/sqrt(size(v0,2)); 
            vm0 = nanmean(v0,2);
            ivals = ~isnan(vm0);
            t{igen} = t0(ivals)*dt;
            vm{igen} = vm0(ivals)'; vstd{igen} = vstd0(ivals)';
            if max(vm{igen}+vstd{igen}) > maxv,
                maxv = max(vm{igen} + vstd{igen});
            end
        end

        % Autoscale histograms?
        if ~params.autoscale,
            maxv = 30; % mm/s
        end

        % PLOT
        if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
        figure(FID); clf;
        set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
        hold on;
        % Plot for legend
        for igen=1:ngens,
            plot(t{igen},vm{igen},'Color',COL1(:,igen)','LineWidth',1);
            hold on;
        end
        [lh,oh] = legend(gen_ident,'Location','NorthEast');
        set(oh,'LineWidth',4); legend('boxoff');
        % Plot mean and s.e.m.
        for igen=1:ngens,
            meastd = vm{igen}+vstd{igen};
            fill([t{igen} t{igen}(end:-1:1)],...
                [vm{igen}-vstd{igen} meastd(end:-1:1)],...
                COL2(:,igen)','EdgeColor','none');
            plot(t{igen},vm{igen},'Color',COL1(:,igen)','LineWidth',4);
        end
        set(gca,'XTick',0:dtp:maxt); set(gca,'XTickLabel',xtick_lab);
        xlabel('time [min]');
        ylim([0 maxv]); ylabel('velocity [mm/s]');
        titl = 'Fly velolcity time series (mean & sem)';
        set(FID,'Name',titl);
        mtit(titl,'FontSize',params.axisfontsize+2,'yoff',.04);
        if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
    end
end