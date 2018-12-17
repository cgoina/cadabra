%% Action_occurrence
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
%% Compute occurence and time spent in aggressive, chasing, and courtship 
%% actions (absolute and percentage over all actions)

function FID = action_occurrence(lunges,tussls,wings,chases,courts,nmovs,dts,gen_ident,gen_ident1,params,FID)
if params.plots.stat.actionrel,
    bool_wingthr = 1;
    bool_chase = numel(find(params.plots.stat.actionrel_vec == 1));
    perctot = [numel(find(params.plots.stat.actionrel_vec == 2)) ...
               numel(find(params.plots.stat.actionrel_vec == 3))];
    
    bool_wing_circ = 0; % 1 = only wing extension; 2 = circling; 0 = both
    
    ngens = length(lunges);
    for ii=1:numel(perctot),
        if perctot(ii),
            if ii == 2, bool_percentage = 0; else bool_percentage = 1; end
            y_aggr.mea = zeros(ngens,1); y_aggr.sem = y_aggr.mea;
            y_aggr.min = y_aggr.mea; y_aggr.max = y_aggr.mea; sy_aggr = y_aggr;
            y_chase = y_aggr; sy_chase = y_aggr;
            y_court = y_aggr; sy_court = y_aggr;
            y_wthraggr = y_aggr; sy_wthraggr = y_aggr;
            y_wthraggrchs = y_aggr; sy_wthraggrchs = y_aggr;
            X_aggr = []; G_aggr = []; X_chase = []; G_chase = []; X_court = []; G_court = [];
            X_wthraggr = []; G_wthraggr = []; X_wthraggrchs = []; G_wthraggrchs = [];
            sX_aggr = []; sG_aggr = []; sX_chase = []; sG_chase = []; sX_court = []; sG_court = [];
            sX_wthraggr = []; sG_wthraggr = []; sX_wthraggrchs = []; sG_wthraggrchs = [];
            dataarr.aggr = cell(ngens,1); dataarr.chase = cell(ngens,1); dataarr.court = cell(ngens,1); 
            dataarr.wingthr = cell(ngens,1); dataarr.wingthr_over_aggrchase = cell(ngens,1);
            dataarr.aggr_l = cell(ngens,1); dataarr.chase_l = cell(ngens,1); dataarr.court_l = cell(ngens,1); 
            dataarr.wingthr_l = cell(ngens,1); dataarr.wingthr_l_over_aggrchase = cell(ngens,1);
            
            for igen=1:ngens,
                aggr = []; chase = []; court = []; aggr_l = []; chase_l = []; court_l = [];
                wingthr = []; wingthr_l = [];
                cl1 = [0 cumsum(lunges{igen}.obj1.number)]; cl2 = [0 cumsum(lunges{igen}.obj2.number)];
                ct = [0 cumsum(tussls{igen}.number)];
                cwthr1 = [0 cumsum(wings.threat{igen}.obj1.number)]; cwthr2 = [0 cumsum(wings.threat{igen}.obj2.number)];
                cch1 = [0 cumsum(chases{igen}.obj1.number)]; cch2 = [0 cumsum(chases{igen}.obj2.number)];
                cwextl1 = [0 cumsum(wings.ext.l{igen}.obj1.number)]; cwextl2 = [0 cumsum(wings.ext.l{igen}.obj2.number)];
                cwextr1 = [0 cumsum(wings.ext.r{igen}.obj1.number)]; cwextr2 = [0 cumsum(wings.ext.r{igen}.obj2.number)];
                ccrt1 = [0 cumsum(courts{igen}.obj1.number)]; ccrt2 = [0 cumsum(courts{igen}.obj2.number)];
                for imov=1:nmovs(igen),
                    % Compute sums of the number of wing threat, aggressive, chasing, and 
                    % courtship actions for each genotype
                    wingthr = [wingthr wings.threat{igen}.obj1.number(imov)+wings.threat{igen}.obj2.number(imov)];
                    aggr = [aggr lunges{igen}.number(imov)+tussls{igen}.number(imov)+ ...
                        wings.threat{igen}.obj1.number(imov)+wings.threat{igen}.obj2.number(imov)];
                    chase = [chase chases{igen}.obj1.number(imov)+chases{igen}.obj2.number(imov)];
                    if bool_wing_circ == 1,
                        court = [court wings.ext.l{igen}.obj1.number(imov)+wings.ext.l{igen}.obj2.number(imov)+ ...
                            wings.ext.r{igen}.obj1.number(imov)+wings.ext.r{igen}.obj2.number(imov)];
                    elseif bool_wing_circ == 2,
                        court = [court courts{igen}.obj1.number(imov)+courts{igen}.obj2.number(imov)];
                    else
                        court = [court wings.ext.l{igen}.obj1.number(imov)+wings.ext.l{igen}.obj2.number(imov)+ ...
                            wings.ext.r{igen}.obj1.number(imov)+wings.ext.r{igen}.obj2.number(imov)+ ...
                            courts{igen}.obj1.number(imov)+courts{igen}.obj2.number(imov)];
                    end

                    % Compute sums of the duration of wing threats, aggressive, chasing, and 
                    % courtship actions for each genotype in minutes
                    swingthr = [wings.threat{igen}.obj1.len(cwthr1(imov)+1:cwthr1(imov+1)) ...
                        wings.threat{igen}.obj2.len(cwthr2(imov)+1:cwthr2(imov+1))];
                    saggr = [lunges{igen}.obj1.len(cl1(imov)+1:cl1(imov+1)) ...
                        lunges{igen}.obj2.len(cl2(imov)+1:cl2(imov+1)) ...
                        tussls{igen}.len(ct(imov)+1:ct(imov+1)) ...
                        wings.threat{igen}.obj1.len(cwthr1(imov)+1:cwthr1(imov+1)) ...
                        wings.threat{igen}.obj2.len(cwthr2(imov)+1:cwthr2(imov+1))];
                    schase = [chases{igen}.obj1.len(cch1(imov)+1:cch1(imov+1)) ...
                        chases{igen}.obj2.len(cch2(imov)+1:cch2(imov+1))];
                    if bool_wing_circ == 1,
                        scourt = [wings.ext.l{igen}.obj1.len(cwextl1(imov)+1:cwextl1(imov+1)) ...
                            wings.ext.l{igen}.obj2.len(cwextl2(imov)+1:cwextl2(imov+1)) ...
                            wings.ext.r{igen}.obj1.len(cwextr1(imov)+1:cwextr1(imov+1)) ...
                            wings.ext.r{igen}.obj2.len(cwextr2(imov)+1:cwextr2(imov+1))];
                    elseif bool_wing_circ == 2,
                        scourt = [courts{igen}.obj1.len(ccrt1(imov)+1:ccrt1(imov+1)) ...
                            courts{igen}.obj2.len(ccrt2(imov)+1:ccrt2(imov+1))];
                    else
                        scourt = [wings.ext.l{igen}.obj1.len(cwextl1(imov)+1:cwextl1(imov+1)) ...
                            wings.ext.l{igen}.obj2.len(cwextl2(imov)+1:cwextl2(imov+1)) ...
                            wings.ext.r{igen}.obj1.len(cwextr1(imov)+1:cwextr1(imov+1)) ...
                            wings.ext.r{igen}.obj2.len(cwextr2(imov)+1:cwextr2(imov+1)) ...
                            courts{igen}.obj1.len(ccrt1(imov)+1:ccrt1(imov+1)) ...
                            courts{igen}.obj2.len(ccrt2(imov)+1:ccrt2(imov+1))];
                    end
                    wingthr_l = [wingthr_l sum(swingthr)*dts{igen}(imov)/60];
                    aggr_l = [aggr_l sum(saggr)*dts{igen}(imov)/60];
                    chase_l = [chase_l sum(schase)*dts{igen}(imov)/60];
                    court_l = [court_l sum(scourt)*dts{igen}(imov)/60];
                end
                                
                % OCCURRENCES
                % Total aggressive and courtship actions
                total = aggr + court;
                % Include chasing?
                if bool_chase, total = total + chase; end
                
                % Aggression
                % Compute absolute numbers or percentages?
                if ~bool_percentage,
                    data = aggr;
                else
                    data = aggr./total; data = data(~isnan(data))*100;
                end
                
                % Compute statistics (mean, sem, min, max); prepare data for boxplots (X,G)
                y_aggr.mea(igen) = mean(data); y_aggr.sem(igen) = std(data)/sqrt(nmovs(igen));
                y_aggr.min(igen) = min(data); y_aggr.max(igen) = max(data);
                X_aggr = [X_aggr data]; G_aggr = [G_aggr ones(1,length(data))*igen];

                % Prepare data for text output
                dataarr.aggr{igen} = data;
                
                % Chasing
                % Compute absolute numbers or percentages?
                if bool_chase,
                    if ~bool_percentage,
                        data = chase;
                    else
                        data = chase./total; data = data(~isnan(data))*100;
                    end
                    
                    % Statistics
                    y_chase.mea(igen) = mean(data); y_chase.sem(igen) = std(data)/sqrt(nmovs(igen));
                    y_chase.min(igen) = min(data); y_chase.max(igen) = max(data);
                    X_chase = [X_chase data]; G_chase = [G_chase ones(1,length(data))*igen];

                    % Prepare data for text output
                    dataarr.chase{igen} = data;
                end
                
                % Courtship
                % Compute absolute numbers or percentages?
                if ~bool_percentage,
                    data = court;
                else
                    data = court./total; data = data(~isnan(data))*100;
                end
                
                % Statistics
                y_court.mea(igen) = mean(data); y_court.sem(igen) = std(data)/sqrt(nmovs(igen));
                y_court.min(igen) = min(data); y_court.max(igen) = max(data);
                X_court = [X_court data]; G_court = [G_court ones(1,length(data))*igen];

                % Prepare data for text output
                dataarr.court{igen} = data;
                
                % Wing threats
                % Compute absolute numbers or percentages?
                if ~bool_percentage,
                    data = wingthr;
                else
                    data = wingthr./aggr; data = data(~isnan(data))*100;
                end
                
                % Statistics
                y_wthraggr.mea(igen) = mean(data); y_wthraggr.sem(igen) = std(data)/sqrt(nmovs(igen));
                y_wthraggr.min(igen) = min(data); y_wthraggr.max(igen) = max(data);
                X_wthraggr = [X_wthraggr data]; G_wthraggr = [G_wthraggr ones(1,length(data))*igen];

                % Prepare data for text output
                dataarr.wingthr{igen} = data;

                % Wing threats over aggression + chasing
                data = wingthr./(aggr+chase); data = data(~isnan(data))*100;
                y_wthraggrchs.mea(igen) = mean(data); y_wthraggrchs.sem(igen) = std(data)/sqrt(nmovs(igen));
                y_wthraggrchs.min(igen) = min(data); y_wthraggrchs.max(igen) = max(data);
                X_wthraggrchs = [X_wthraggrchs data]; G_wthraggrchs = [G_wthraggrchs ones(1,length(data))*igen];

                % Prepare data for text output
                dataarr.wingthr_over_aggrchase{igen} = data;


                % DURATIONS
                % Total aggressive and courtship action durations
                total = aggr_l + court_l;
                % Include chasing?
                if bool_chase, total = total + chase_l; end
                
                % Aggression
                % Compute absolute numbers or percentages?
                if ~bool_percentage,
                    data = aggr_l;
                else
                    data = aggr_l./total; data = data(~isnan(data))*100;
                end
                
                % Compute statistics (mean, sem, min, max); prepare data for boxplots (X,G)
                sy_aggr.mea(igen) = mean(data); sy_aggr.sem(igen) = std(data)/sqrt(nmovs(igen));
                sy_aggr.min(igen) = min(data); sy_aggr.max(igen) = max(data);
                sX_aggr = [sX_aggr data]; sG_aggr = [sG_aggr ones(1,length(data))*igen];
                                
                % Prepare data for text output
                dataarr.aggr_l{igen} = data;
                
                % Chasing
                % Compute absolute numbers or percentages?
                if bool_chase,
                    if ~bool_percentage,
                        data = chase_l;
                    else
                        data = chase_l./total; data = data(~isnan(data))*100;
                    end
                    
                    % Statistics
                    sy_chase.mea(igen) = mean(data); sy_chase.sem(igen) = std(data)/sqrt(nmovs(igen));
                    sy_chase.min(igen) = min(data); sy_chase.max(igen) = max(data);
                    sX_chase = [sX_chase data]; sG_chase = [sG_chase ones(1,length(data))*igen];    
                    
                    % Prepare data for text output
                    dataarr.chase_l{igen} = data;
                end

                
                % Courtship
                % Compute absolute numbers or percentages?
                if ~bool_percentage,
                    data = court_l;
                else
                    data = court_l./total; data = data(~isnan(data))*100;
                end
                
                % Statistics
                sy_court.mea(igen) = mean(data); sy_court.sem(igen) = std(data)/sqrt(nmovs(igen));
                sy_court.min(igen) = min(data); sy_court.max(igen) = max(data);
                sX_court = [sX_court data]; sG_court = [sG_court ones(1,length(data))*igen];
                
                % Prepare data for text output
                dataarr.court_l{igen} = data;
                
                % Wing threats
                % Compute absolute numbers or percentages?
                if ~bool_percentage,
                    data = wingthr_l;
                else
                    data = wingthr_l./aggr_l; data = data(~isnan(data))*100;
                end
                
                % Statistics
                sy_wthraggr.mea(igen) = mean(data); sy_wthraggr.sem(igen) = std(data)/sqrt(nmovs(igen));
                sy_wthraggr.min(igen) = min(data); sy_wthraggr.max(igen) = max(data);
                sX_wthraggr = [sX_wthraggr data]; sG_wthraggr = [sG_wthraggr ones(1,length(data))*igen];

                % Prepare data for text output
                dataarr.wingthr_l{igen} = data;

                % Wing threats over aggression + chasing
                data = wingthr_l./(aggr_l+chase_l); data = data(~isnan(data))*100;
                sy_wthraggrchs.mea(igen) = mean(data); sy_wthraggrchs.sem(igen) = std(data)/sqrt(nmovs(igen));
                sy_wthraggrchs.min(igen) = min(data); sy_wthraggrchs.max(igen) = max(data);
                sX_wthraggrchs = [sX_wthraggrchs data]; sG_wthraggrchs = [sG_wthraggrchs ones(1,length(data))*igen];

                % Prepare data for text output
                dataarr.wingthr_l_over_aggrchase{igen} = data;                
            end
            % Detemine genotype pairs to be statistically compared
            utest = utestpairs(gen_ident1);
            
            % PLOTS
            if ~bool_percentage,
                if ~params.tim,
                    maxy = 200;
                else
                    maxy = 1;
                end
            else
                maxy = 100;
            end
            if bool_wing_circ == 1,
                wing_circ_court = 'wing extension';
            elseif bool_wing_circ == 2,
                wing_circ_court = 'circling';
            else
                wing_circ_court = 'courtship';
            end

            if bool_chase, txt = 'all actions'; else txt = ['aggression and ' wing_circ_court]; end
            if ~params.tim,
                ylab = 'aggression';
                if bool_percentage, 
                    titl = ['Proportion of aggressive actions over ' txt]; 
                    ylab = [ylab ' / ' txt ' [%]']; 
                    txtsuffix = ['over ' txt ' [percent]']; 
                else
                    titl = 'Occurrence of aggressive actions'; 
                    txtsuffix = 'number';
                end
                if params.bool_xls, textoutput(dataarr.aggr,gen_ident,txtsuffix,'Aggressive-actions',params); end
                FID = plot_bar_msem(y_aggr,utest,X_aggr,G_aggr,maxy,ylab,titl,FID,gen_ident,params);
                if bool_chase,
                    ylab = 'chasing';
                    if bool_percentage, 
                        titl = ['Proportion of chasing over ' txt]; 
                        ylab = [ylab ' / ' txt ' [%]']; 
                        txtsuffix = ['over ' txt ' [percent]']; 
                    else
                        titl = 'Occurrence of chasing';
                        txtsuffix = 'number';
                    end
                    if params.bool_xls, textoutput(dataarr.chase,gen_ident,txtsuffix,'Chasing-actions',params); end
                    FID = plot_bar_msem(y_chase,utest,X_chase,G_chase,maxy,ylab,titl,FID,gen_ident,params);
                end
                ylab = wing_circ_court;
                if bool_percentage, 
                    titl = ['Proportion of ' wing_circ_court ' over ' txt]; 
                    ylab = [ylab ' / ' txt ' [%]']; 
                    txtsuffix = ['over ' txt ' [percent]'];
                else
                    titl = ['Occurrence of ' wing_circ_court]; 
                    txtsuffix = 'number';
                end
                if params.bool_xls, textoutput(dataarr.court,gen_ident,txtsuffix,'Courtship-actions',params); end
                FID = plot_bar_msem(y_court,utest,X_court,G_court,maxy,ylab,titl,FID,gen_ident,params);

                if bool_wingthr,
                    ylab = 'wing threat';
                    if ~bool_percentage,
                        titl = 'Occurrence of wing threat'; 
                        txtsuffix = 'number';
                        if params.bool_xls, textoutput(dataarr.wingthr,gen_ident,txtsuffix,'Wingthreat-actions',params); end
                        FID = plot_bar_msem(y_wthraggr,utest,X_wthraggr,G_wthraggr,maxy,ylab,titl,FID,gen_ident,params);
                    else
                        titl = 'Proportion of wing threat over aggression and chasing'; 
                        if params.bool_xls, textoutput(dataarr.wingthr_over_aggrchase,gen_ident,'over aggr chasing [percent]','Wingthreat-actions',params); end
                        ylab = ['wing threat / (aggression + chasing) [%]'];
                        FID = plot_bar_msem(y_wthraggrchs,utest,X_wthraggrchs,G_wthraggrchs,maxy,ylab,titl,FID,gen_ident,params);
                    end
                end
            else
                ylab = 'aggression';
                if bool_percentage, 
                    titl = ['Proportion of time spent in aggression over ' txt]; 
                    ylab = [ylab ' / ' txt ' [%]']; 
                    txtsuffix = ['over ' txt ' time [percent]']; 
                else
                    titl = 'Time spent in aggressive actions'; 
                    ylab = [ylab ' [min]']; 
                    txtsuffix = 'time spent [min]';
                end
                if params.bool_xls, textoutput(dataarr.aggr_l,gen_ident,txtsuffix,'Aggressive-actions',params); end
                FID = plot_bar_msem(sy_aggr,utest,sX_aggr,sG_aggr,maxy,ylab,titl,FID,gen_ident,params);
                if bool_chase,
                    ylab = 'chasing';
                    if bool_percentage,
                        titl = ['Proportion of time spent in chasing over ' txt]; 
                        ylab = [ylab ' / ' txt ' [%]']; 
                        txtsuffix = ['over ' txt ' time [percent]'];
                    else
                        titl = 'Time spent in chasing'; 
                        ylab = [ylab ' [min]'];
                        txtsuffix = 'time spent [min]';
                    end
                    if params.bool_xls, textoutput(dataarr.chase_l,gen_ident,txtsuffix,'Chasing-actions',params); end
                    FID = plot_bar_msem(sy_chase,utest,sX_chase,sG_chase,maxy,ylab,titl,FID,gen_ident,params);
                end
                ylab = wing_circ_court;
                if bool_percentage, 
                    titl = ['Proportion of time spent in ' wing_circ_court ' over ' txt]; 
                    ylab = [ylab ' / ' txt ' [%]']; 
                    txtsuffix = ['over ' txt ' time [percent]']; 
                else
                    titl = ['Time spent in ' wing_circ_court]; 
                    ylab = [ylab ' [min]']; 
                    txtsuffix = 'time spent [min]';
                end
                if params.bool_xls, textoutput(dataarr.court_l,gen_ident,txtsuffix,'Courtship-actions',params); end
                FID = plot_bar_msem(sy_court,utest,sX_court,sG_court,maxy,ylab,titl,FID,gen_ident,params);

                if bool_wingthr,
                    if ~bool_percentage,
                        titl = 'Time spent in wing threat'; 
                        txtsuffix = 'time spent [min]';
                        if params.bool_xls, textoutput(dataarr.wingthr_l,gen_ident,txtsuffix,'Wingthreat-actions',params); end
                        ylab = 'wing threat [min]';
                        FID = plot_bar_msem(sy_wthraggr,utest,sX_wthraggr,sG_wthraggr,maxy,ylab,titl,FID,gen_ident,params);
                    else
                        titl = 'Proportion of time spent in wing threat'; 
                        if params.bool_xls, textoutput(dataarr.wingthr_l_over_aggrchase,gen_ident,'over aggr chasing time [percent]','Wingthreat-actions',params); end
                        titl = [titl ' and chasing']; ylab = ['wing threat / (aggression + chasing) [%]'];
                        FID = plot_bar_msem(sy_wthraggrchs,utest,sX_wthraggrchs,sG_wthraggrchs,maxy,ylab,titl,FID,gen_ident,params);
                    end
                end
            end
        end
    end

    
    % COMPUTE TIME SERIES INFORMATION ON AGGRESSIVE,
    % CHASING, AND COURTSHIP ACTIONS
    if numel(find(params.plots.stat.actionrel_vec == 4)),
        FID = beh_occ_ts(wings,lunges,tussls,chases,courts,gen_ident1,ngens,nmovs,dts,params,FID);
    end
end
