%% Plot_loser_winner_stat
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
%% Plot statistics of provided genotype action data, sorted by
%% loser/winner or by fly identities

% PLOT GENOTYPE STATISTICS (BARPLOTS OR BOXPLOTS)
function FID = plot_loser_winner_stat(behavior,FID,maxy,gen_ident,params,ylab,dts,titl,bool_sumobj,lunges,utest,behavior2)

ngens = length(lunges);
y.mea = zeros(ngens,2); y.sem = y.mea; 
y.min = y.mea; y.max = y.mea; sy = y; X = []; G = []; X2 = []; G2 = []; sX = []; sG = []; sX2 = []; sG2 = [];
data = cell(ngens,1); sdata = data; data1 = data; data2 = data; sdata1 = data; sdata2 = data;

% COLLECT ACTION NUMBER AND DURATION INFORMATION
for igen=1:ngens, 
    cs1 = [0 cumsum(behavior{igen}.obj1.number)]; 
    cs2 = [0 cumsum(behavior{igen}.obj2.number)];
    if exist('behavior2'),
        cs12 = [0 cumsum(behavior2{igen}.obj1.number)]; 
        cs22 = [0 cumsum(behavior2{igen}.obj2.number)]; 
    end
    for imov=1:numel(lunges{igen}.mov),        
        % COLLECT ACTION NUMBER AND DURATION INFORMATION FOR EACH
        % MOVIE / FLY PAIR
        % Distinguish between winner / loser or fly 1 / 2?
        if (lunges{igen}.obj1.number(imov) >= lunges{igen}.obj2.number(imov)) || params.courtship  || ~params.winlos,
            % Winner / loser
            if exist('behavior2'),
                % Add a second action into the counting (i.e., right wing to left wing data)
                tmp1 = behavior{igen}.obj1.number(imov) + behavior2{igen}.obj1.number(imov);
                tmp2 = behavior{igen}.obj2.number(imov) + behavior2{igen}.obj2.number(imov);
                stmp1 = [behavior{igen}.obj1.len(cs1(imov)+1:cs1(imov+1)) behavior2{igen}.obj1.len(cs12(imov)+1:cs12(imov+1))];
                stmp2 = [behavior{igen}.obj2.len(cs2(imov)+1:cs2(imov+1)) behavior2{igen}.obj2.len(cs22(imov)+1:cs22(imov+1))];
            else
                tmp1 = behavior{igen}.obj1.number(imov);
                tmp2 = behavior{igen}.obj2.number(imov);
                stmp1 = behavior{igen}.obj1.len(cs1(imov)+1:cs1(imov+1));
                stmp2 = behavior{igen}.obj2.len(cs2(imov)+1:cs2(imov+1));
            end
        else
            % Fly 1 / 2
            if exist('behavior2'),
                tmp1 = behavior{igen}.obj2.number(imov) + behavior2{igen}.obj2.number(imov);
                tmp2 = behavior{igen}.obj1.number(imov) + behavior2{igen}.obj1.number(imov);
                stmp1 = [behavior{igen}.obj2.len(cs2(imov)+1:cs2(imov+1)) behavior2{igen}.obj2.len(cs22(imov)+1:cs22(imov+1))];
                stmp2 = [behavior{igen}.obj1.len(cs1(imov)+1:cs1(imov+1)) behavior2{igen}.obj1.len(cs12(imov)+1:cs12(imov+1))];
            else
                tmp1 = behavior{igen}.obj2.number(imov);
                tmp2 = behavior{igen}.obj1.number(imov);
                stmp1 = behavior{igen}.obj2.len(cs2(imov)+1:cs2(imov+1));
                stmp2 = behavior{igen}.obj1.len(cs1(imov)+1:cs1(imov+1));
            end
        end
        data1{igen} = [data1{igen} tmp1]; sdata1{igen} = [sdata1{igen} sum(stmp1)*dts{igen}(imov)/60];
        data2{igen} = [data2{igen} tmp2]; sdata2{igen} = [sdata2{igen} sum(stmp2)*dts{igen}(imov)/60];
    end

    % STATISTICS FOR BARPLOTS
    % Occurrences
    y.mea(igen,1) = mean(data1{igen});
    y.sem(igen,1) = std(data1{igen}) / sqrt(numel(data1{igen})); 
    y.min(igen,1) = min(data1{igen}); y.max(igen,1) = max(data1{igen});
    y.mea(igen,2) = mean(data2{igen});
    y.sem(igen,2) = std(data2{igen}) / sqrt(numel(data2{igen})); 
    y.min(igen,2) = min(data2{igen}); y.max(igen,2) = max(data2{igen});

    % Time spent
    sy.mea(igen,1) = mean(sdata1{igen});
    sy.sem(igen,1) = std(sdata1{igen}) / sqrt(numel(sdata1{igen})); 
    sy.min(igen,1) = min(sdata1{igen}); sy.max(igen,1) = max(sdata1{igen});
    sy.mea(igen,2) = mean(sdata2{igen});
    sy.sem(igen,2) = std(sdata2{igen}) / sqrt(numel(sdata2{igen})); 
    sy.min(igen,2) = min(sdata2{igen}); sy.max(igen,2) = max(sdata2{igen});

    % Sum over both flies
    if bool_sumobj,
        % Statistics for barplots
        data{igen} = data1{igen} + data2{igen};
        ys.mea(igen,1) = mean(data{igen});
        ys.sem(igen,1) = std(data{igen}) / sqrt(numel(data{igen}));
        ys.min(igen,1) = min(data{igen}); ys.max(igen,1) = max(data{igen});
        sdata{igen} = sdata1{igen} + sdata2{igen};
        sys.mea(igen,1) = mean(sdata{igen});
        sys.sem(igen,1) = std(sdata{igen}) / sqrt(numel(sdata{igen}));
        sys.min(igen,1) = min(sdata{igen}); sys.max(igen,1) = max(sdata{igen});        
        % Data for boxplots
        X = [X data{igen}]; G = [G ones(1,numel(data{igen}))*igen];
        sX = [sX sdata{igen}]; sG = [sG ones(1,numel(sdata{igen}))*igen];
    end

    % DATA FOR BOXPLOTS
    X2 = [X2 data1{igen} data2{igen}]; 
    G2 = [G2 ones(1,length(data1{igen}))*(2*igen-1) ones(1,length(data2{igen}))*2*igen];
    sX2 = [sX2 sdata1{igen} sdata2{igen}]; 
    sG2 = [sG2 ones(1,length(sdata1{igen}))*(2*igen-1) ones(1,length(sdata2{igen}))*2*igen];
end

% TEXT OUTPUT
if params.bool_xls,
    if ~params.tim
        % Occurrence number
        textoutput([data1 , data2],gen_ident,'number',titl,params);
    else
        % Time spent
        textoutput([sdata1 , sdata2],gen_ident,'time spent [min]',titl,params);
    end
    % Time points of action occurrences
    texttimeoutput(behavior,gen_ident,titl,params);
end


% PLOTS
% Both flies separated
% Occurrence statistics
if ~params.tim && ~bool_sumobj,
    FID = plot_bar_msem(y,utest,X2,G2,maxy,ylab,['Number of ' titl],FID,gen_ident,params);
end
% Time spent statistics
if params.tim && ~bool_sumobj,
    FID = plot_bar_msem(sy,utest,sX2,sG2,maxy,[ylab ' [min]'],['Total Time Spend for ' titl ' [min]'],FID,gen_ident,params);
end

% Both flies together
% Occurrence statistics
if ~params.tim && bool_sumobj,
    FID = plot_bar_msem(ys,utest,X,G,maxy,ylab,['Number of ' titl],FID,gen_ident,params);
end
% Time spent statistics
if params.tim && bool_sumobj,
    FID = plot_bar_msem(sys,utest,sX,sG,maxy,[ylab ' [min]'],['Total Time Spend for ' titl ' [min]'],FID,gen_ident,params);
end

end