%% Analyze_copulation
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
%% Analyze detected copulation and compute statistics

function FID = analyze_copulation(wings,copu,dts,gen_ident,gen_ident1,params,FID,mov_count)

if params.courtship && ~params.oneobj && params.plots.stat.copulation && ~params.plots.movieclips,
    % Plot copulation begin, duration, end
    if params.max_frames <= 1200, params.max_frames = 30*60; end % 30 minutes
    ngens = length(wings.ext.r);
    maxt = params.max_frames / 60; %[min]
    if numel(find(params.plots.stat.copul_vec == 1)),
        FID = barplot_frame(copu.pre,dts,'Begin Copulation',maxt,gen_ident,mov_count,params.cross_to,FID,params);
    end
    if numel(find(params.plots.stat.copul_vec == 2)),
        FID = barplot_frame(copu.int,dts,'Copulation Duration',maxt,gen_ident,mov_count,params.cross_to,FID,params);
    end
    if numel(find(params.plots.stat.copul_vec == 3)),
        FID = barplot_frame(copu.post,dts,'End Copulation',maxt,gen_ident,mov_count,params.cross_to,FID,params);
    end

    if sum(wings.ext.r{1}.obj1.number) > 0 && numel(find(params.plots.stat.copul_vec == 4)),
        % Latency to the first wing extension
        first_wing = cell(1,ngens);
        for igen=1:ngens,
            mov_r = [0 cumsum(wings.ext.r{igen}.obj1.number)];
            mov_l = [0 cumsum(wings.ext.l{igen}.obj1.number)];
            dat_r = wings.ext.r{igen}.obj1.t; dat_l = wings.ext.l{igen}.obj1.t;
            data = [];
            for j=1:numel(mov_r)-1,
                tmp_r = dat_r(mov_r(j)+1:mov_r(j+1)); tmp_l = dat_l(mov_l(j)+1:mov_l(j+1));
                if numel(tmp_r) && numel(tmp_l),
                    data = [data min([tmp_r(1) tmp_l(1)])/dts{igen}(j)];
                else
                    data = [data 0];
                end
            end
            first_wing{igen} = data;
        end
        FID = barplot_frame(first_wing,dts,'Latency to 1st Wing Extension',5,gen_ident,mov_count,params.cross_to,FID,params);
    end

    if numel(find(params.plots.stat.copul_vec == 5)),
        % Copulation index
        x = 1:ngens; y.mea = zeros(ngens,1); y.sem = y.mea; X = []; G = [];
        for igen=1:ngens,
            data = copu.pre{igen} > 0;
            if numel(data)>0,
                y.mea(igen) = mean(data)*100;
                y.sem(igen) = std(data)/sqrt(numel(data))*100;
                X = [X data]; G = [G ones(1,length(data))*igen];
            end
        end
        utest = utestpairs(gen_ident1);
        ytit = 'copulations [%]'; titl = 'Copulation Index';
        FID = plot_bar_msem(y,utest,X,G,100,ytit,titl,FID,gen_ident,params);
    end
    
    %JL08182009 Add the textouput per Eric Hoopfer's request 
    textoutput4copu([copu.pre, copu.post], gen_ident, 'pre_post', 'copulation', params);
end