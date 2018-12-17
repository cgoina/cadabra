%% Analyze_circling
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
%% Analyze circling

function FID = analyze_circling(courts,lunges,dts,nframes,dxy,gen_ident,gen_ident1,params,chamber,FID)
% CIRCLING
if params.plots.heatmaps.circling,
    % CIRCLING HEATMAPS AND RASTERPLOTS
    if params.tim,
        max_val = params.range.heatmaps.tim.circl;
    else
        max_val = params.range.heatmaps.occ.circl;
    end    
    plot_oppon = 0;
    if numel(find(params.plots.heatmaps.circl_vec == 1)),
        titl = 'Histograms of Fly Circling Positions';
        FID = plot_heatmaps(courts,plot_oppon,titl,max_val,courts{1}.min_bout,courts{1}.max_gap,chamber,gen_ident,dxy,dts,nframes,FID,params);
    end
    if numel(find(params.plots.heatmaps.circl_vec == 2)),
        titl = 'Circling Distribution';
        FID = plot_2ddistr(courts,titl,max_val,gen_ident,dxy,dts,nframes,FID,params);
    end
end
if params.plots.stat.circling && ~params.plots.movieclips,
    % CIRCLING 1D HISTOGRAM OVER GENOTYPES
    utest = utestpairs(gen_ident1);
    titl = 'Circling'; ylab = 'circling';
    if params.tim,
        max_val = params.range.stat.tim.circl;
    else
        max_val = params.range.stat.occ.circl;
    end
    bool_sumobj = 0;
    FID = plot_loser_winner_stat(courts,FID,max_val,gen_ident,params,ylab,dts,titl,bool_sumobj,lunges,utest);
    bool_sumobj = 1;
    FID = plot_loser_winner_stat(courts,FID,max_val,gen_ident,params,ylab,dts,titl,bool_sumobj,lunges,utest);
end