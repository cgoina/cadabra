%% Analyze_chasing
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
%% Analyze chasing

function FID = analyze_chasing(chases,lunges,dts,nframes,dxy,gen_ident,gen_ident1,params,chamber,FID)

% CHASING
if ~params.oneobj,
    if params.plots.heatmaps.chasing,
        % CHASING HEATMAPS AND RASTERPLOTS
        if params.tim,
            max_val = params.range.heatmaps.tim.chase;
        else
            max_val = params.range.heatmaps.occ.chase;
        end
        min_bout = -99; max_gap = -99;
        plot_oppon = 1;
        if numel(find(params.plots.heatmaps.chase_vec == 1)),
            titl = 'Histograms of Fly Chasing Positions';   
            FID = plot_heatmaps(chases,plot_oppon,titl,max_val,min_bout,max_gap,chamber,gen_ident,dxy,dts,nframes,FID,params);
        end
        if numel(find(params.plots.heatmaps.chase_vec == 2)),
            titl = 'Chasing Distribution';
            FID = plot_2ddistr(chases,titl,max_val,gen_ident,dxy,dts,nframes,FID,params);
        end
    end
    if params.plots.stat.chasing && ~params.plots.movieclips,
        % CHASING 1D HISTOGRAM OVER GENOTYPES
        bool_sumobj = 0;
        titl = 'Chases'; ylab = 'chases';
        if params.tim,
            max_val = params.range.stat.tim.chase;
        else
            max_val = params.range.stat.occ.chase;
        end
        FID = plot_loser_winner_stat(chases,FID,max_val,gen_ident,params,ylab,dts,titl,bool_sumobj,lunges,0);

        utest = utestpairs(gen_ident1); bool_sumobj = 1;
        titl = 'Chases'; ylab = 'chases';
        FID = plot_loser_winner_stat(chases,FID,max_val,gen_ident,params,ylab,dts,titl,bool_sumobj,lunges,utest);
    end
end
