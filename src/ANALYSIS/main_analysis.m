%% Main_analysis
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
%% MAIN PROGRAM, called by the ANALYSIS GUI (analysis.m)

function FID = main_analysis(path,sorted_names,mov_count,gen_ident,gen_ident1,genotypes,chamber,params)

% PARAMETERS
% Barplot Parameters
params.flycol.win = [.87 .76 .54]; 
params.flycol.los = [.71 .82 .49];
params.flycol.ext = [0.5 .8 .76];
params.barwidth = .3; params.twobarswidth = .18;
params.bardispl = .2; params.semwidth = .08;
params.whisker = 1.8;

% Transition Matrix Threshold for considered transitions into another behavior
params.trans_thr = 10;  % [s]

% ANALYZE DATA (detect and extract behavioral data)
[obj1,obj2,lunges,wings,chases,courts,tussls,copu,dists,proxi,nmovs,dts,nframes] = ...
    analyze_data(path,sorted_names,gen_ident,genotypes,params);

% Replacing genotype-code by other name/code
if exist('gen_ident1'), gen_ident = gen_ident1; end  %#ok<EXIST>

% Border for wing extension, circling heatmaps (zoom into the arenas' center)
dxy.x = 0; dxy.y = 0; % [mm]

% Miscellaneous
if ispc, params.slash = '\'; else params.slash = '/'; end
FID = params.fid;

% =======================================================

% FLY PROXIMITY MAPS
FID = analyze_proximities(obj1,obj2,dists,lunges,tussls,wings,chases,courts,dts,gen_ident,genotypes,sorted_names,params,chamber,FID);

% POSITIONS
FID = analyze_positions(obj1,obj2,dts,nframes,gen_ident,params,chamber,FID);

% VELOCITY HEATMAPS
FID = plot_velocity_heatmaps(obj1,obj2,gen_ident,params,chamber,path,FID);

% CHASING
FID = analyze_chasing(chases,lunges,dts,nframes,dxy,gen_ident,gen_ident1,params,chamber,FID);

% CIRCLING
FID = analyze_circling(courts,lunges,dts,nframes,dxy,gen_ident,gen_ident1,params,chamber,FID);

% DISTANCES
FID = analyze_distances(obj1,obj2,lunges,copu,dists,proxi,nmovs,dts,gen_ident,gen_ident1,params,chamber,FID);

% MOVE-STOP ANALYSIS
FID = plot_move_stop_trans(obj1,obj2,lunges,gen_ident,params,dts,FID);

% FLY BODY SIZE ANALYSIS
FID = analyze_flybodysize(obj1,obj2,lunges,nmovs,gen_ident,gen_ident1,params,FID);

% AGGRESSION ANALYSIS
FID = analyze_aggression(obj1,obj2,lunges,wings,chases,courts,tussls,copu,nmovs,dts,nframes,gen_ident,gen_ident1,params,chamber,FID);

% WING ANALYSIS
FID = analyze_wings(lunges,wings,copu,dists,dts,nframes,dxy,path,sorted_names,gen_ident,genotypes,params,chamber,FID);

% VELOCITY ANALYSIS
FID = analyze_velocities(obj1,obj2,lunges,proxi,dts,nframes,gen_ident,params,FID);

% COPULATION ANALYSIS
FID = analyze_copulation(wings,copu,dts,gen_ident,gen_ident1,params,FID,mov_count);

% PERFORMANCE INDICES
FID = action_occurrence(lunges,tussls,wings,chases,courts,nmovs,dts,gen_ident,gen_ident1,params,FID);

% TIME SERIES ANALYSIS
FID = analyze_timeseries(obj1,obj2,nmovs,gen_ident,params,FID);

% ETHOGRAMS
FID = ethograms(path,gen_ident,gen_ident1,params,FID);

% CLUSTER ANALYSIS
FID = analyze_clusters(path,gen_ident,params,FID);

% MOVIE CLIPS
plot_clips(path,sorted_names,gen_ident1,genotypes,params);

%%
end