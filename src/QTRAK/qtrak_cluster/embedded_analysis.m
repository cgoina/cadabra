%% embedded_analysis
% Copyright (C) 2010 Heiko Dankert, California Institute of Technology
%                    Jinyang Liu, Howard Hughes Medical Institute

% This file is part of qtrak_cluster


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
% * Implementation by Heiko Dankert, Jinyang Liu
%
%% MAIN PROGRAM, called by the ANALYSIS GUI (analysis.m)

function embedded_analysis(path,name,params)

% WIDTH OF BORDER THAT MAY BE EXCLUDED FROM ANALYSIS
% analysis.radius = -99; % default (no circular arena)
% analysis.tune.correct_orient = 0; % correct fly orientation
% analysis.tune.correct_positions = 0; % find most likely fly path
params.tune = default_parameters;
params.tune.reanalyze = 1;
params.tune.correct_orient = params.correct_orient;
params.tune.correct_positions = params.correct_positions;
params.tune.lunging.thresh = params.tuningthreshold;
params.analyze_new = 1;

if ispc, params.slash = '\'; else params.slash = '/'; end


%analysis GUI parameters
%params.max_frames = params.timelimit;
params.border = 0;  %include border



% PARAMETERS
% Barplot Parameters
% params.flycol.win = [.87 .76 .54]; 
% params.flycol.los = [.71 .82 .49];
% params.flycol.ext = [0.5 .8 .76];
% params.barwidth = .3; params.twobarswidth = .18;
% params.bardispl = .2; params.semwidth = .08;
% params.whisker = 1.8;

% Transition Matrix Threshold for considered transitions into another behavior
params.trans_thr = 10;  % [s]

% ANALYZE DATA (detect and extract behavioral data)
analyze_data(path,name,params);


%%
end