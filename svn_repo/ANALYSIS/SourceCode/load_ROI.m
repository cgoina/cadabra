%% Load_ROI
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
%% Load ROI from tracked movie

% LOAD ROI AND SCALE FROM TRACKING DATA
function [ROI,scale,fname] = load_ROI(feat_filen)

fname = feat_filen(1:end-5);
in = strfind(fname,'_'); fname = fname(1:in(end)-1);
fid = fopen([feat_filen(1:end-5) '_roi.mat'],'r');
load([feat_filen(1:end-5) '_roi.mat']);
fclose(fid);
end