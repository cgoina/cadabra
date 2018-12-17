%% Analyze_clusters
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
%% Cluster analysis based on transition matrices

% CLUSTER ANALYSIS
function FID = analyze_clusters(path,gen_ident,params,FID)
params.plots.clusteranalysis = params.plots.stat.dendro;

if params.plots.clusteranalysis && length(gen_ident)>1,
    % MERGE TRANSITIONS MATRICES?
    flyno = 0; % (0 = f1->f1 + f2->f2, 1 = f1, 2 = f2, 12 = f1->f2, 21 = f2->f1, 1221 = f1->f2 + f2->f1)
    addname = 'analysis/';
    params.trans_t = params.max_frames; %[s] transition time series, chunk size

    for igen=1:length(gen_ident),
        if iscell(path), path1 = path{igen}; else path1 = path; end
        genfname = gen_ident{igen}; 
        genfname = genfname(setdiff(1:length(genfname),strfind(genfname,'_')));
        ind = strfind(genfname,'/'); if numel(ind), genfname(ind) = '_'; end

        % LOAD TRANSITION ANALYSIS DATA
        fid = fopen([path1 addname genfname '_transitions.mat'],'r');
        if (fid < 0) || params.analyze_new,
            [h,lab] = transitions(path,gen_ident,igen,flyno,params,1);
        else
            load([path1 addname genfname '_transitions.mat']);
        end

        % MERGE TRANSITION MATRICES REGARDING 'flyno'
        [h2,lab2] = two2oneflymatr(h,lab,flyno);
        
        if ~exist('trans'),
            trans = zeros(length(gen_ident),numel(h2));
        end
        % Reshape transition data into a vector for each genotype
        trans(igen,:) = reshape(h2,1,numel(h2));
    end
    % MEASURE EUCLIDEAN DISTANCES
    C = pdist(trans,'euclidean');
    % figure(1); colormap('hot'); imagesc(squareform(C));
    % axis square; xlabel('sequence no.'); ylabel('sequence no.');
    % title('correlation-distance matrix'); colorbar;
    
    % COMPUTE AND PLOT DENDROGRAM
    L = linkage(C,'average');
    if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
    figure(FID); clf; [H,T,P] = dendrogram(L,'colorthreshold','default'); hold on;
    set(H,'LineWidth',2);
    gen_id_perm = gen_ident(P); set_xtick_label(gen_id_perm,45,params.cross_to);
    titl = 'Dendrogram from Transition Matrices';
    set(FID,'Name',titl); mtit(titl,'FontSize',params.axisfontsize+2,'yoff',.04);
    if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end  
end