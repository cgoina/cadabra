%% Ethograms
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
%% Ethogram main routine

% ETHOGRAMS
function FID = ethograms(path,gen_ident,gen_ident1,params,FID)

if params.plots.ethogram,
    for ilen=1:numel(params.plots.ethogram_vec),
        if (params.plots.ethogram_vec(ilen) == 2) || (params.plots.ethogram_vec(ilen) == 4) || (params.plots.ethogram_vec(ilen) == 6),
            % WHICH REACTIONS?
            flyno = 1221; % (12 = f1->f2, 21 = f2->f1, 1221 = f1->f2 + f2->f1)
        else
            % MERGE TRANSITIONS OF FLIES?
            flyno = 0; % (0 = f1->f1 + f2->f2, 1 = f1, 2 = f2)
        end
        
        % PLOT TRANSITION MATRICES?
        if (params.plots.ethogram_vec(ilen) == 5) || (params.plots.ethogram_vec(ilen) == 6),
            bool_trans = 1;
        else
            bool_trans = 0;
        end

        % TRANSITION TIME SERIES OR GLOBAL?
        if (params.plots.ethogram_vec(ilen) == 3) || (params.plots.ethogram_vec(ilen) == 4),
            % Transition time series;  split movie into chunks of 'trans_t' seconds
            params.trans_t = 300; % [s] 
            nloops = ceil(params.max_frames / params.trans_t);
        else
            % Compute all transitions over full movie length
            params.trans_t = params.max_frames;            
            nloops = 1;
        end
        
        % OCCURENCE SCALING FOR TRANSITION MATRICES
        max_v = 20;
        % NUMBER OF PLOTS PER PAGE
        nplt = 4;

        k1 = floor(sqrt(nplt)); k2 = ceil(nplt/k1);
        icnt = 1;
        for igen=1:length(gen_ident),
            for iloop=1:nloops,
                if ~mod(icnt-1,4),
                    if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
                    icnt = 1;
                    figure(FID); clf;
                    set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
                end
                if iscell(path), path1 = path{igen}; else path1 = path; end
                genfname = gen_ident{igen};
                genfname = genfname(setdiff(1:length(genfname),strfind(genfname,'_')));
                ind = strfind(genfname,'/'); if numel(ind), genfname(ind) = '_'; end

                % ANALYZE TRANSITIONS/REACTIONS
                [h,lab,chains] = transitions(path,gen_ident,igen,flyno,params,iloop);
                
                % CHAIN ANALYSIS
%                 chain_analysis(chains);

                % MERGING OF TRANSITION MATRICES
                [h2,lab2] = two2oneflymatr(h,lab,flyno);
                
                % PLOTS
                figure(FID); subplot(k1,k2,icnt);

                % PLOT TRANSITION MATRIX OR ETHOGRAM?
                if bool_trans,
                    etho_matrix_plot(h2,lab2,max_v,8);
                else
                    ethogram_plot(h2,lab2);
                end

                tit = gen_ident1{igen}; in = strfind(tit,'_'); tit(in) = '-';
                if nloops == 1,
                    title(tit);
                else
                    title([tit ', ' num2str((iloop-1)*params.trans_t) '-' ...
                           num2str(iloop*params.trans_t) ' s']);
                end
                
                if params.pdf,
                    if ~mod(icnt,nplt) || (igen == length(gen_ident) && iloop == nloops),
                        if ~flyno,
                            titl = 'f1->f1 + f2->f2';
                        elseif flyno == 1,
                            titl = 'f1->f1';
                        elseif flyno == 2,
                            titl = 'f2->f2';
                        elseif flyno == 12,
                            titl = 'f1->f2';
                        elseif flyno == 21,
                            titl = 'f2->f1';
                        elseif flyno == 1221,
                            titl = 'f1->f2 + f2->f1';
                        end
                        if flyno >= 12,
                            titl = ['Reactions, ' titl];
                        else
                            titl = ['Transitions, ' titl];
                        end
                        if ~bool_trans, titl = ['Ethograms, ' titl]; end

                        titl = [titl ', t <= ' num2str(params.trans_thr) ' s'];
                        set(FID,'Name',titl); mtit(titl,'FontSize',16','yoff',.02);
                        print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN);
                    end
                end
                icnt = icnt + 1;
            end
        end
    end
end