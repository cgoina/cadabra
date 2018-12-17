%% Extract_wing_orient
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
%% Extract relativ orientation and position of one fly towards 
%% the other fly, while the other fly is extending a wing

function [ob1,ob2] = extract_wing_orient(igen,wings,extfli,lrb,add,shift,phi0,genotypes,sorted_names,copu)
    indc1 = eval(['cumsum(wings.' extfli '.' lrb '{' num2str(igen) '}.obj1.number)']);
    indm1 = eval(['wings.' extfli '.' lrb '{' num2str(igen) '}.obj1.number']);
    indc1 = [0 indc1(indm1>0)];
    indc2 = eval(['cumsum(wings.' extfli '.' lrb '{' num2str(igen) '}.obj2.number)']);
    indm2 = eval(['wings.' extfli '.' lrb '{' num2str(igen) '}.obj2.number']);
    indc2 = [0 indc2(indm2>0)];
    ob1.tim = cell(1,numel(indm1)); ob1.do1 = ob1.tim; ob1.do2 = ob1.tim; ob1.dis = ob1.tim;
    ob1.tim = cell(1,numel(indm1)); ob1.do1 = ob1.tim; ob1.do2 = ob1.tim; ob1.dis = ob1.tim;
    ob2 = ob1;
    gen_ind = find(genotypes == igen);    
    for ifile = 1:numel(indm1),
        if indm1(ifile) || indm2(ifile),
            if iscell(sorted_names), FileN = sorted_names{gen_ind(ifile)}; else FileN = sorted_names; end
%             fprintf('-----------------------------------');
%             fprintf('\n'); fprintf('%s',FileN); fprintf('\n');
            fid = fopen([FileN(1:end-5) '_feat.mat'],'r');
            if (fid < 0),
                tic;
                [fly_feat,NFrms] = read_feat(FileN);
                save([FileN(1:end-5) '_feat.mat'],'fly_feat','NFrms');
                ttt=toc; fprintf(1,'left %2.0f:%2.0f mm:ss\n',floor(ttt * (length(gen_ind)-ifile+1)/60),mod(ttt * (length(gen_ind)-ifile+1),60));
            else
                load([FileN(1:end-5) '_feat.mat']);
                fclose(fid);
            end
            if copu.pre{igen}(ifile) > 1,
                tmp = fly_feat(1).time(2:NFrms(1)) - fly_feat(1).time(1:NFrms(1)-1);
                dt = tmp(10);            
                copu_start = copu.pre{igen}(ifile) * dt;
            else
                copu_start = fly_feat.time(end);
            end
            if indm1(ifile),
                icnt1 = sum(indm1(1:ifile)>0);
                ob1.tim{ifile} = cell(1,indm1(ifile));
                ob1.do1{ifile} = ob1.tim{ifile}; ob1.do2{ifile} = ob1.tim{ifile};
                ob1.dis{ifile} = ob1.tim{ifile}; 
                ob1.x1{ifile} = ob1.tim{ifile}; ob1.y1{ifile} = ob1.tim{ifile};
                ob1.x2{ifile} = ob1.tim{ifile}; ob1.y2{ifile} = ob1.tim{ifile};
                for j=1:indm1(ifile),
                    st = eval(['sum(wings.' extfli '.' lrb '{' num2str(igen) '}.obj1' ...
                        '.len(1:indc1(icnt1)+j-1))+1']);
                    en = eval(['sum(wings.' extfli '.' lrb '{' num2str(igen) '}.obj1' ...
                        '.len(1:indc1(icnt1)+j))']);
                    lab = eval(['wings.' extfli '.' lrb '{' num2str(igen) '}.obj1.lab(st:en)']);
                    phi = eval(['wings.' extfli '.' lrb '{' num2str(igen) '}.obj1' ...
                        '.phi' lrb '(indc1(icnt1)+j)']);
                    if (phi >= phi0.min) && (phi < phi0.max) && (fly_feat.time(lab(end)) < copu_start),
                        % Extract period before, after, or during wing extension
                        if shift<0,
                            lab = 2*lab(1)-lab(end)-1:lab(1)-1;
                        elseif shift>0,
                            lab = lab(end)+1:2*lab(end)-lab(1)+1;
                        else
                            lab = lab(1)-add:lab(end)+add;
                        end
                        lab = lab(lab>0 & lab<=NFrms);
                        % Extract: time of wing extension of fly 1 (tim), 
                        % head direction of fly 1,2
                        % (do1, do2), distance between flies (dis), 
                        % position of fly 1, 2 (x,y)
                        ob1.tim{ifile}{j} = fly_feat.time(lab);
                        ob1.do1{ifile}{j} = fly_feat.obj1.headdir(lab);
                        ob1.do2{ifile}{j} = fly_feat.obj2.headdir(lab);
                        ob1.dis{ifile}{j} = fly_feat.distc(lab);
                        ob1.x1{ifile}{j} = fly_feat.obj1.pos_x(lab);
                        ob1.y1{ifile}{j} = fly_feat.obj1.pos_y(lab);
                        ob1.x2{ifile}{j} = fly_feat.obj2.pos_x(lab);
                        ob1.y2{ifile}{j} = fly_feat.obj2.pos_y(lab);
                    end
                end
            end
            if indm2(ifile), 
                icnt2 = sum(indm2(1:ifile)>0);
                ob2.tim{ifile} = cell(1,indm2(ifile));
                ob2.do1{ifile} = ob2.tim{ifile}; ob2.do2{ifile} = ob2.tim{ifile};
                ob2.dis{ifile} = ob2.tim{ifile};
                ob2.x1{ifile} = ob2.tim{ifile}; ob2.y1{ifile} = ob2.tim{ifile};
                ob2.x2{ifile} = ob2.tim{ifile}; ob2.y2{ifile} = ob2.tim{ifile};
                for j=1:indm2(ifile),
                    st = eval(['sum(wings.' extfli '.' lrb '{' num2str(igen) '}.obj2' ...
                        '.len(1:indc2(icnt2)+j-1))+1']);
                    en = eval(['sum(wings.' extfli '.' lrb '{' num2str(igen) '}.obj2' ...
                        '.len(1:indc2(icnt2)+j))']);
                    lab = eval(['wings.' extfli '.' lrb '{' num2str(igen) '}.obj2.lab(st:en)']);
                    phi = eval(['wings.' extfli '.' lrb '{' num2str(igen) '}.obj2' ...
                        '.phi' lrb '(indc2(icnt2)+j)']);
                    if (phi >= phi0.min) && (phi < phi0.max),
                        if shift<0,
                            lab = 2*lab(1)-lab(end)-1:lab(1)-1;
                        elseif shift>0,
                            lab = lab(end)+1:2*lab(end)-lab(1)+1;
                        else
                            lab = lab(1)-add:lab(end)+add;
                        end
                        lab = lab(lab>0 & lab<=NFrms);
                        % Extract: time of wing extension of fly 2 (tim), 
                        % head direction of fly 1,2
                        % (do1, do2), distance between flies (dis), 
                        % position of fly 1, 2 (x,y)
                        ob2.tim{ifile}{j} = fly_feat.time(lab);
                        ob2.do1{ifile}{j} = fly_feat.obj2.headdir(lab);
                        ob2.do2{ifile}{j} = fly_feat.obj1.headdir(lab);
                        ob2.dis{ifile}{j} = fly_feat.distc(lab);
                        ob2.x1{ifile}{j} = fly_feat.obj2.pos_x(lab);
                        ob2.y1{ifile}{j} = fly_feat.obj2.pos_y(lab);
                        ob2.x2{ifile}{j} = fly_feat.obj1.pos_x(lab);
                        ob2.y2{ifile}{j} = fly_feat.obj1.pos_y(lab);
                    end
                end
            end
        end
    end
end