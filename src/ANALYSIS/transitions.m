%% Transitions
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
%% Compute transition matrix

% TRANSITION MATRIX
function [h,lab,chains] = transitions(path,gen_ident,igen,flyno,params,iloop)

% LOAD DATA
addname = 'analysis/';
if iscell(path), path1 = path{igen}; else path1 = path; end
genfname = gen_ident{igen};
genfname = genfname(setdiff(1:length(genfname),strfind(genfname,'_')));
ind = strfind(genfname,'/'); if numel(ind), genfname(ind) = '_'; end
load([path1 addname genfname '_aggress_court.mat']);

% LABELS OF ACTIONS FOR THE TRANSITION MATRIX PLOT
lab = {'lunge f1', 'lunge f2', 'tussle f1f2', 'wingthr f1', 'wingthr f2', 'chase f1', 'chase f2', ...
       'wingext f1', 'wingext f2', 'circl f1', 'circl f2', 'circw f1', 'circw f2', 'nobeh f1f2'};

nmov = max(tussl.mov);
num_beh = length(lab);
h = zeros(num_beh,num_beh); ichain = 0; chains = [];
if  flyno > 2, nfl = 1; else nfl = 2; end
for imov=1:nmov,
    for ifly=1:nfl,
        % COLLECT ACTION START & END TIME (X) AND ACTION CODE (I) FOR EACH
        % MOVIE AND FLY; EACH ACTION AND FLY HAS A DIFFERENT CODE
        I = []; X.A = []; X.E = [];
        if ifly == 1 || flyno > 2,
            if cop.pre(imov) == 0, [X,I] = transm_prep(lunge.obj1.t, lunge.obj1.number, lunge.obj1.t*0+1, cop, dt(imov), 1, X, I, imov, params.trans_t, iloop); end
            [X,I] = transm_prep(wing.threat.obj1.t, wing.threat.obj1.number, wing.threat.obj1.len, cop, dt(imov), 4, X, I, imov, params.trans_t, iloop);
        end
        if ifly == 2 || flyno > 2,
            if cop.pre(imov) == 0, [X,I] = transm_prep(lunge.obj2.t, lunge.obj2.number, lunge.obj2.t*0+1, cop, dt(imov), 2, X, I, imov, params.trans_t, iloop); end
            [X,I] = transm_prep(wing.threat.obj2.t, wing.threat.obj2.number, wing.threat.obj2.len, cop, dt(imov), 5, X, I, imov, params.trans_t, iloop);
        end
        if cop.pre(imov) == 0, [X,I] = transm_prep(tussl.t, tussl.number, tussl.len, cop, dt(imov), 3, X, I, imov, params.trans_t, iloop); end
        if ifly == 1 || flyno > 2,
            [X,I] = transm_prep(chase.obj1.t, chase.obj1.number, chase.obj1.len, cop, dt(imov), 6, X, I, imov, params.trans_t, iloop);
            [X,I] = transm_prep(wing.ext.l.obj1.t, wing.ext.l.obj1.number, wing.ext.l.obj1.len, cop, dt(imov), 8, X, I, imov, params.trans_t, iloop);
            [X,I] = transm_prep(wing.ext.r.obj1.t, wing.ext.r.obj1.number, wing.ext.r.obj1.len, cop, dt(imov), 8, X, I, imov, params.trans_t, iloop);
            [X,I] = transm_prep(court.obj1.t, court.obj1.number, court.obj1.t, cop, dt(imov), 10, X, I, imov, params.trans_t, iloop);
        end
        if ifly == 2 || flyno > 2,
            [X,I] = transm_prep(chase.obj2.t, chase.obj2.number, chase.obj2.len, cop, dt(imov), 7, X, I, imov, params.trans_t, iloop);
            [X,I] = transm_prep(wing.ext.l.obj2.t, wing.ext.l.obj2.number, wing.ext.l.obj2.len, cop, dt(imov), 9, X, I, imov, params.trans_t, iloop);
            [X,I] = transm_prep(wing.ext.r.obj2.t, wing.ext.r.obj2.number, wing.ext.r.obj2.len, cop, dt(imov), 9, X, I, imov, params.trans_t, iloop);
            [X,I] = transm_prep(court.obj2.t, court.obj2.number, court.obj2.t, cop, dt(imov), 11, X, I, imov, params.trans_t, iloop);
        end
        % SORT (X,I) BY ACTION START TIMES (X.A)
        [X.A,ix] = sort(X.A); X.E = X.E(ix); I = I(ix)';
        I0 = I;
        NX = length(X.A);
        if NX > 1,
            % CHECK TIME DIFFERENCE BETWEEN START TIME OF ONE ACTION AND
            % END TIME OF PREVIOUS ACTION (D), START TIME OF ONE ACTION AND
            % START TIME OF PREVIOUS ACTION (DAA), END TIME OF ONE ACTION AND
            % END TIME OF PREVIOUS ACTION (DEE), 
            D = X.A(2:end) - X.E(1:end-1);
            DAA = X.A(2:end) - X.A(1:end-1);
            DEE = X.E(2:end) - X.E(1:end-1);

            % FIND TRANSITIONS INTO COMBINED WING EXTENSION + CIRCLING
            ind = find(DAA > 0 & D < 0 & DAA <= params.trans_thr);
            in = ind((I(ind) == 8 & I(ind+1) == 10) | (I(ind) == 10 & I(ind+1) == 8));
            if numel(in), D(in) = 2; I(in+1) = 12; end
            in = ind((I(ind) == 9 & I(ind+1) == 11) | (I(ind) == 11 & I(ind+1) == 9));
            if numel(in), D(in) = 2; I(in+1) = 13; end

            % COMPUTE TRANSITION MATRIX
            ind = find(D > 0 & D <= params.trans_thr);
            for i=1:numel(ind),
                h(I(ind(i)),I(ind(i)+1)) = h(I(ind(i)),I(ind(i)+1)) + 1;
                % TREAT SPECIAL CASES
                % A new action starts before an ongoing action is finished
                if DAA(ind(i)) > 0 && DEE(ind(i)) < 0 && I(ind(i)) > 7 && ind(i)+2 <= numel(ind),
                    h(I(ind(i)+1),I(ind(i))) = h(I(ind(i)+1),I(ind(i))) + 1;
                    I(ind(i)+1) = I(ind(i));
                end
                % A combined wing extension + circling ends partially,
                % resulting in pure wing extension or circling
                if D(ind(i)) == 2 && DAA(ind(i)) > 0 && DEE(ind(i)) > 0 && I(ind(i)) > 7,
                    h(I(ind(i)+1),I0(ind(i)+1)) = h(I(ind(i)+1),I0(ind(i)+1)) + 1;
                    I(ind(i)+1) = I0(ind(i)+1);
                end
                % COMPUTE 3-DIGIT CODE FOR CHAINS OF 3 ACTIONS               
                if i > 1,
                    ichain = ichain + 1;
                    chains.i(ichain) = str2double([num2str(I(ind(i)-1),'%02g') num2str(I(ind(i)),'%02g') num2str(I(ind(i)+1),'%02g')]);
                    chains.t(ichain) = X.A(ind(i));
                    chains.f(ichain) = ifly;
                    chains.m(ichain) = imov;
                end
            end
            % TRANSITION OF ACTIONS INTO "NO ACTION"
            ind = find(D > params.trans_thr);
            for i=1:numel(ind),
                h(I(ind(i)),num_beh) = h(I(ind(i)),num_beh) + 1;
                h(num_beh,I(ind(i)+1)) = h(num_beh,I(ind(i)+1)) + 1;
            end
        end
    end
end
% NORMALIZE THE NUMBER OF TRANSITIONS BY 
% NUMBER OF FLY PAIRS/MOVIES
h = h / nmov;

% SAVE TRANSITION MATRIX AND LABELS
save([path1 addname genfname '_transitions.mat'],'h','lab','chains');
