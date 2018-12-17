%% Utestpairs
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
%% Find genotype pairs to be tested pairwise by Mann-Whitney U-test in 
%% 'plot_utests'

% GENOTYPE PAIRS TO BE STATISTICALLY COMPARED
function utest = utestpairs(gen_ident)
ngens = numel(gen_ident);
% FIND GENOTYPE AND "OVER", i.e. "genotype/kir"
% This part may be made more universal
gen = cell(ngens,1); over = gen;
for igen=1:ngens,
    ind = strfind(gen_ident{igen},'/');
    if numel(ind), 
        gen{igen} = gen_ident{igen}(1:ind(1)-1); 
        over{igen} = gen_ident{igen}(ind(1)+1:end); 
    else
        gen{igen} = gen_ident{igen}; 
        over{igen} = 'kir'; 
    end
end
% FIND CONTROL LINE
utest.ind_kirop = find(strcmp(gen_ident,'kir/+'));
if ~numel(utest.ind_kirop), utest.ind_kirop = find(strcmp(gen_ident,'w-Exel')); end
if ~numel(utest.ind_kirop), utest.ind_kirop = find(strcmp(gen_ident,'w-exel')); end
if ~numel(utest.ind_kirop), utest.ind_kirop = find(strcmp(gen_ident,'CS-Tully')); end
if ~numel(utest.ind_kirop), utest.ind_kirop = find(strcmp(gen_ident,'CS')); end
utest.ind_over = setdiff(find(strcmp(over,'kir')),utest.ind_kirop);
% COMPARE EACH GENOTYPE WITH THE REST
ugen = unique(gen); utest.ind_utest = [];
for iugen=1:numel(ugen),
    ind = find(strcmp(gen,ugen{iugen}));
    if numel(ind) > 1,
        for i=1:numel(ind)-1,
            for ii=i+1:numel(ind),
                utest.ind_utest = [utest.ind_utest ; ind(i) ind(ii)];
            end
        end
    end
end

end