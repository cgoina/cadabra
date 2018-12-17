%% Gen_name
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
%% Automatic filter for genotype names, encoded in the file name

% FILE NAME CONVENTION:
% "YYMMDD_Lgenotpe_Rgenotype_SX"; genotype Left/Right chamber; X = assay number
% "YYMMDD_genotype_SX" for same-genotype or single-chamber movies
% do not use "_" within genotype names, use "-" instead
% THERE IS NO NAME DECODING FOR ASSAYS WITH TWO DIFFERENT GENOTYPES 

% DECODE GENOTYPE NAME FROM FILE NAME
function gen = genname(fname,n)

ind = strfind(fname,'_');
if ~isempty(ind),
    strt = ind(1)+1;
    if numel(ind) > 1,
        strt2 = ind(end)-1;
        L = strfind(fname,'L'); R = strfind(fname,'R');
        if numel(L) && numel(R),
            L = L(1); R = R(1);
            if (n == 1),
                gen = fname(L+1:R-1);
            elseif (n == 2),
                l = strfind(fname(R:end),'_'); if ~isempty(l), l = R + l(1)-1; else l = length(fname); end
                gen = fname(R+1:l);
            end
        else
            gen = fname(strt:strt2);
        end
%         if length(gen) > 1,
%             ind = strfind(gen,'_');
%             if numel(ind), gen = gen(1:ind(1)-1); end
% %             if strcmp(gen(end-1),'_'), gen = gen(1:end-2); end
% %             if strcmp(gen(end),'_'),
% %                 gen = gen(1:end-1);
% %             end
%         end
    else
        gen = fname(strt:end);
    end
else
    gen = fname;
end
