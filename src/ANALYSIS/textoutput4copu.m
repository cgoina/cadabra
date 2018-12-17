%% Textoutput
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
%% Output analysis data into a formated text file, i.e. for 
%% Excel import and further analysis

% OUTPUT ANALYSIS DATA INTO A FORMATED TEXT FILE
function textoutput4copu(data,gen_ident,titl,action_name,params,xtitl,xvals)

ind = strfind(params.PSFileN,params.slash); path = params.PSFileN(1:ind(end));
[s,mess,messid] = mkdir([path 'tables' params.slash]); %#ok<NASGU>
fname0 = action_name; in = strfind(fname0,' '); fname0 = fname0(setdiff(1:length(fname0),in));

% FIND MAXIMUM NUMBER OF FLY PAIRS PER GENOTYPE
ngens = size(data,1); pair_no = 0; maxlen = 0;
for igen=1:ngens,
    tmp = numel(data{igen,1}); namlen = length(gen_ident{igen});
    if tmp > pair_no, pair_no = tmp; end
    if namlen > maxlen, 
        maxlen = namlen; 
    end
end

if nargin < 7, 
    xtitl = 'pair no.'; 
    xvals = 1:pair_no;
end

% WRITE A FORMATED TABLE
% Columns = fly pairs
% Rows = genotypes and fly (winner/loser or 1/2)
fname = [path 'tables' params.slash fname0 '_' titl(~isspace(titl)) '.txt']; 
FID = fopen(fname,'w');
fprintf(FID,[action_name ' - ' titl '\n\n']);
if size(data,2) > 1
    fprintf(FID,['%' num2str(maxlen) 's       '],[xtitl ' -->']);
else
    fprintf(FID,['%' num2str(maxlen) 's'],[xtitl ' -->']);
end
if numel(strfind(titl,'number')),
    fprintf(FID,'%9g',xvals(1)); accur = '%9g'; 
else
    fprintf(FID,'%9g',xvals(1)); accur = '%9.3f'; 
end
fprintf(FID,'%9g',xvals(2:end)); fprintf(FID,'\n');

varname = {'Pre', 'Post'};
for igen=1:ngens,
    for ifly=1:size(data,2),
        if size(data,2) > 1,
            fprintf(FID,['%' num2str(maxlen) 's %5s'],cell2mat(gen_ident(igen)), varname{ifly});
        else
            fprintf(FID,['%' num2str(maxlen) 's'],cell2mat(gen_ident(igen)));
        end
        fprintf(FID,accur,data{igen,ifly}); fprintf(FID,'\n');
    end
end

fclose(FID);
