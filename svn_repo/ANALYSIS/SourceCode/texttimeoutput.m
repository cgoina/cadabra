%% Texttimeoutput
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
%% Output the action occurrence dates into a formated 
%% text file, i.e. for Excel import and further analysis

% OUTPUT ACTION OCCURRENCE DATES INTO A FORMATED TEXT FILE
% Dates are relative to the start time point of a movie
function texttimeoutput(behavior,gen_ident,action_name,params)

ngens = length(behavior);
ind = strfind(params.PSFileN,params.slash); path = params.PSFileN(1:ind(end));
[s,mess,messid] = mkdir([path 'tables' params.slash]); %#ok<NASGU>
fname0 = action_name; in = strfind(fname0,' '); fname0 = fname0(setdiff(1:length(fname0),in));

for igen=1:ngens,
    npairs = numel(behavior{igen}.obj1.number); 
    maxnum = max(max(behavior{igen}.obj1.number),max(behavior{igen}.obj2.number));
    fname = [path 'tables' params.slash fname0 '_dates_' ...
             cell2mat(gen_ident(igen)) '.txt'];
    FID = fopen(fname,'w');
    fprintf(FID,[action_name ' - ' cell2mat(gen_ident(igen)) '\n']);
    fprintf(FID,['Dates in [s] from movie start\n\n']);
    fprintf(FID,'%12s','pair no. -->');
    fprintf(FID,'%11g',1); fprintf(FID,'%19g',2:npairs); 
    fprintf(FID,'\n');
    fprintf(FID,'%12s','fly -->');
%     if params.winlos,
%         onetwos = repmat({'win' 'los'},1,npairs);
%     else
        onetwos = repmat({'1 ' '2 '},1,npairs); 
%     end
    fprintf(FID,'%7s',onetwos{1});
    donetwos = repmat({'%9s' '%10s'},1,npairs);
    for ii=2:length(onetwos),
        fprintf(FID,donetwos{ii},onetwos{ii});
    end
    fprintf(FID,'\n');
    
    csum1 = [0 cumsum(behavior{igen}.obj1.number)];
    csum2 = [0 cumsum(behavior{igen}.obj2.number)];
    data = zeros(maxnum,2*npairs);
    for ipair=1:npairs,
        st = csum1(ipair)+1; en = csum1(ipair+1);
        if (en-st)>0, 
            dat1 = behavior{igen}.obj1.t(st:en);
            data(1:numel(dat1),2*ipair-1) = dat1;
        end        
        st = csum2(ipair)+1; en = csum2(ipair+1);
        if (en-st)>0, 
            dat2 = behavior{igen}.obj2.t(st:en);
            data(1:numel(dat2),2*ipair) = dat2;
        end
    end
    
    for ievent=1:maxnum,
        fprintf(FID,'%12s',' ');
        for ipair=1:npairs,
            if data(ievent,2*ipair-1), 
                dat1 = num2str(data(ievent,2*ipair-1),'%9.3f');
            else
                dat1 = '         ';
            end
            if data(ievent,2*ipair), 
                dat2 = num2str(data(ievent,2*ipair),'%9.3f');
            else
                dat2 = '         ';
            end
            fprintf(FID,'%9s %9s',dat1,dat2);
        end
            fprintf(FID,'\n');
    end
    fclose(FID);
end
