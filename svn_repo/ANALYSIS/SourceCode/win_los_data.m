%% Win_los_data
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
%% Process detected actions using provided operator and find out the
%% loser/winner for each pair

% SORT ACTION DATA INTO WINNER/LOSER OR FLY 1/2 AND PEFORM OPERATION
function data = win_los_data(field,operat,obj1,obj2,lunges,pre_copu,dts,params)
% field: .. to be processed
% operat: operation such as max, mean, median

ngens = length(obj2);
data = cell(ngens,1);
for igen=1:ngens,
    % FIND START/END INDICES OF MOVIES/FLY PAIRS
    if ~strcmp(field,'movdists'),
        ind = [0 find(obj1{igen}.t(2:end)-obj1{igen}.t(1:end-1)<0) length(obj1{igen}.t)];
    else
        ind = [0 find(obj1{igen}.movdists(2:end)-obj1{igen}.movdists(1:end-1)<0) length(obj1{igen}.movdists)];
    end
    if ~ind, ind = [ind length(obj1{igen}.t)]; end
    d.win = []; d.los = [];
    for j=1:numel(ind)-1,
        % FOR EACH MOVIE/FLY PAIR SORT ACTION DATA
        tim = obj1{igen}.t(ind(j)+1:ind(j+1));
        if (lunges{igen}.obj1.number(j) >= lunges{igen}.obj2.number(j)) || ~params.winlos,
            dat1 = eval(['obj1{' num2str(igen) '}.' field '(ind(j)+1:ind(j+1))']);
            dat2 = eval(['obj2{' num2str(igen) '}.' field '(ind(j)+1:ind(j+1))']);
        else
            dat1 = eval(['obj2{' num2str(igen) '}.' field '(ind(j)+1:ind(j+1))']);
            dat2 = eval(['obj1{' num2str(igen) '}.' field '(ind(j)+1:ind(j+1))']);
        end
        % In case of copulation take only time span berfore copulation
        % start into account
        limt = pre_copu{igen}(j)*dts{igen}(j);
        indt = find(tim < limt);
        if numel(indt),
            dat1 = dat1(indt); dat2 = dat2(indt);
        end
        % PERFORM OPERATION ON DATA
        d.win = [d.win eval([operat '(dat1)'])]; d.los = [d.los eval([operat '(dat2)'])];
    end
    data{igen} = d; 
end

end