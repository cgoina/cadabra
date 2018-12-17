%% Auto_Cal
% Copyright (C) 2008 Heiko Dankert, California Institute of Technology

% This file is part of QTRAK
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
% * Heiko Dankert 12/2007
%
%%
% This function performs an autocalibration for the Heisenberg Chamber
% (single or double) by applying the Canny edge detector and morphological
% operations.

function coord = auto_cal(mean_image)
border = 50;
con = ones(5,5); con(:,2:4) = 0;
a = bwmorph(edge(ordfilt2(mean_image,5,con),'canny',0.05),'clean',3) + ...
    bwmorph(edge(ordfilt2(mean_image,5,con'),'canny',0.05),'clean',3);
a(1:border,:) = 0; a(size(a,1)-border:end,:) = 0; a(:,1:border) = 0; a(:,size(a,2)-border:end) = 0;

inda = find(mean(a(size(a,1)/2-10:size(a,1)/2+10,:),1)>.5);
inda = inda([inda(1) inda(2:end)-inda(1:end-1)]>2);
b = mean(mean_image(size(mean_image,1)/2-10:size(mean_image,1)/2+10,:));
in = get_ind(b,inda,border);
coord.x = inda(in);

coord = 0;
if ~mod(numel(in),2),
    coord.x = inda(in);
    b = mean(mean_image(:,(coord.x(1)+coord.x(2))/2-10:(coord.x(1)+coord.x(2))/2+10),2)';
    inda = find(mean(a(:,(coord.x(1)+coord.x(2))/2-10:(coord.x(1)+coord.x(2))/2+10),2)>.5)';
    if numel(inda),
        in = get_ind(b,inda,border);
        coord.y = inda(in);
    else
        coord = 0;
    end
end
end

function in = get_ind(b,inda,border)
b(1:border) = max(b); b(length(b)-border:end) = max(b);
indb = find(b < 0.2*(max(b(border:end-border))-min(b(border:end-border)))+min(b(border:end-border)));
d = zeros(1,numel(inda)); for i=1:numel(inda), d(i) = min(abs(indb-inda(i))); end
dd = [inda(1) inda(2:end) - inda(1:end-1)];
in = []; icnt = min(d);
while numel(in) < 4 && icnt < min(d)+30,
    in = find(d <= icnt & dd > 5);
    icnt = icnt + 1;
end
end