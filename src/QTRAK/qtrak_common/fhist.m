%% FHist
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
% * Written by Edwin Soedarmadji
%
%% 
% This histogram function is optimized specifically for use with QTracker. 
% As much as possible, it is programmed with compatibility to Matlab's 
% regular histogram function _hist_. 

function [no,xo] = fhist(y, x)

%% 
% Determine the edges from the minimum and maximum data values if the 
% bin vector is not specified (and instead only the number of bins is 
% specified by the user. Otherwise, calculate the edges from the bin 
% vector itself. 

if length(x) == 1
    ind = ~isnan(y);
    miny = min(y(ind));
    maxy = max(y(ind));
    if (isempty(miny))
      miny = NaN;
      maxy = NaN;
    end
    if miny == maxy,
      miny = miny - floor(x/2) - 0.5; 
      maxy = maxy + ceil(x/2) - 0.5;
    end
    binwidth = (maxy - miny) ./ x;
    xx = miny + binwidth*(0:x);
    xx(length(xx)) = maxy;
    x = xx(1:length(xx)-1) + binwidth/2;
else
    xx = x(:)';
    binwidth = [diff(xx) 0];
    xx = [xx(1)-binwidth(1)/2 xx+binwidth/2];
end
  
%%
% Shift bins so the interval is ( ] instead of [ ).

xx = full(real(xx)); 
y = full(real(y)); 
bins = xx + eps(xx);

%%
% Call the optimized histogram MEX-function

nn = zeros( 1, length(bins) + 1 );
fhistc( y,[-inf bins],nn );
    
%% 
% Combine the first bin with the 2nd bin and 
% the last bin with the next-to-last bin

nn(2) = nn(2)+nn(1);
nn(end-1) = nn(end-1)+nn(end);
no = nn(2:end-1);
xo = x;
