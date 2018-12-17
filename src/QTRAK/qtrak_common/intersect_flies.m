%% Intersect_flies
% Copyright (C) 2008 Heiko Dankert, California Institute of Technology

% This file is part of QTRAK.

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
% * Written by Heiko Dankert 2006/2007
% * Documentation by Edwin Soedarmadji 10/2008
%
%% Resolve Fly Intersections
% Having obtained the pixel indices _f01_ and _f02_ for flies 1
% and 2, the code in (16a) calculates the overlap between _f01_ and
% _f02_. The function |[c, if01, if02] = intersect(f01, f02)| below
% calculates three pixel indices:
% * _c_ the intersection pixels
% * _if01_, the pixel indices of f01 such that |c = f01(if01)|
% * _if01_, the pixel indices of f01 such that |c = f01(if01)|
%
% If the intersection is not empty (16b), then mark those
% intersecting pixels in _f01_ and _f02_ with -99, and extract into
% _f01_ and _f02_ only those pixels that are not marked. The two
% objects are then merged with the body pixels _f011_ and _f022_
% (16c). Finally, we do one more test to determine if _f022_ meets
% the conditions to be considered a fly.

function [f01,f02] = intersect_flies(f01,f02,f011,f022)

% ...................(16a)
[c,if01,if02] = intersect(f01,f02);
lc = numel(c);

% ...................(16b)
if lc,
    lf01 = numel(f01);
    lf02 = numel(f02);
    if (lf01 > lc),
        f01(if01) = -99;
        f01 = f01(f01 > -99);
    end
    if (lf02 > lc),
        f02(if02) = -99;
        f02 = f02(f02 > -99);
    end
end

% ...................(16c)
f01 = union(f01,f011);
f02 = union(f02,f022);
end