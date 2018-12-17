%% Seg_fullfly
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
% * Optimization and documentation by Edwin Soedarmadji 10/2008
%
%% Segmentation of Pixels of Full Fly (incl. body and wings)
% If possible, expand the pixel indices in several directions (1a),
% and get the enclosing rectangle (and its corners) of this area (1b).
% Finally, depending on how we get to this point, we calculate a
% segmentation for _f01_. If the flies started out nicely separated,
% we apply thresholding to turn _f01_ pixels into five segments.

function [f01,fly_ind1] = seg_fullfly(img2,f01,chamber,params,scale,ind_head)

% ...................(1a)
if ( max(f01) + chamber.s_arr(1) < chamber.s_arr(1) * chamber.s_arr(2)),
    f01 = [ f01; ...
        f01 - (params.cs-1) ; ...
        f01 + (params.cs-1) ; ...
        f01 - (params.cs-1) * chamber.s_arr(1) ; ...
        f01 + (params.cs-1) * chamber.s_arr(1) ; ...
        f01 - (params.cs-1) * (1+chamber.s_arr(1)) ; ...
        f01 + (params.cs-1) * (1+chamber.s_arr(1)) ];
end

% ...................(1b)
xa = min(chamber.y_arr(f01)/scale.y);
xb = max(chamber.y_arr(f01)/scale.y);
fly_ind1.r = uint16(xa : xb);
xa = min(chamber.x_arr(f01)/scale.x);
xb = max(chamber.x_arr(f01)/scale.x);
fly_ind1.c = uint16(xa : xb);

% ...................(1c)
if ind_head,
    tmp = img2(f01);
    f01 = unique(f01(tmp > 0.05));
else
    img_f01 = segfly(img2(fly_ind1.r,fly_ind1.c));
    xPixels = find( img_f01 > 0 );
    xDividend = floor( xPixels / length(fly_ind1.r) );
    xRemainder = mod( xPixels, length(fly_ind1.r) );
    f01 = ( xDividend  + double(min(fly_ind1.c)) - 1 )*size(img2,1) + ...
        ( xRemainder + double(min(fly_ind1.r)) - 1 );
end
