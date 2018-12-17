%% txtImage
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
%% Stamps a string 'txt' into a given raster image I, part of 'write_clip'
%% and 'plot_clips'

function I = txtImage(I,txt)
imtext = (1-text2im(txt));
imSize = size(I); imtextSize = size(imtext);
ratio = ceil(max(imtextSize ./ imSize));
if ratio > 1, 
    I = imresize(I,ratio+1,'bilinear');
    imSize = size(I);   
end
textImage = zeros(imSize);
textImage(1:imtextSize(1),1:imtextSize(2)) = imtext;
I(logical(textImage)) = round(max(max(I))/2);
if ratio > 1, I = imresize(I,2/(ratio+1),'bilinear'); end
