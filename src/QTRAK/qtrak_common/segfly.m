%% SegFly
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
% * Original implementation by Heiko Dankert
% * Optimization and documentation by Edwin Soedarmadji
%
%%
% This function computes the threshold-based segmentation of an image.
% The main workhorse of this function is the Otsu's algorithm contained in
% the |otsu.m| file. 

function img_f01 = segfly(img_f01)

    %%
    % First, the image is adjusted by adding a bias level.
    % Next, the contrast is enhanced by taking the fourth power.
    % Image noise (i.e. due to compression) and edges are 
    % smoothed out by the convolution kernel.
    % Finally, segmentation is done using Otsu's algorithm
    % with par+1 segments requested.

    par = 4;                                
    kernel = [.25 .5 .25]';                 
    img_f01 = img_f01 + (abs(min(min(img_f01)))+1.2);
    img_f01 = img_f01.^4;
    img_f01 = conv2(img_f01,kernel,'same');
    
    img_f01 = fotsu(img_f01,par+1); 

    a = 0; 
    b = 0; 
    c = numel(img_f01); 
    i = par; 
    spar = 1/par/10;
    while (b<=0.8*c) && (a<=0.4*c) && (i>=0),
        a = numel(find(img_f01 >= i/par-spar & img_f01 <= i/par+spar)); 
        b = b + a; 
        i = i - 1;
    end

    %%
    % Perform morphological operations on the binary, thresholded version
    % of the segmented image. The effect of 'open' is to trim sharp edges. 
    % Next, connected components in the resulting binary image are labeled.
    % If there are more than two components (l01), the component with the
    % largest number of pixel is identified and designated as the "fly". 

    img_f01 = bwmorph(img_f01>(i+1)/par+spar,'open');
    [l11,l01] = bwlabel(img_f01);
    if l01 > 1, % remove noise
        m01 = zeros(1,l01);
        for i=1:l01, 
            m01(i) = numel(find(l11 == i)); 
        end
        ind = find(m01 == max(m01));
        img_f01(l11 ~= ind(1)) = 0;
    end

end