%% ProcessFrameMeanWin
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
% * Original implementation by Heiko Dankert
% * Optimization and documentation by Edwin Soedarmadji
%
%%
% Just like _ProcessFrameCapWin_, this Matlab function is called from the 
% frame grabber. The input arguments are identical, although the variable 
% _image_ is really undefined because when this function is executed, the 
% frame grabber has not been supplied with the mean image, and thus the 
% chamber pixels are not computed. 

function ProcessFrameMeanWin(data, images, width, height, frameNr, time) %#ok<INUSD,INUSL>
global mean_image std_image nframes_mean frmindx;
persistent h1;

    try
        %% 
        % If it's the first frame, then initialize the progress bar 
        % control, as well as the mean and std images with zeros. 
        % If it's the last frame in this mean image calculation, close 
        % the progress bar control, and otherwise, update the mean and 
        % std images and the progress bar control.
        
        frmindx = frmindx + 1;
        if (frmindx == 1)
            h1 = waitbar(0,'Computing background image');
            mean_image = zeros(height,width,3);
            std_image = zeros(height,width,3);
%         elseif (frmindx >= nframes_mean)
%             close(h1);
        elseif (frmindx >= 1)
            mean_image = mean_image + data;
            std_image = std_image + data.^2;
            waitbar(frmindx/nframes_mean,h1);
        end
    catch err
        err.disp;
    end

end