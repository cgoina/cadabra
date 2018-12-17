%% FHistFit
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
% This function computes the histogram of _inimg_, and returns two images:
% a binary (logical) image _binimg_ containing ROUGHLY the suspected fly 
% objects and a copy of _inimg_ masked by _binimg_ (called _outimg_). 
% Unlike its counterpart _histfit_, this function uses the optimized 
% histograming function _fhist_. The function _histfit_ has a different 
% set of return values, and thus retained for use in other parts of the code. 

function [ outimg, binimg ] = fhistfit( inimg )
global CurrentK;
persistent tmp;

    %%
    % Start by initializing the bins, iteration variables, and dimensions
    
    BINS   = -0.10 : 0.02 : 1.00;
    BIN_SZ = BINS(2)-BINS(1);

    %% 
    % Next, use the bins to produce a histogram (tmp)
    % This histogram is then smoothed by the kernel smo.
    % Note that _fhistfit_ uses the optimized _fhist_ function. 
    
    if isempty(tmp) || mod(CurrentK,10) == 0,
        tmp = fhist( inimg(:), BINS );
    end
    smo = [1 2 1];
    smo = smo ./ sum(smo);
    hist1 = conv2(tmp,smo,'same');

    %%
    % The left hand side of the expression below is just the maximum 
    % center coordinate of the non-empty bins. If this coordinate is 
    % less than or equal to 0.25, the fly most likely cannot be detected. 
    
    if max( BINS( hist1 > 0 ) ) <= 0.25,  
        outimg = inimg .* 0;
        binimg = outimg > 0;
    else
      
        %%
        % The background image is hereby defined as the pixel population 
        % with values less than the histogram maxima + 2 bins. 
        % If the peak population meets a certain condition, refine
        % the bins to accomodate analysis of non-background pixels. 
        
        ind_max = find(hist1 == max(hist1));
        if ((1/BIN_SZ)*hist1(ind_max)/sum(hist1) > 13),
            BINS = BIN_SZ * (ind_max+2) : 0.01 : 1;
        end

        %%
        % The bin value at (ind_cross - shift) is the thresholding value 
        % used to mask the input image.
        
        binimg = inimg > BINS(1);
        outimg = inimg .* binimg;

    end

end