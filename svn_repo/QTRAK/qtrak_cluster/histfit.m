%% HistFit
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
% a binary (logical) image _binimg_ containing the suspected fly objects,  
% and a copy of _inimg_ masked by _binimg_ (called _outimg_). Unlike its 
% counterpart _fhistfit_, this function uses the standard histograming 
% function _hist_, and computes the mean and variance vector _mu_ and 
% _cov_ that contains the parameters of three Gaussian distributions used 
% to fit the distribution. 

function [ outimg, binimg, mu, cov ] = histfit( inimg, shift, mu, cov )
persistent zeroimg

    if (nargin < 2),
        shift = 0;
    end
    if (nargin < 3),
        mu = [0 3 5];
        cov = [1 1 1].^2;
    end

    %%
    % Start by initializing the bins, iteration variables, and dimensions.
    % Next, use the bins to produce a histogram (tmp).
    % This histogram is then smoothed by the kernel smo.
    
    N_IT = 10;
    N_COMP = length(mu);
    BINS   = -0.10 : 0.02 : 1.00;
    BIN_SZ = BINS(2)-BINS(1);

    v = inimg(:);
    tmp = fhist(v,BINS);
    smo = [1 2 1];
    smo = smo ./ sum(smo);
    hist1 = conv2(tmp,smo,'same');

    %%
    % The left hand side of the expression below is just the maximum 
    % center coordinate of the non-empty bins. If this coordinate is 
    % less than or equal to 0.25, the fly most likely cannot be detected. 
    
    if max( BINS( hist1 > 0 ) ) <= 0.25,  
        zeroimg = inimg * 0;
        outimg = zeroimg;
        binimg = zeroimg;
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
        % Shift is a flag indicating that a more detailed histogram  
        % analysis based on the Gaussian mixture model is requested. 
        
        if (shift > -99),
            
            %%
            % Compute the histogram and smooth out the result. 
            % Remove the content of the first three bins, which may
            % be spillover background pixels. 
            
            hist1 = conv2(fhist(v,BINS),smo,'same');
            hist1(1:3) = 0;
            shist = sum(hist1);
            if shist > 0;

                % Compute the adjusted histogram
                
                histo = (1/BIN_SZ)*hist1/shist;
                old_lik = 0;
                lik = zeros(N_COMP,length(BINS));

                %%
                % The following section identifies the parameters of 
                % the three-component Gaussian Mixture Model. 
                % The loop terminates when the maximum number of iteration
                % is reached, or when the parameters converge.
                
                for it=1:N_IT,
                    
                    %%
                    % First, update the component probability densities,
                    % and and the mixture density. Next, normalize the 
                    % component densities by the mixture density. 
                    
                    for c=1:N_COMP,
                        lik(c,:) = gauss( mu(c), cov(c), BINS' );
                    end;
                    tot_lik = sum(lik,1);
                    lik_norm = lik ./ repmat(tot_lik,N_COMP,1);
                    
                    %%
                    % The parameters are then updated, ensuring that the 
                    % standard deviation exceeds the bin step size.
                    % Expectation values are computed using the product of 
                    % the component densities and the actual histogram
                    
                    for c=1:N_COMP,
                        mu(c)  = sum( lik_norm(c,:).*histo.* BINS ) / ... 
                                 sum( lik_norm(c,:).*histo );
                        cov(c) = sum( lik_norm(c,:).*histo.* (BINS-mu(c)).^2 ) / ... 
                                 sum( lik_norm(c,:).*histo );
                        cov(c) = max(cov(c),0.02^2);
                    end;
                    
                    %%
                    % Check for convergence 
                    
                    tlik = sum(tot_lik);
                    if abs(tlik - old_lik) < 0.1, break; end
                    old_lik = tlik;
                    
                end;

                %%
                % Compute the indices of the maximum of the three 
                % unnormalized component Gaussian densities 
                
                indm = zeros(N_COMP,1);
                for c=1:N_COMP,
                    tmp = lik(c,:) * sum(lik_norm(c,:).*histo) * BIN_SZ;
                    if (tmp > -Inf), 
                        in = find( tmp == max(tmp) );
                        indm(c) = in(1);
                    end
                end

                %%
                % Define a range spanning the minimum and maximum of the 
                % component density peak location, shrunk by one bin on
                % each side. Calculate the absolute value of the difference 
                % between the second and the last component density vector 
                % weighted by lik_norm.*histo inside this range. 
                % Find ind_cross, the index (relative to the beginning of 
                % this range) where the difference is minimized. 
                % When shifted relative to the first bin (i.e., by
                % min(indm) ), ind_cross is the index at which the 2nd and 
                % last component density intersect. 
                
                diff = [];
                if indm,
                    diff = abs( ...
                        lik( 2, min(indm)+1 : max(indm)-1 ) * ... 
                        sum( lik_norm( 2, min(indm)+1 : max(indm)-1 ) .* ...
                             histo( min(indm)+1 : max(indm)-1 ) ) ...
                      - lik( N_COMP, min(indm)+1 : max(indm)-1 ) * ... 
                        sum( lik_norm( N_COMP, min(indm)+1 : max(indm)-1 ) .* ... 
                             histo( min(indm)+1 : max(indm)-1) ) ); 
                end
                if diff, 
                    ind_cross = find( diff == min(diff) ) + min(indm); shift = shift + 10;
                else
                    ind_cross = 1;
                end

            else
                ind_cross = 1; shift = 0; 
            end
        else
            ind_cross = 1; shift = 0; 
        end
        
        %%
        % The bin value at (ind_cross - shift) is the thresholding value 
        % used to mask the input image.
        
        binimg = inimg > BINS(ind_cross(1)-min(ind_cross(1)-1,shift));
        outimg = inimg .* binimg;

    end

end