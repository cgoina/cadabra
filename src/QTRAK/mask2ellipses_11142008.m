%% Mask2Ellipses
% Copyright (C) 2008 Lihi Zelnik, California Institute of Technology

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

function [Lfinal,STATSfinal] = mask2ellipses(L,Img,SolidityThresh,minRegionSize)

%  [new_mask,STATS] = mask2ellipse(mask,SolidityThresh,minRegionSize)
%
% Given a mask image represent it with best fitted ellipses.
% mask is assumed to be a labeling matrix with 0 marking background.
% SolidityThresh is the solidity threshold for discarding regions.
% minRegionSize is minimum region size. Smaller ones are discarded.
%
% Lihi Zelnik, May 2006, Caltech

if nargin < 2, Img = 0; end

if nargin < 3, SolidityThresh = 0.98; end

if nargin < 4, minRegionSize = 30; end

N = max(L(:)); % number of regions


%% Fit a single ellipse to each region
STATS  = regionprops(L, 'Solidity','PixelList','PixelIdxList');

if SolidityThresh==0
    Lfinal = L;
    STATSfinal = STATS;
end

%% go over regions and split into two those with low solidity
% also - remove too small regions (less than minRegionSize pixels)
Lnew = L;
Nnew = N;
% addpath('../netlab');
% these are the common options for all the optimizations
foptions(18) = 0; foptions(2:3) = 0.0001;  foptions(17) = 0.0001;
opts = foptions; opts(1) = 0; opts(14) = 50; % for the EM
init_opts = foptions; opts(1) = 0; opts(14) = 150; % for the kmean
for i=1:N,
    if length(STATS(i).PixelIdxList) < minRegionSize,
%         disp(['Discarding of too small ellipse ' num2str(i)]);
        Idx = STATS(i).PixelIdxList;
        Lnew(Idx) = 0; % set pixels of small region to background
        Nnew = Nnew - 1;
    else
        if STATS(i).Solidity < SolidityThresh
            %% Bad ellipses are splitted into K smaller ones using EM
            K = 2;     % new number of ellipses
            X = STATS(i).PixelList;

            %% centralize and scale the data between -1 and 1
%             X = X - repmat(mean(X),size(X,1),1);
            X = X./repmat(max(abs(X)),size(X,1),1);
            X = X/max(max(abs(X)));

            %% EM
            if numel(Img) > 1,
                tmp = Img(STATS(i).PixelIdxList);
                X = [X,tmp - min(tmp)];
                mix_std0 = gmm(3, K, 'full'); % create the EM object with two centers
            else
                mix_std0 = gmm(2, K, 'diag');
            end
            mix_std0 = gmminit(mix_std0, X, init_opts); % initial estimate from k-means
            mix_std1 = gmmem(mix_std0, X, opts); %, L, STATS(i)); % em fitting
            a = gmmactiv(mix_std1,X); %% activation level
            [tmp,assignment] = max(a'); %% assign pixels to nearest Gaussian

            % set segmentation result in Lnew
            Idx = STATS(i).PixelIdxList;
            
            for k=1:max(assignment),
                Nnew = Nnew + 1; % add a new region
                inclass = find(assignment==k);
                Lnew(Idx(inclass)) = Nnew;
            end
        end
    end
end

%% Fit again a single ellipse to each region
STATSnew  = regionprops(Lnew, 'Solidity','PixelIdxList');

%% go over new ellipses and keep only those with good solidity
Lfinal = Lnew;
% Nfinal = Nnew;
% for i=1:length(STATSnew)
%     if STATSnew(i).Solidity < SolidityThresh
% %        disp(['Discarding of ellipse ' num2str(i)]);
%         Idx = STATSnew(i).PixelIdxList;
%         Lfinal(Idx) = 0; % set pixels of bad region to background
%         Nfinal = Nfinal - 1;
%     end
% end

%% Fit yet again a single ellipse to each region
STATSfinal  = regionprops(Lfinal, 'Centroid','Orientation','MajorAxisLength','MinorAxisLength','Solidity','PixelIdxList');


