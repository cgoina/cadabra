%% OTSU 
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
% Gray-level image segmentation using Otsu's method.
%
% * Originally written by Damien Garcia 2007/08
% * Optimized for QTracker by Edwin Soedarmadji 2008
% * Changed by Heiko Dankert 11/07/2008
% 
%% 
%
% Iseg = OTSU(I,n) computes a segmented image (Iseg) containing n classes
% by means of Otsu's n-thresholding method (Otsu N, A Threshold Selection
% Method from Gray-Level Histograms, IEEE Trans. Syst. Man Cybern.
% 9:62-66;1979). Thresholds are computed to maximize a separability
% criterion of the resultant classes in gray levels.
%
% OTSU(I) is equivalent to OTSU(I,2). By default, n=2 and the
% corresponding Iseg is therefore a binary image. The pixel values for
% Iseg are [0 1] if n=2, [0 0.5 1] if n=3, [0 0.333 0.666 1] if n=4, ...
%
% [Iseg,sep] = OTSU(I,n) returns the value (sep) of the separability
% criterion within the range [0 1]. Zero is obtained only with images
% having less than n gray level, whereas one (optimal value) is obtained
% only with n-valued images.
%
% Notes:
% -----
% It should be noticed that the thresholds generally become less credible
% as the number of classes (n) to be separated increases (see Otsu's
% paper for more details).
%
% If n=2 or 3, the separability criterion is directly maximized by simply
% using the MAX function. For n values >= 4, a minimization method is
% used by means of the FMINSEARCH function.
%
% The OTSU function works with I of any size. An RGB image I will thus be
% considered as a gray-level 3D-array!
%   
% Example:
% -------
% load clown
% subplot(221)
% imshow(X/max(X(:)))
% title('Original','FontWeight','bold')
% for n = 2:4
%     Iseg = otsu(X,n);
%     subplot(2,2,n), colormap(gray)
%     imshow(Iseg)
%     title(['n = ' int2str(n)],'FontWeight','bold')
% end

function [Iseg,sep] = fotsu(I,n)

%% 
% Checking n (the number of classes). 

if nargin==1
    n = 2;
elseif n==1;
    Iseg = NaN(size(I),'single');
    sep = 0;
    return
elseif n~=abs(round(n)) || n==0
    error('n must be a strictly positive integer!')
end

%%% change made by H. Dankert, 11072008
% nbins = min(numel(I),256);
% K0 = 0.1;
% K = K0*(diag(-2*ones(n-1,1))+diag(ones(n-2,1),1)+diag(ones(n-2,1),-1));
% x0 = zeros(n-1,1); x0(n-1) = -K0;
%%%

%% 
% Calculating the probability distribution.

I = single(I); % the HIST function only accepts single or double
%%% change made by H. Dankert, 11072008
unI = sort(unique(I)); nbins = min(length(unI),256);
%%%
if nbins==n
    Iseg = single(ones(size(I)));
    graycol = linspace(0,1,n);
    for i = 1:n-1
        Iseg(I==unI(i)) = graycol(i);
    end
    sep = 1;
    return
elseif nbins<n
    Iseg = NaN(size(I),'single');
    sep = 0;
    return
%%% change made by H. Dankert, 11072008
elseif nbins<256
    [histo,pixval] = fhist(I(:),unI);
else
    [histo,pixval] = fhist(I(:),256);
% else
%     [histo,pixval] = fhist(I(:),nbins);
%%%
end
P = histo/sum(histo);

%% 
% Compute the zeroth- and first-order cumulative moments

w = cumsum(P);
mu = cumsum((1:nbins).*P);

%% 
% Compute the maximal sigmaB^2 and segmented image

if n==2
    sigma2B =...
        (mu(end)*w(2:end-1)-mu(2:end-1)).^2./w(2:end-1)./(1-w(2:end-1));
    
    [maxsig,k] = max(sigma2B);
        
    % segmented image
    Iseg = single(ones(size(I)));
    Iseg(I<=pixval(k+1)) = 0;
    
    % separability criterion
    sep = maxsig/sum(((1:nbins)-mu(end)).^2.*P);
    
elseif n==3
    w0 = w;
    w2 = fliplr(cumsum(fliplr(P)));
    [w0,w2] = ndgrid(w0,w2);
   
    mu0 = mu./w;
    mu2 = fliplr(cumsum(fliplr((1:nbins).*P))./cumsum(fliplr(P)));
    [mu0,mu2] = ndgrid(mu0,mu2);
    
    w1 = 1./(1-w0-w2);
    w1(w1<=0) = 0;

    tmp0 = mu0-mu(end); tmp2 = mu2-mu(end); w00 = w0.*tmp0; w22 = w2.*tmp2;
    sigma2B = w00.*tmp0 + w22.*tmp2 + (w00 + w22).^2.*w1;
    sigma2B(isnan(sigma2B)) = 0; % zeroing if k1>= k2;
    
    [maxsig,k] = max(sigma2B(:));
    [k1,k2] = ind2sub([nbins nbins],k);
    
    % segmented image
    Iseg = single(ones(size(I)));
    Iseg(I<=pixval(k1)) = 0;
    Iseg(I>pixval(k1) & I<=pixval(k2)) = 0.5;
    
    % separability criterion
    sep = maxsig/sum(((1:nbins)-mu(end)).^2.*P);

%% 
% Threshold positions are adjusted using a horizontal mass-spring
% system. Threshold positions can thus be adjusted by modulating the
% force exerted on each mass.

else
    %%% change made by H. Dankert, 11072008
    K0 = 0.1;
    K = K0*(diag(-2*ones(n-1,1))+diag(ones(n-2,1),1)+diag(ones(n-2,1),-1));
    x0 = zeros(n-1,1); x0(n-1) = -K0;
    %%%
    
    [F,y,exitflag] = otsuminsearch( ones(n-1,1), n, nbins, K, P, x0 ); %#ok<NASGU>
    
    k = K\(x0+F-1);
    k = round(k*(nbins-1)+1);
    k = [0;k;nbins];
    
    % segmented image
    Iseg = single(ones(size(I)));
    graycol = linspace(0,1,n);
    for i = 1:n-1
        Iseg(I>=pixval(k(i)+1) & I<=pixval(k(i+1))) = graycol(i);
    end

    % separability criterion
    sep = 1-y; 
    
end
end

%% 
% This is the actual function to be minimized if n>=4

function y = sig_func(F, nclass, nbins, K, P, x0)

k = K\(x0+F-1);
k = round(k*(nbins-1)+1);

if any(k<1 | k>nbins) || any(diff(k)<1)
    y = 1;
    return
end

muT = sum((1:nbins).*P);
sigma2T = sum(((1:nbins)-muT).^2.*P);

k = [0;k;nbins];
sigma2B = 0;
for i = 1:nclass
    kv = (k(i)+1:k(i+1)); 
    Pkv = P(kv);
    wi = sum(Pkv);
    if wi==0
        y = 1;
        return
    end
    mui = sum(kv.*Pkv)/wi;
    sigma2B = sigma2B + wi*(mui-muT)^2;
end

y = 1-sigma2B/sigma2T; % within the range [0 1]
end

%% FFMINSEARCH 
% Optimized Multidimensional unconstrained nonlinear minimization (Nelder-Mead).
%
% * Optimization by Edwin Soedarmadji 2008
% 
% Reference: Jeffrey C. Lagarias, James A. Reeds, Margaret H. Wright,
% Paul E. Wright, "Convergence Properties of the Nelder-Mead Simplex
% Method in Low Dimensions", SIAM Journal of Optimization, 9(1):
% p.112-147, 1998.
% Copyright 1984-2007 The MathWorks, Inc.

function [x,fval,exitflag] = otsuminsearch(x, nclass, nbins, K, P, x0 )

% Check for non-double inputs
if ~isa(x,'double')
  error('MATLAB:fminsearch:NonDoubleInput', ...
         'FMINSEARCH only accepts inputs of data type double.')
end

% Initialize parameters
n = numel(x);
tolx = 1e-4;
tolf = 1e-4;
maxfun = 25;
maxiter = 20;
rho = 1; chi = 2; psi = 0.5; sigma = 0.5;
onesn = ones(1,n);
two2np1 = 2:n+1;
one2n = 1:n;

% Set up a simplex near the initial guess.
xin = x(:); % Force xin to be a column vector
v = zeros(n,n+1); fv = zeros(1,n+1);
v(:,1) = xin;    % Place input guess in the simplex! (credit L.Pfeffer at Stanford)
x(:) = xin;    % Change x to the form expected by sig_func
fv(:,1) = sig_func(x, nclass, nbins, K, P, x0 );
itercount = 0;

% Following improvement suggested by L.Pfeffer at Stanford
usual_delta = 0.05;             % 5 percent deltas for non-zero terms
zero_term_delta = 0.00025;      % Even smaller delta for zero elements of x
for j = 1:n
    y = xin;
    if y(j) ~= 0
        y(j) = (1 + usual_delta)*y(j);
    else
        y(j) = zero_term_delta;
    end
    v(:,j+1) = y;
    x(:) = y; 
    f = sig_func(x, nclass, nbins, K, P, x0 );
    fv(1,j+1) = f;
end

% sort so v(1,:) has the lowest function value
[fv,j] = sort(fv);
v = v(:,j);

%how = 'initial simplex';
itercount = itercount + 1;
func_evals = n+1;

% Main algorithm: iterate until 
% (a) the maximum coordinate difference between the current best point and the 
% other points in the simplex is less than or equal to TolX. Specifically,
% until max(||v2-v1||,||v2-v1||,...,||v(n+1)-v1||) <= TolX,
% where ||.|| is the infinity-norm, and v1 holds the 
% vertex with the current lowest value; AND
% (b) the corresponding difference in function values is less than or equal
% to TolFun. (Cannot use OR instead of AND.)
% The iteration stops if the maximum number of iterations or function evaluations 
% are exceeded
while func_evals < maxfun && itercount < maxiter
    if max(abs(fv(1)-fv(two2np1))) <= max(tolf,10*eps(fv(1))) && ...
            max(max(abs(v(:,two2np1)-v(:,onesn)))) <= max(tolx,10*eps(max(v(:,1))))
        break
    end
    
    % Compute the reflection point
    
    % xbar = average of the n (NOT n+1) best points
    xbar = sum(v(:,one2n), 2)/n;
    xr = (1 + rho)*xbar - rho*v(:,end);
    x(:) = xr; 
    fxr = sig_func(x, nclass, nbins, K, P, x0 );
    func_evals = func_evals+1;
    
    if fxr < fv(:,1)
        % Calculate the expansion point
        xe = (1 + rho*chi)*xbar - rho*chi*v(:,end);
        x(:) = xe; 
        fxe = sig_func(x, nclass, nbins, K, P, x0 );
        func_evals = func_evals+1;
        if fxe < fxr
            v(:,end) = xe;
            fv(:,end) = fxe;
            %how = 'expand';
        else
            v(:,end) = xr;
            fv(:,end) = fxr;
            %how = 'reflect';
        end
    else % fv(:,1) <= fxr
        if fxr < fv(:,n)
            v(:,end) = xr;
            fv(:,end) = fxr;
            %how = 'reflect';
        else % fxr >= fv(:,n)
            % Perform contraction
            if fxr < fv(:,end)
                % Perform an outside contraction
                xc = (1 + psi*rho)*xbar - psi*rho*v(:,end);
                x(:) = xc; 
                fxc = sig_func(x, nclass, nbins, K, P, x0 );
                func_evals = func_evals+1;
                
                if fxc <= fxr
                    v(:,end) = xc;
                    fv(:,end) = fxc;
                    how = 'contract outside';
                else
                    % perform a shrink
                    how = 'shrink';
                end
            else
                % Perform an inside contraction
                xcc = (1-psi)*xbar + psi*v(:,end);
                x(:) = xcc; 
                fxcc = sig_func(x, nclass, nbins, K, P, x0 );
                func_evals = func_evals+1;
                
                if fxcc < fv(:,end)
                    v(:,end) = xcc;
                    fv(:,end) = fxcc;
                    how = 'contract inside';
                else
                    % perform a shrink
                    how = 'shrink';
                end
            end
            if strcmp(how,'shrink')
                for j=two2np1
                    v(:,j)=v(:,1)+sigma*(v(:,j) - v(:,1));
                    x(:) = v(:,j); 
                    fv(:,j) = sig_func(x, nclass, nbins, K, P, x0 );
                end
                func_evals = func_evals + n;
            end
        end
    end
    [fv,j] = sort(fv);
    v = v(:,j);
    itercount = itercount + 1;
end   % while

x(:) = v(:,1);
fval = fv(:,1);

if func_evals >= maxfun
    exitflag = 0;
elseif itercount >= maxiter
    exitflag = 0;
else
    exitflag = 1;
end
end