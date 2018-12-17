%% Fly_bodymeas
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
% * Written by Heiko Dankert 2006/2007
% * Documentation by Edwin Soedarmadji 10/2008
%
%% B.7 Fly Body Measurements
% At this point, the program has determined the number of flies in the
% chamber, updated the appropriate flags ( _ind_twoobj_, etc.) and
% segmented the fly body with pixel indices (f011). 
% In this section, an ellipse is fit to the segemented pixels of the 
% fly body to do measurements such as lengths, areas, positions, etc.

function [obj1,obj_1] = fly_bodymeas(obj1,obj_1,img,f011,fly_ind1,chamber,params,scale,arr_lim,ind_head)

if numel(f011) * scale.xy > arr_lim,
    tmp = img * 0; tmp(f011) = img(f011);
    tmp(fly_ind1.r,fly_ind1.c) = bwmorph( bwmorph( ...
        histfit(tmp(fly_ind1.r,fly_ind1.c))>0, 'dilate',2), 'close');
    tmp1 = f011;
    f011 = find(tmp);
    if numel(f011)<5,
        f011 = tmp1;
    end
end

xe = chamber.x_arr(f011);
ye = chamber.y_arr(f011);
obj1.length = sqrt((max(xe)-min(xe)).^2+(max(ye)-min(ye)).^2);
obj1.Area = numel(f011) * scale.xy;

xe = chamber.x_arr(f011);
ye = chamber.y_arr(f011);
obj1.exc = obj_1.exc;
obj1.eyc = obj_1.eyc;
obj1.A = obj_1.A;
obj1.B = obj_1.B;
obj1.phi = obj_1.phi;
if ( sum(abs(xe(2:end)-xe(1:end-1))) > 0.01 ) && ...
        ( sum(abs(ye(2:end)-ye(1:end-1))) > 0.01 ),
    [Xce1,Yce1,A01,B01,Phi01,P] = ellipsefit(xe,ye);
    if numel(P),
        if (abs(Xce1) <= chamber.s_arr(2)*scale.x) && ...
                (abs(Yce1) <= chamber.s_arr(1)*scale.y) && ...
                (abs(A01) < 3) && (abs(B01) < 3),
            obj1.exc = real(Xce1);
            obj1.eyc = real(Yce1);
            obj1.A = real(A01);
            obj1.B = real(B01);
            if obj1.A/obj_1.A > 2.5 && ~isinf(obj1.A/obj_1.A), 
                obj1.A = obj_1.A; obj1.B = obj_1.B; 
            end
            obj1.phi = real(Phi01);
        end
    end
end
obj1.xc = obj1.exc;
obj1.yc = obj1.eyc;
obj1.area = pi * obj1.A * obj1.B;
obj1.ra = sqrt( ...
    (obj1.xc - chamber.s_arr(2)/2 * scale.x)^2 + ...
    (obj1.yc - chamber.s_arr(1)/2 * scale.y)^2 );
if isempty(obj_1),
    obj_1.xc = obj1.xc;
    obj_1.yc = obj1.yc;
    obj_1.head = -99;
end
end

%% C.1 EllipseFit 
% Stable Direct Least Squares Ellipse Fit to Data.
% [Xc,Yc,A,B,Phi,P]=ELLIPSEFIT(X,Y) finds the least squares ellipse that
% best fits the data in X and Y. X and Y must have at least 5 data points.
% _Xc_ and _Yc_ are the x- and y-axis center of the ellipse respectively.
% _A_ and _B_ are the major and minor axis of the ellipse respectively.
% _Phi_ is the radian angle of the major axis with respect to the x-axis.
% _P_ is a vector containing the general conic parameters of the ellipse.
function [varargout] = ellipsefit( x , y )
global params;

x=x(:);                                 % convert data to column vectors
y=y(:);
if numel(x)~=numel(y) || numel(x)<5
   error('X and Y Must be the Same Length and Contain at Least 5 Values.')
end

D1=[x.*x x.*y y.*y];                    % quadratic terms
D2=[x y ones(size(x))];                 % linear terms
S1=D1'*D1;
S2=D1'*D2;

[Q2,R2]=qr(D2,0);
T=-R2\(R2'\S2');                        % -inv(S3) * S2'

M=S1+S2*T;
CinvM=[M(3,:)/2; -M(2,:); M(1,:)/2];
[V,na]=eig(CinvM); %#ok<NASGU>
c=4*V(1,:).*V(3,:) - V(2,:).^2;
A1=V(:,c>0);
P=[A1; T*A1];

                                        % correct signs if needed
if numel(P),
    P=sign(P(1))*P;

    Phi=atan(P(2)/(P(3)-P(1)))/2;
    c=cos(Phi);
    s=sin(Phi);

                                        % rotate the ellipse parallel to x-axis
                                        
    Pr=zeros(6,1);
    Pr(1)=P(1)*c*c - P(2)*c*s + P(3)*s*s;
    Pr(2)=2*(P(1)-P(3))*c*s + (c^2-s^2)*P(2);
    Pr(3)=P(1)*s*s + P(2)*s*c + P(3)*c*c;
    Pr(4)=P(4)*c - P(5)*s;
    Pr(5)=P(4)*s + P(5)*c;
    Pr(6)=P(6);

                                        % extract other data
                                        
    XcYc=[c s;-s c]*[-Pr(4)/(2*Pr(1));-Pr(5)/(2*Pr(3))];
    Xc=XcYc(1);
    Yc=XcYc(2);
    F=-Pr(6) + Pr(4)^2/(4*Pr(1)) + Pr(5)^2/(4*Pr(3));
    AB=sqrt(F./Pr(1:2:3));
    A=AB(1);
    B=AB(2);
    Phi=-Phi;
    if A<B                              % x-axis not major axis, so rotate it pi/2
        Phi=Phi-sign(Phi)*params.hpi;
        A=AB(2);
        B=AB(1);
    end
    S.Xc=Xc;
    S.Yc=Yc;
    S.A=A;
    S.B=B;
    S.Phi=Phi;
    S.P=P;
else
    S.Xc=[]; S.Yc=S.Xc; S.A=S.Xc; S.B=S.Xc; S.Phi=S.Xc; S.P=S.Xc;
end
if nargout==1
    varargout{1}=S;
else
    outcell=struct2cell(S);
    varargout=outcell(1:nargout);
end
end
