%% Analyze_wing_mov
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
%% Analyze wing movements

function [img,phivelo,rotdir,dis_arr] = analyze_wing_mov(obj,phivelo,img,igen,add,tstep,min_vortho,min_vparal,min_vphi,xv,yv,chamber,oneobj)

rotdir = []; dis_arr = [];
maxit = 4;
if tstep>0, 
    nt = add/tstep*2;
    if ~numel(img{igen}),
        for j=1:nt, img{igen,j} = zeros(numel(xv),numel(yv)); end
    end
else
    nt = 1;
    if ~numel(img), img = zeros(numel(xv),numel(yv)); end
end
cnt1 = 0;
for j=1:length(obj{igen}.tim),
    if numel(obj{igen}.tim{j}),
        for l=1:length(obj{igen}.tim{j}),
            lobj = length(obj{igen}.do1{j}{l});
            % Consider in case a wing extension lasting longer than 
            % 10 frames
            if lobj > 10,
                for l1=1:nt,
                    if (l1-1)*tstep+add <= lobj,
                        cnt1 = cnt1 + 1;
                        % Compute either time series of histograms/plots or
                        % integrate over the full assay period (tstep = 0)
                        if tstep>0,
                            inddat = (l1-1)*tstep+1:(l1-1)*tstep+add;
                        else
                            inddat = 1:lobj;
                        end
                        % Extract fly positions (x,y), compute azimuthal (v_az) and 
                        % parallel (v_pa) fly velocities, extract relative fly orientation 
                        % (phi), compute angular velocities (phi_velo)
                        x = [obj{igen}.x1{j}{l}(inddat) ; obj{igen}.x2{j}{l}(inddat)]; 
                        y = [obj{igen}.y1{j}{l}(inddat) ; obj{igen}.y2{j}{l}(inddat)];
                        [v_az,v_pa,phi] = v_az_pa(x,y,chamber,oneobj);                        
                        if oneobj,
                            r = sqrt((x(1,:)'-chamber.width/2).^2 + (y(1,:)'-chamber.height/2).^2);
                        else
                            r = obj{igen}.dis{j}{l}(inddat)';
                        end
                        data = obj{igen}.do1{j}{l}(inddat(2:end-1));
                        [data,phi_velo] = smooth_angle(data,maxit);
                        phi = phi(maxit+2:end-maxit-1);
                        % Compute fly distances (dis); normalize fly distances to fly
                        % distance at center frame (dis_arr)
                        if (nargout > 3),
                            r = r(maxit+2:end-maxit-1);
                            dis = conv2(r,[.25 .5 .25],'same')';
                            dis_arr = [dis_arr ; dis / dis(round(length(dis)/2))];
                        end
                        v_az = v_az(1,maxit+2:end-maxit-1);
                        v_pa = v_pa(1,maxit+2:end-maxit-1);
                        % Range/azimuth movement of one fly relative to the other
                        % one, without considering the body-orientation.
                        % Consider only events with minimum azimuthal and parallel 
                        % velocity and with a angular velocity between 10-30 degrees
                        in = find((abs(v_az) > min_vortho) & (abs(v_pa) > min_vparal));
                        if numel(in),
                            data = data(in); phi_velo = phi_velo(in);
                            phi = phi(in); r = r(in);
                            in = find(abs(phi_velo) < 30);
                            if numel(in)>10,
                                data = data(in); phi_velo = phi_velo(in);
                                phi = phi(in); r = r(in);
                                % Prepare 2d histograms of positions of a fly relative to
                                % the opponent who is extending a wing. The wing-extending
                                % fly is located in the center of the histogram.
                                dori = phi - data + 90;
                                in = find(dori > 180); if numel(in), dori(in) = dori(in) - 360; end
                                in = find(dori < -180); if numel(in), dori(in) = dori(in)+ 360; end
                                dori = dori' * pi/180;
                                xy = [r .* cos(dori) r .* sin(dori)];
                                if tstep,
                                    img{igen,l1} = img{igen,l1} + hist3(xy,{xv yv});
                                else
                                    img = img + hist3(xy,{xv yv});
                                end
                                % Compute number of positive and negative angular
                                % velocities to analyze clockwise and counterclockwise
                                % circling
                                phi_velo = phi_velo(3:end-2);
                                phi_velo = phi_velo(abs(phi_velo) >= min_vphi);
                                phivelo{igen,1}{cnt1} = phi_velo;
                                np = sum(phi_velo(phi_velo>0)); nn = sum(phi_velo(phi_velo<0));
                                phivelo{igen,2}(cnt1,:) = [np nn]; %mean(phi_velo);
                            end
                        end
                    end
                end
            end
        end
    end
end
if cnt1,
    rotdir = sum(phivelo{igen,2},1);
end

end