%% Analyze_proximities
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
%% Analyze fly positions in relation to the other fly (coordinate 
%% system is shifted onto the other fly with head facing North)

function FID = analyze_proximities(obj1,obj2,dists,lunges,tussls,wings,chases,courts,dts,gen_ident,genotypes,sorted_names,params,chamber,FID) %#ok<*INUSL>

dr = .1; % [mm]
h = 5; %chamber.height/2;
w = h; %chamber.width/2;
dax = h/2;
% NUMBER OF PLOTS PER PAGE
nplt = 4;

if params.plots.heatmaps.wingrel && ~params.oneobj,
    for ilen=1:min(5,numel(params.plots.heatmaps.wingrel_vec)),
        if params.plots.heatmaps.wingrel_vec(ilen) >= 4,
            ngens = numel(gen_ident);
            % GLOBAL OR TIME SERIES?
            if (params.plots.heatmaps.wingrel_vec(ilen) == 4) || ...
                (params.plots.heatmaps.wingrel_vec(ilen) == 6),
                % Compute all positions over full movie length
                params.trans_t = params.max_frames;
                nloops = 1;
            elseif (params.plots.heatmaps.wingrel_vec(ilen) == 5) || ...
                    (params.plots.heatmaps.wingrel_vec(ilen) == 7),
                % Transition time series;  split movie into chunks of 'trans_t' seconds
                params.trans_t = 300; % [s]
                nloops = ceil(params.max_frames / params.trans_t);
                k1 = floor(sqrt(nplt)); k2 = ceil(nplt/k1);
            end
            
            xv = -w:dr:w; yv = -h:dr:h;
            xaxis = -w:dax:w;
            xtick_lab = cell(1,length(xaxis)); for i=1:length(xaxis), xtick_lab{i} = num2str(-w+(i-1)*dax); end
            icnt = 1;
            
            for igen=1:ngens,
                dt = median(dts{igen});
                if params.tim, fac = dt; else fac = 1; end
                for iloop=1:nloops,
                    if (nloops > 1 && ~mod(icnt-1,nplt)) || (nloops == 1 && igen == 1),
                        if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
                        icnt = 1;
                        figure(FID); clf;
                        set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
                    end
                    if nloops == 1,
                        h1 = subplot('Position',[mod(igen-1,params.k(2))/(params.k(2)*1.1)+.05, ...
                            0.9-fix((igen-1)/params.k(2))/(params.k(1)*1.1)-1/(params.k(1)*1.5), ...
                            1/(params.k(2)*1.5), 1/(params.k(1)*1.5)]); hold on;
                    else
                        h1 = subplot('Position',[mod(icnt-1,k2)/(k2*1.1)+.05, ...
                            0.9-fix((icnt-1)/k2)/(k1*1.1)-1/(k1*1.5), ...
                            1/(k2*1.5), 1/(k1*1.5)]); hold on;
                    end
                    
                    if (params.plots.heatmaps.wingrel_vec(ilen) == 6) || ...
                        (params.plots.heatmaps.wingrel_vec(ilen) == 7),
                        dy = wings.threat{igen}.obj1.y1-wings.threat{igen}.obj1.y2;
                        dx = wings.threat{igen}.obj1.x1-wings.threat{igen}.obj1.x2;
                        theta1 = atan2(dy,dx)'; dc1 = sqrt(dx.^2 + dy.^2);
                        dy = wings.threat{igen}.obj2.y1-wings.threat{igen}.obj2.y2;
                        dx = wings.threat{igen}.obj2.x1-wings.threat{igen}.obj2.x2;
                        theta2 = atan2(dy,dx)'; dc2 = sqrt(dx.^2 + dy.^2);                        
                        l1 = [0 wings.threat{igen}.obj1.len];
                        l2 = [0 wings.threat{igen}.obj2.len];
                    else
                        dy = obj1{igen}.y - obj2{igen}.y; dx = obj1{igen}.x - obj2{igen}.x;
                        theta = atan2(dy,dx)';
                    end
                                        
                    ind = [0 find(obj1{igen}.t(2:end)-obj1{igen}.t(1:end-1)<0) length(obj1{igen}.t)];
                    nmovs = numel(ind)-1;
                    img = zeros(numel(xv),numel(yv),nmovs);

                    cwthr1 = [0 cumsum(wings.threat{igen}.obj1.number)];
                    cwthr2 = [0 cumsum(wings.threat{igen}.obj2.number)];
                    
                    
                    for imov=1:nmovs,
                        indt = ind(imov)+1:ind(imov+1);
                        
                        if (params.plots.heatmaps.wingrel_vec(ilen) == 6) || ...
                                (params.plots.heatmaps.wingrel_vec(ilen) == 7),
                           % Wing threats only                           
                           st = sum(l1(1:cwthr1(imov)+1))+1; en = sum(l1(1:cwthr1(imov+1)+1));
                           indt1 = st:en;
                           st = sum(l2(1:cwthr2(imov)+1))+1; en = sum(l2(1:cwthr2(imov+1)+1));
                           indt2 = st:en;                           
                           if nloops > 1,
                               indtp = (wings.threat{igen}.obj1.tim(indt1) > (iloop-1)*params.trans_t) & ...
                                   (wings.threat{igen}.obj1.tim(indt1) <= iloop*params.trans_t);
                               indt1 = indt1(indtp);
                               indtp = (wings.threat{igen}.obj2.tim(indt2) > (iloop-1)*params.trans_t) & ...
                                   (wings.threat{igen}.obj2.tim(indt2) <= iloop*params.trans_t);
                               indt2 = indt2(indtp);
                           end
                           if numel(indt1),
                               dori = theta1(indt1) - wings.threat{igen}.obj1.do2(indt1)'*pi/180 + pi/2;
                           else
                               dori = [];
                           end
                           if numel(indt2),
                               dori = [dori ; theta2(indt2) - wings.threat{igen}.obj2.do2(indt2)'*pi/180 + pi/2]; %#ok<AGROW>
                           end
                           in = find(dori > pi); if numel(in), dori(in) = dori(in) - 2*pi; end
                           in = find(dori < -pi); if numel(in), dori(in) = dori(in)+ 2*pi; end
                           r = [dc1(indt1)' ; dc2(indt2)'];
                        else
                            % All fly positions
                            if nloops > 1,
                                indtp = (obj1{igen}.t(indt) > (iloop-1)*params.trans_t) & ...
                                    (obj1{igen}.t(indt) <= iloop*params.trans_t);
                                indt = indt(indtp);
                            end
                            dori = theta(indt) - obj2{igen}.do(indt)'*pi/180 + pi/2;
                            in = find(dori > pi); if numel(in), dori(in) = dori(in) - 2*pi; end
                            in = find(dori < -pi); if numel(in), dori(in) = dori(in)+ 2*pi; end
                            r = dists{igen}.distc(indt)';                            
                        end
                        
                        in = find(r);
                        r1 = r(in); dori1 = dori(in);
                        data = [r1 .* cos(dori1) , r1 .* sin(dori1)];
                        img(:,:,imov) = hist3(data,{xv yv}) * fac;
                    end
                    img = mean(img,3)';
                    
                    % Time spent or occurrence?
                    if params.tim,
                        max_val = params.range.heatmaps.tim.wing ;
                    else
                        max_val = params.range.heatmaps.occ.wing;
                    end
                    if (params.plots.heatmaps.wingrel_vec(ilen) == 6) || ...
                            (params.plots.heatmaps.wingrel_vec(ilen) == 7),
                        max_val = max_val / 2;
                    end
                    
                    %                 dhh = dists{igen}.disth(indt);
                    %                 dtt = dists{igen}.distt(indt);
                    %                 dh1t2 = dists{igen}.disth1t2(indt);
                    %                 dh2t1 = dists{igen}.disth2t1(indt);
                    % f1 faces towards side of f2
                    %     in = find(dhh < 2*dtt & ...
                    %               dhh < 2*dh1t2 & ...
                    %               dh2t1 < 2*dtt & ...
                    %               dh2t1 < 2*dh1t2 & ...
                    %               abs(dhh-dh2t1)<.5);
                    % f1 parallel to f2
                    %     in = find(abs(dhh-dtt)<.5 & ...
                    %               dhh < 2);
                    % f1 faces f2
                    %     in = find(abs(dh1t2-dh2t1)<.1 &...
                    %               dhh < 2*dtt & ...
                    %               dhh < 2);
                    %     r1 = r(in); dori1 = dori(in);
                    %     data = [r1 .* cos(dori1) , r1 .* sin(dori1)];
                    %     img = hist3(data,{-w:dr:w -h:dr:h});
                    
                    % f1 parallel to f2
                    %     in = find(abs(dh2t1-dh1t2)<.5 & ...
                    %               abs(dhh-dtt)<.5);
                    % f1 back f2
                    %     in = find(abs(dh1t2-dh2t1)<.1 &...
                    %               dtt < 2*dhh & ...
                    %               dtt < 2);
                                        
                    xlim([-w w]); ylim([-h h]); 
                    title(cell2mat(gen_ident(igen)),'FontSize',params.axisfontsize); 
                    axis on;
                    imagesc(xv, yv, img, [0 max_val]); drawnow; axis square;
                    
                    plot([0 0], [0 w/8],'w');
                    plot([-w/8 w/8], [0 0],'w');
                    set(h1,'FontSize',params.axisfontsize); set(gca,'FontSize',params.axisfontsize);
                    set(gca,'XTick',xaxis); set(gca,'XTickLabel',xtick_lab);
                    set(gca,'YTick',xaxis); set(gca,'YTickLabel',xtick_lab);
                    if (nloops == 1 && ~mod(igen-1,params.k(2))) || ...
                            ((nloops > 1 && ~mod(icnt-1,k2)) || (icnt == 1 && igen == ngens)),
                        ylabel('distance [mm]','FontSize',params.axisfontsize);
                    end
                    if (nloops == 1 && (igen > (params.k(1)-1)*params.k(2))) || ...
                            (nloops > 1 && (icnt > (k1-1)*k2)),
                        xlabel('distance [mm]','FontSize',params.axisfontsize);
                    end
                    colormap('default'); colormap([0 0 0 ; jet]);
                    
                    tit = gen_ident{igen}; in = strfind(tit,'_'); tit(in) = '-';
                    if nloops == 1,
                        title(tit);
                    else
                        title([tit ', ' num2str((iloop-1)*params.trans_t) '-' ...
                            num2str(iloop*params.trans_t) ' s']);
                    end
                    
                    if params.pdf,
                        if ((nloops == 1) && (igen == ngens)) || ...
                                ((nloops > 1) && ~mod(icnt,nplt) || (igen == length(gen_ident) && iloop == nloops)),
                            pos = get(gca,'Position');
                            ch = colorbar('location','EastOutside');
                            if params.tim,
                                set(get(ch,'YLabel'),'String','cumulative time [s]');
                            else
                                set(get(ch,'YLabel'),'String','occurence [-]');
                            end
                            set(ch,'Position',[0.93 pos(2) 0.01 pos(4)]);
                            
                            if (params.plots.heatmaps.wingrel_vec(ilen) == 6) || ...
                                    (params.plots.heatmaps.wingrel_vec(ilen) == 7),
                                titl = 'Fly head position during Wing Threats';
                            else
                                titl = 'Fly HEAD position';
                            end
                            if nloops > 1,
                                titl = [titl ' per ' num2str(params.trans_t/60) ' minute bin']; %#ok<AGROW>
                            end
                            set(FID,'Name',titl); mtit(titl,'FontSize',16','yoff',.03);
                            print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN);
                        end
                    end
                    icnt = icnt + 1;
                end
            end
        end
    end    
end