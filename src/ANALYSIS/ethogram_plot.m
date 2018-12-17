%% Ethogram_plot
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
%% Plot an ethogram

% PLOT ETHOGRAM
function ethogram_plot(h,lab)

% PARAMETERS FOR DIFFERENT ETHOGRAM STYLES
bool_val = 0; % display transition probability values
bool_heat = 0; % display transition prob. (arrows) color-coded  (=1) or by thickness (=0)
simplify = 0; % constant action circle diameter (=1) or logarithmically scaled (=0)?
circ_dia = 0.22;
half_arrows = 0; % show half-arrows
thresh = 0.1; % minimum transition probability for an arrow to get plotted
fntsize = 8;
norm_arrow = 8; % graphical normalization of arrows
circ_selfprob = 0; % plot self transition as boxes (0) or circles (1)

hold on;
% Postion of action circles
coords.tussle = [.97 .78];
coords.wingthr = [.62 1.06];
coords.lunge = [.16 .98];

coords.chase = [-.15 .54];

coords.wingext = [.06 .12];
coords.circl = [.9 .3];
coords.circw = [.51 -0.08];

% ENABLE THIS LINE TO GET 'NO BEHAVIORS' PLOTTED!
% coords.nobeh = [.5 .5];

% Field names
fieldn = fieldnames(coords);
% Replacement for field names
fieldn1= {'tussle','wing threat','lunge','chase','wing extension','circling','wing & circle','no behavior'};

% COLOR-CODING OF CIRCLES, ARROWS AND TEXT
if bool_heat,
    step = 10; 
    colheat = colormap(jet(100/step)); 
    yvec = zeros(1,100/step+1); for i=1:100/step+1, yvec(i) = (i-1)*step/100; end
    if bool_heat, colorbar('Location','East','YTick',1:100/step+1,'YTickLabel',yvec); end
    col = ones(3,length(fieldn))*.8;
    colt = ones(3,length(fieldn))*.6;
else
    COL = lines(length(fieldn))';
    col = [0 0 1; .9 0 .9; 1 .1 .1; .3 .7 1; 0 1 0; 1 .6 .1; .6 .6 .4; .2 .2 .2]';
    coltxt = [1 1 1; 1 1 1; 1 1 1; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 1 1 1]';
    colt= 0.2 + 0.6 * col;
end
norm = 30;

icnt = 0; txt.x = []; txt.y = []; txt.p = []; txt.c = [];

minmax = [-0.5 1.5];

% PLOT TRANSITION ARROWS
for i=1:length(fieldn),
    indx = find(strncmp(lab,fieldn{i},5));
    colsum = sum(sum(h(indx,:)+0)); 
    % Action circles with constant diameter or logarithmically scaled?
    if simplify, w = circ_dia; else w = log2(ceil(colsum)*10)/norm; end
    xy = coords.(fieldn{i});
    for ii=1:length(fieldn),
        indx2 = find(strncmp(lab,fieldn{ii},5));
        trans_prob = (sum(sum(h(indx,indx2)+0)))/colsum;
        if trans_prob >= thresh && ~isnan(trans_prob),
            colsum2 = sum(sum(h(indx2,:)+0)); 
            % Action circles with constant diameter or logarithmically scaled?
            if simplify, w2 = circ_dia; else w2 = log2(ceil(colsum2)*10)/norm; end
            if i ~= ii,
                xy2 = coords.(fieldn{ii});
            else
                if xy(1) >= .5 && xy(2) >=.5,
                   fac = [1 1]/sqrt(2);
                elseif xy(1) < .5 && xy(2) >= .5,
                    fac = [-1 1]/sqrt(2);
                elseif xy(1) < .5 && xy(2) < .5,
                    fac = [-1 -1]/sqrt(2);
                elseif xy(1) >= .5 && xy(2) < .5,
                    fac = [1 -1]/sqrt(2);
                end
                if circ_selfprob,
                    xy2 = xy; w = 0;
                else
                    xy2 = xy + (sqrt(2)*w + .15) .* fac;
                    w = w/sqrt(2);
                end
            end
            theta = atan2(xy2(2)-xy(2),xy2(1)-xy(1));
            s = sin(theta); c = cos(theta);
            if half_arrows, 
                shift = 1/500; %1/20;
            else
                shift = 1/20;
            end
            if i == ii, shift = 0; end
%             if i ~= ii,
                if ~bool_heat, 
                    if bool_val,
                        th = 0.1/norm_arrow; 
                    else
                        th = trans_prob/norm_arrow; %sum(sum(h(indx,indx2)+0))/norm/4;
                    end
                    co = col(:,i)';
                else
                    th = .02;
                    co = colheat(ceil(length(colheat)*trans_prob),:);
                end
                if circ_selfprob && i == ii, th = 0; end
                arrow.width = log10(th*4000)/80;
                arrow.length = log10(th*4000)/20;
                alc = arrow.length * c; als = arrow.length * s;
                xy2 = xy2 - [alc als];
                xs = -1 * s * shift; ys = c * shift;
                xth1 = -s * th; yth1 = c * th;
                xth2 = s * th; yth2 = -c * th;
                fac = -.1;

                if half_arrows,
                    arrow.x1 = xy2(1)-c*w2/2+xs;
                    arrow.y1 = xy2(2)-s*w2/2+ys;
                    arrow.x2 = xy2(1)-c*w2/2+xs+xth1+alc*fac - arrow.width * cos(pi/2-theta);
                    arrow.y2 = xy2(2)-s*w2/2+ys+yth1+als*fac + arrow.width * sin(pi/2-theta);
                    arrow.x3 = xy2(1)-c*w2/2+xs + arrow.length * cos(theta);
                    arrow.y3 = xy2(2)-s*w2/2+ys + arrow.length * sin(theta);

                    x1 = xy(1)+c*w/2+xs+xth1; x2 = xy2(1)-c*w2/2+xs+xth1;
                    y1 = xy(2)+s*w/2+ys+yth1; y2 = xy2(2)-s*w2/2+ys+yth1;
                    patch([x1 x2 xy2(1)-c*w2/2+xs ...
                        xy(1)+c*w/2+xs xy(1)+c*w/2+xs+xth1],...
                        [y1 y2 xy2(2)-s*w2/2+ys ...
                        xy(2)+s*w/2+ys xy(2)+s*w/2+ys+yth1],co,'EdgeColor','none'); hold on;
                    patch([arrow.x1 arrow.x2 arrow.x3 arrow.x1],...
                        [arrow.y1 arrow.y2 arrow.y3 arrow.y1],co,'EdgeColor','none');
                else
                    arrow.x1 = xy2(1)-c*w2/2+xs+xth2+alc*fac + arrow.width * cos(pi/2-theta);
                    arrow.y1 = xy2(2)-s*w2/2+ys+yth2+als*fac - arrow.width * sin(pi/2-theta);
                    arrow.x2 = xy2(1)-c*w2/2+xs+xth1+alc*fac - arrow.width * cos(pi/2-theta);
                    arrow.y2 = xy2(2)-s*w2/2+ys+yth1+als*fac + arrow.width * sin(pi/2-theta);
                    arrow.x3 = xy2(1)-c*w2/2+xs + arrow.length * cos(theta);
                    arrow.y3 = xy2(2)-s*w2/2+ys + arrow.length * sin(theta);

                    x1 = xy(1)+c*w/2+xs+xth1; x2 = xy2(1)-c*w2/2+xs+xth1;
                    y1 = xy(2)+s*w/2+ys+yth1; y2 = xy2(2)-s*w2/2+ys+yth1;
                    patch([x1 x2 xy2(1)-c*w2/2+xs+xth2 ...
                        xy(1)+c*w/2+xs+xth2 xy(1)+c*w/2+xs+xth1],...
                        [y1 y2 xy2(2)-s*w2/2+ys+yth2 ...
                        xy(2)+s*w/2+ys+yth2 xy(2)+s*w/2+ys+yth1],co,'EdgeColor','none'); hold on;
                    if i ~= ii,
                        patch([arrow.x1 arrow.x2 arrow.x3 arrow.x1],...
                            [arrow.y1 arrow.y2 arrow.y3 arrow.y1],co,'EdgeColor','none');
                    end
                end
%             end

            icnt = icnt + 1;
            if i ~= ii,
                fac = -2;
                xs = fac * c * shift; ys = fac * s * shift;
                txt.x(icnt) = (x1+x2+xs)/2;
                txt.y(icnt) = (y1+y2+ys)/2;
            else
                fac = 1/3; cst = 0.02;
                if xy(1) >= .5 && xy(2) >=.5,
                    pos.x = (2*w*fac+cst); pos.y = (2*w*fac+cst);
                elseif xy(1) < .5 && xy(2) >= .5,
                    pos.x = -(2*w*fac+cst); pos.y = (2*w*fac+cst);
                elseif xy(1) < .5 && xy(2) < .5,
                    pos.x = -(2*w*fac+cst); pos.y = -(2*w*fac+cst);
                elseif xy(1) >= .5 && xy(2) < .5,
                    pos.x = (2*w*fac+cst); pos.y = -(2*w*fac+cst);
                end
                txt.x(icnt) = xy(1)+pos.x;
                txt.y(icnt) = xy(2)+pos.y;
            end
            txt.c(:,icnt) = colt(:,i);
            txt.p(icnt) = trans_prob;
        end
    end
end


% PLOT ACTION CIRCLES AND TEXT
% if ~bool_heat, colormap([lines ; 1 1 1]); end
for i=1:length(fieldn),
    indx = find(strncmp(lab,fieldn{i},5));
    colsum = sum(sum(h(indx,:)+0)); 
    % Action circles with constant diameter or logarithmically scaled?
    if simplify, w = circ_dia; else w = log2(ceil(colsum)*10)/norm; end
    if ceil(colsum) >= 1,
        if w == 0,
            w = 0.001;
        end
        xy = coords.(fieldn{i});

        trans = sum(sum(h(indx,indx)+0))/sum(sum(h(indx,:)+0));
        
        if trans > thresh,
            if ~bool_heat,
                w_self = trans/norm_arrow;
                co = col(:,i)';
            else
                w_self = .1;
                co = colheat(ceil(length(colheat)*trans_prob),:);
            end

            d1 = w_self*2;
            shfac = 1/sqrt(2)/2; shfac2 = (1 - 1/sqrt(2))/2;
            if xy(1) >= .5 && xy(2) >=.5,
                pos.x = w*shfac-d1*shfac2; pos.y = w*shfac-d1*shfac2;
            elseif xy(1) < .5 && xy(2) >= .5,
                pos.x = -w*shfac-d1+d1*shfac2; pos.y = w*shfac-d1*shfac2;
            elseif xy(1) < .5 && xy(2) < .5,
                pos.x = -w*shfac-d1+d1*shfac2; pos.y = -w*shfac-d1+d1*shfac2;
            elseif xy(1) >= .5 && xy(2) < .5,
                pos.x = w*shfac-d1*shfac2; pos.y = -w*shfac-d1+d1*shfac2;
            end

            if circ_selfprob,
                rectangle('Position',[xy(1)+pos.x,xy(2)+pos.y,d1,d1],...
                    'Curvature',[1,1],'EdgeColor',co,'LineWidth',3);
            end
        end
        
        rectangle('Position',[xy(1)-w/2,xy(2)-w/2,w,w],...
            'Curvature',[1,1],'FaceColor',col(:,i),'EdgeColor',colt(:,i),'LineWidth',.5)
        shift = []; shift.x = -.04; shift.y = 0.04;
        if colsum,
            if xy(1) < .5
                text(xy(1)-w/2+shift.x,xy(2)+shift.y,fieldn1{i},'HorizontalAlignment','Right','Fontsize',fntsize);
%                 if ~simplify,
%                     text(xy(1)-w/2+shift.x,xy(2)-shift.y,num2str(ceil(colsum)),'HorizontalAlignment','Right','Fontsize',fntsize);
%                 end
            else
                text(xy(1)+w/2-shift.x,xy(2)+shift.y,fieldn1{i},'HorizontalAlignment','Left','Fontsize',fntsize);
%                 if ~simplify,
%                     text(xy(1)+w/2-shift.x,xy(2)-shift.y,num2str(ceil(colsum)),'HorizontalAlignment','Left','Fontsize',fntsize);
%                 end
            end
%             if simplify, 
                text(xy(1),xy(2),num2str(ceil(colsum)),...
                  'HorizontalAlignment','Center','Fontsize',fntsize,'Color',coltxt(:,i)); 
%             end
        end
    end
end
% PLOT TRANSITION PROBABILITIES
if bool_val,
    for i=1:icnt,
        text(txt.x(i),txt.y(i),...
            num2str(txt.p(i),'%4.2f'),'Color','k',... %txt.c(:,i),...
            'HorizontalAlignment','Center','Fontsize',fntsize);
    end
end
axis([minmax minmax]); axis square; axis off;

% PLOT LEGEND (reference arrow with p = 0.5)
xy0 = [1 0]; dl = .3;
th = [0.5]/norm_arrow;
for i=1:length(th),
    dy = (i-1)*.2;
    patch([xy0(1) xy0(1)+dl xy0(1)+dl xy0(1) xy0(1)],...
          [xy0(2) xy0(2) xy0(2)-th(i)*2 xy0(2)-th(i)*2 xy0(2)]-dy,[.5 .5 .5],'EdgeColor','none'); hold on;
    text(xy0(1)+dl/2, xy0(2)-dy-th(i), num2str(th(i)*norm_arrow),'Color','w','HorizontalAlignment','Center','Fontsize',fntsize);
end

