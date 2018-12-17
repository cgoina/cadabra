%% Etho_matrix_plot
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
%% Plot transition matrix

% PLOT TRANSITION MATRIX
function etho_matrix_plot(h,lab,max_v,fnts)

% DATA VALUE SCALE
if min(min(h)) < 0, 
    min_v = -max_v;
    col = hot; for j=1:4, col(ceil(length(col)/2)+j-2,:) = [0 0 0]; end
else
    min_v = 0; 
    col = hot;
end
num_beh = length(lab);
colormap(col); 

% PLOT MATRIX
imagesc(h,[min_v max_v]); 

% PLOT OCCURRENCE VALUES INTO MATRIX FIELDS
h(abs(h)>=max_v) = 0; ind = find(abs(h) > 0); 
txtclr = hot(max_v); xi = repmat(1:num_beh,num_beh,1); yi = xi';
for i=1:numel(ind),
    val = h(ind(i)); dx = .32;
    text(xi(ind(i))-dx,yi(ind(i))+.06,num2str(val,'%.1f'),...
        'Color',txtclr(ceil(max_v-abs(h(ind(i)))),:),'FontSize',fnts);
end
for i=1:length(lab),
    in = strfind(lab{i},' ');
    lab{i} = lab{i}(1:in(end)-1);
end
axis square; set(gca,'XTick', 1:num_beh); xticklabel_rotate([1:num_beh],45,lab,'interpreter','none');
set(gca,'YTick', 1:num_beh); set(gca,'YTickLabel',lab);
set(gca,'FontSize',fnts);
ch = colorbar('location','EastOutside'); set(get(ch,'YLabel'),'String','occurence [-]')
