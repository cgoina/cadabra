%% Set_xtick_label
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
%% rotated xtick labels

% PRINT THE XTICK LABELS AT AN ANGLE INSTEAD OF HORIZONTALLY
function set_xtick_label(tick_labels, angle, axis_label,axisfontsize)
% SET_XTICK_LABEL Print the xtick labels at an angle instead of horizontally
% set_xtick_label(tick_labels, angle, axis_label)
%
% angle default = 90
% axis_label default = ''
%
% This is derived from Solution Number: 5375 on mathworks.com
% See set_xtick_label_demo for an example

if nargin < 2, angle = 90; end
if nargin < 3, axis_label = []; end

% Reduce the size of the axis so that all the labels fit in the figure.
pos = get(gca,'Position');
set(gca,'Position',[pos(1), .2, pos(3) .65])

ax = axis; % Current axis limits
axis(axis); % Fix the axis limits
Yl = ax(3:4); % Y-axis limits

%set(gca, 'xtick', 1:length(tick_labels));
set(gca, 'xtick', 0.7:1:length(tick_labels));
Xt = get(gca, 'xtick');

% Place the text labels
if nargin < 4, 
    t = text(Xt,Yl(1)*ones(1,length(Xt)),tick_labels);
else
    t = text(Xt,Yl(1)*ones(1,length(Xt)),tick_labels,'FontSize',axisfontsize);
end
set(t,'HorizontalAlignment','right','VerticalAlignment','top', 'Rotation', angle);

% Remove the default labels
set(gca,'XTickLabel','')

% Get the Extent of each text object. This
% loop is unavoidable.
for i = 1:length(t)
ext(i,:) = get(t(i),'Extent');
end

% Determine the lowest point. The X-label will be
% placed so that the top is aligned with this point.
LowYPoint = min(ext(:,2));

% Place the axis label at this point
if ~isempty(axis_label)
Xl = get(gca, 'Xlim');
XMidPoint = Xl(1)+abs(diff(Xl))/2;
tl = text(XMidPoint,LowYPoint, axis_label, 'VerticalAlignment','top', ...
'HorizontalAlignment','center');

end