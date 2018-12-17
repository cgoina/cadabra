%% Callback
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
function callback
global params

    switch get(gcbo,'Label')
        case 'Open'
            openavi
        case 'Quit [Esc]'
            params.QuitPrg = 1;
        case 'Show Tracking'
            if strcmp(get(gcbo, 'Checked'),'on')
                set(gcbo, 'Checked', 'off');
                params.bool_plottrack = 0; 
            else 
                set(gcbo, 'Checked', 'on');
                params.bool_plottrack = 1; 
            end
        case 'Show Counting'
            if strcmp(get(gcbo, 'Checked'),'on')
                set(gcbo, 'Checked', 'off');
                params.bool_plotcount = 0; 
            else 
                set(gcbo, 'Checked', 'on');
                params.bool_plotcount = 1; 
            end
        case 'Show Bounding Box'
            if strcmp(get(gcbo, 'Checked'),'on')
                set(gcbo, 'Checked', 'off');
                params.bool_boundbox = 0; 
            else 
                set(gcbo, 'Checked', 'on');
                params.bool_boundbox = 1; 
            end
        case 'About'
            msg = GNU_message;
            h = msgbox(msg);
            uiwait(h);
    end

return
