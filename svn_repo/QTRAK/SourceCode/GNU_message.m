%% GNU message generator
% Copyright (C) 2009 JFRC/HHMI, Caltech

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
function msg = GNU_message(short_msg)

if nargin < 1,
    short_msg = 0;
end

if ~short_msg,
    msg = cell(19,1);
    msg{1} = 'CADABRA v1.1 (QTRAK v1.1)';
    
    msg{3} = 'Copyright (c) 2009, Heiko Dankert    Jinyang Liu, Caltech & JFRC/HHMI';

    msg{5} = 'This program is part of QTRAK and the "Caltech Automated Drosophila';
    msg{6} = 'Aggression-Courtship Behavioral Repertoire Analysis (CADABRA)".';
    
    msg{8} = 'This program is free software: you can redistribute it and/or modify';
    msg{9} = 'it under the terms of the GNU General Public License as published by';
    msg{10} = 'the Free Software Foundation, either version 3 of the License, or';
    msg{11} = '(at your option) any later version.';

    msg{13} = 'This program is distributed in the hope that it will be useful,';
    msg{14} = 'but WITHOUT ANY WARRANTY; without even the implied warranty of';
    msg{15} = 'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.';
    msg{16} = 'See the GNU General Public License for more details.';

    msg{18} = 'You should have received a copy of the GNU General Public License';
    msg{19} = 'along with this program.  If not, see <http://www.gnu.org/licenses/>.';
else
    msg = cell(4,1);
    msg{1} = 'Copyright © 2009 , Heiko Dankert    Jinyang Liu, Caltech & JFRC/HHMI';
    msg{3} = 'You should have received a copy of the GNU General Public License';
    msg{4} = 'along with this program.  If not, see <http://www.gnu.org/licenses/>.';
end