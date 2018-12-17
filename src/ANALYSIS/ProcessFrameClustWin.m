%% ProcessFrameClustWin
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
%% Routine called by 'mmread' to process frames that were read with
%% 'mexDDGrab' and them to 'write_clip'

function ProcessFrameClustWin(data,width,height,frameNr,time)
global ImageBuffer;
global dframe;
global start_frm;
global end_frm;
global position;
global timeindex;
global indx;
global indxg;
global b;
global addName;
global cluster;
global onefly;
global gen;
global h1;

persistent warned;

if nargin,
    try
        if (mod(frameNr,1000) == 0), 
            waitbar((frameNr-start_frm(1))/(end_frm(end)-start_frm(1)),h1);
        end
        if (indx <= length(start_frm)) && (frameNr >= start_frm(indx)),
            nelem = width*height*3;
            data = 0.11*data(1:3:nelem) + 0.59*data(2:3:nelem) + 0.3*data(3:3:nelem);
            data = rot90(reshape(data, width, height));
            [indx,indxg,b,ImageBuffer] = write_clip(data,dframe,start_frm,end_frm,position,timeindex,indx,indxg,b,ImageBuffer,addName,cluster,frameNr,onefly,gen);
        end
    catch
        disp(lasterr);
    end
end