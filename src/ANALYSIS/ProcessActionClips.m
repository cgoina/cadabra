%% ProcessActionClips
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
%% "DDGrab" / "FFGrab" and pass them to "write_clip"

function ProcessActionClips(data,width,height,frameNr,time)
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
      
            scanline = ceil(width*3/4)*4; % the scanline size must be a multiple of 4.
            
            % some times the amount of data doesn't match exactly what we expect...
            if (numel(data) ~= scanline*height)
                if (numel(data) > 3*width*height)
                    if (isemtpy(warned))
                        warning('mmread:general','dimensions do not match data size. Guessing badly...');
                        warned = true;
                    end
                    scanline = width*3;
                    data = data(1:3*width*height);
                else
                    error('dimensions do not match data size. Too little data.');
                end
            end
            
            % if there is any extra scanline data, remove it
            data = reshape(data,scanline,height);
            data = data(1:3*width,:);
            
            % the data ordering is wrong for matlab images, so permute it
            data = permute(reshape(data, 3, width, height),[3 2 1]);
            
            % CONVERT FRAME INTO GRAYSCALE
            %data = uint8((0.11*data(:,:,1) + 0.59*data(:,:,2) + 0.3*data(:,:,3))*255); %EH100423
            
            % CALL "write_clip"
            [indx,indxg,b,ImageBuffer] = write_clip(data,dframe,start_frm,end_frm,position,timeindex,indx,indxg,b,ImageBuffer,addName,cluster,frameNr,onefly,gen);
        end
    catch
        disp(lasterr);
    end
end