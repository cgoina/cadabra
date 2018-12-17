function processFrame( data, images, width, height, frameNr, time )
% processFrame(data,width,height,frameNr)
%
% This is the function prototype to be used by the matlabCommand option of
% mmread.
% INPUT
%   data        the raw captured frame data, the code below will put it
%               into a more usable form
%   width       the width of the image
%   height      the height of the image
%   frameNr     the frame # (counting starts at frame 1)
%   time        the time stamp of the frame (in seconds)
%
% Warning, the way that this is implemented requires buffering about 1
% seconds worth of video.  If you don't have enough memory in your system
% to hold 1 second's worth of video, this will crash on you.  Either make
% your movies smaller, or buy more memory.
% 
% EXAMPLES
%   Process all frames in a movie using this function:
%   mmread('mymovie.mpg',[],[],false,false,'processFrame');
%
%   Process only frames 1 through 10 using this function:
%   mmread('mymovie.mpg',1:10,[],false,false,'processFrame');
% 
% Copyright 2008 Micah Richert
% 
% This file is part of mmread.
% 
% mmread is free software; you can redistribute it and/or modify it
% under the terms of the GNU Lesser General Public License as
% published by the Free Software Foundation; either version 3 of
% the License, or (at your option) any later version.
% 
% mmread is distributed WITHOUT ANY WARRANTY.  See the GNU
% Lesser Lesser General Public License for more details.
% 
% You should have received a copy of the GNU Lesser General Public
% License along with this program.  If not, see
% <http://www.gnu.org/licenses/>.

persistent warned;

try
    % now do something with the data...
    image(data);
    title(['frame ' num2str(frameNr) ' ' num2str(time) 's ']);
    drawnow;

% stop early
%     if (frameNr == 10)
%         error('processFrame:STOP','STOP!!!');
%     end
    
catch
    % if we don't catch the error here we will loss the stack trace.
    err = lasterror;
    disp([err.identifier ': ' err.message]);
    for i=1:length(err.stack)
        disp([err.stack(i).file ' (' err.stack(i).name ') Line: ' num2str(err.stack(i).line)]);
    end
    rethrow(err);
end
