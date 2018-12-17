%% writeClips
% writeClips grabbs some users specified frames from a moive and save the clips to avi files. 
% 
% Usage:
%    writeClips(movie, frames, prePost, trakcing)
%         movie is the full name of a moive to be processed
%         frames is a frame array which include all the user specified frame number
%         prePost is a two element array [preFrame postFrame]. PreFrame is 
%            the number of frames to be saved before each user sepecified frame number and
%            postFrame is the number of the frames to be saved after it. 
%         tracking is logical value. When tracking equals one, the save clips 
%            include the tracking informaiton, it not, then the trakcing information 
%            is not included in the clips. 
% Example:  
%      writeClips('U:\SourceCode\oneChamberMoive\test7_p1.avi', [5 100], [5 15], 0);

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

% Copyright (C) 2009 Jinyang Liu, Janelia Farm Research Campus
%
%%
% This is the main entry point of the entire tracking application. Most
% actions performed in this module is self-explanatory, and are explained
% in inline comments. The function makes heavy use of global variables, a
% feature that can be revised in the next iteration of this code.

function writeClips(movie, frames, prePost, tracking)

[pathstr, name, ext] = fileparts(movie); 
flyMovObj = mmreader(movie);
clipLength = prePost(2) + prePost(1) + 1;
%nFrames = flyMovObj.NumberOfFrames;
vidHeight = flyMovObj.Height;
vidWidth = flyMovObj.Width;

hfig = figure('Visible', 'off');
axes ('Visible', 'off','Units', 'normalized', 'Position', [0 0 1 1]);

for frameCount = 1:length(frames)
    
    % Preallocate movie structure.
    mov(1:clipLength) = ...
        struct('cdata', zeros(vidHeight, vidWidth, 3, 'uint8'),...
        'colormap', []);
    
    avifileName = [name '_' num2str(frames(frameCount)) '.avi'];
    aviFullFileName = fullfile(pathstr, avifileName);
    aviobj = avifile(aviFullFileName,'fps',flyMovObj.FrameRate/2);
    
    % Read one frame at a time.
    for currentFrame = 1 : clipLength
        frameNum = frames(frameCount) - prePost(1) + currentFrame -1;
        mov(currentFrame).cdata = read(flyMovObj, frameNum);
        image(mov(currentFrame).cdata);
        axis image;

        textInfo = ['FrameNum: ' num2str(frameNum)];
        text(10,10, textInfo);
        mov(currentFrame) = capturescreen(gcf);
    end
    
    aviobj = addframe(aviobj,mov);
    aviobj = close(aviobj);
    
end

close(hfig);
