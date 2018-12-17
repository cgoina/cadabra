%% ProcessFrameCapWin
% Copyright (C) 2008 Heiko Dankert, California Institute of Technology

% This file is part of QTRAK.

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
%%
% This is the Matlab procedure called from the C routine. 
% The input parameter data is a 24-bit RGB image with the 8-bit (R,G,B) 
% values arranged row-wise. The image coordinate system is such that 
% the bottom left corner is (0,0). 

function ProcessFrameCapWin(data, images, width, height, frameNr, time)
global ic dt params CurrentK;
global Files FigureHandle Panels;

    if params.QuitPrg, error('STOP PROCESSING'); return; end
    
    try
 
        %% 
        % Convert the image into normalized double precision RGB image. 
        % Save the frame number into the global variable and calculatecc the 
        % cumulative average of interframe separation.
        
        CurrentK = frameNr;
        dt = time / frameNr;
        

        %% 
        % If requested, update the frame and time counter, and/or 
        % update the tracking markers on the screen. Note that this 
        % feature incurs a heavy performance penalty (relative to its 
        % function) and should be used only when needed. 
        
        if params.bool_plotcount, 
            minu = floor((frameNr-2)*dt/60);
            sec = floor((frameNr-2)*dt-minu*60);
            hun=floor(((frameNr-2)*dt-minu*60-sec)*100);
            figure(FigureHandle);
            axes(Panels.HInfo);
            cla(Panels.HInfo);
            image(Panels.BlankInfoPanel);
            set(Panels.HInfo,'XTick',[],'YTick',[]);
            text(5,15,['Frame : ' num2str((frameNr-2),'%06d')],'FontSize',7,'color','w');
            text(80,15,['Time : ' num2str(minu,'%02d') ' : ' num2str(sec,'%02d') ... 
                ' : ' num2str(hun,'%02d')],'FontSize',7,'color','w');
        end
        
        %% 
        % Finally, analyze the input image by starting the outer and
        % inner loop of the analysis function. 
        
        PlotImage_Captioning( data, images, width, height, FigureHandle, frameNr );            

    catch err
        fprintf(Files.ErrorFID{ic}, '\nFrameNumber: %06d %s', frameNr, err.message);
        rethrow(err);
    end

end