%% OpenAVI
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
%% In the interactive mode, this function opens a dialog box and allows 
%% the user to specify the movie file(s) to process. In the non-interactive
%% mode, it simply prepares file paths for further processing. 

function openavi( interactive )
global params Files NFiles strInVideoFNameArray;

%%
% Open a file open dialog and let the user selects one or more .wmv files. 
% If user does not click cancel or abort the file selection dialog,
% then params.findex is set. Remove trailing slash in the path name, 
% Extract the extension from the first entry of the filename array,
% and count the number of file names in the array. 

if interactive,
[strInVideoFNameArray, Files.strInVideoPath, params.findex] = ...
    uigetfile( { '*.avi;*.wmv;*.mov;*.mpeg' , 'Movie file (*.avi,*.wmv,*.mov,*.mpeg)' ; ...
                 '*.*' , 'All files'}, ... 
                 'Open movie file', 'MultiSelect', 'on' );
end

if params.findex,
    Files.strInVideoPath = Files.strInVideoPath(1:end-1);
    if iscell(strInVideoFNameArray), 
        tmp = cell2mat(strInVideoFNameArray(1)); 
        NFiles = numel(strInVideoFNameArray);
    else
        tmp = strInVideoFNameArray;
        NFiles = 1;
    end
    Files.strVideoFExt = tmp(end-2:end);
    params.InPause = 0;
    params.filechg = 1; 
end
return