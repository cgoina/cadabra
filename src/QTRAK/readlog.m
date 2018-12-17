%% Readlog
% Copyright (C) 2008 Heiko Dankert, California Institute of Technology

% This file is part of QTRAK.
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
%%
% Read log information about finished jobs, read actual file list, 
% remove finished jobs from list

function [strInVideoFNameArray,NFiles,Files] = readlog(infofile,finished)

logname = 'jobs';
if ispc, slashstr = '\'; else slashstr = '/'; end

fileext = '';
if isdir(infofile), 
    basepath = infofile; 
else
    indslsh = strfind(infofile,slashstr);
    if numel(indslsh),
        indslsh = indslsh(end);
        basepath = infofile(1:indslsh);
        fileext = infofile(indslsh+1:end);
    else
        basepath = './';
    end
end

% Search for all existing movies in basepath
if ~numel(fileext),
    files = [dir([basepath '*.avi']) ; dir([basepath '*.wmv'])];
else
    files = dir(infofile) ;
end

if ~isempty(files),
    
    % Update log file in case a file was completely processed
    if finished,
        % Read current log file
        [strInVideoFNameArray_compressed, proc] = read_logfile([basepath logname '.log']);
        
        % Update processing index and write back log file
        swap = 0;
        fid = fopen([basepath logname '.log'],'w');
        for i=1:numel(strInVideoFNameArray_compressed),
            if ~proc(i) && ~swap, proc(i) = 1; swap = 1; end
            fprintf(fid,'%40s %1g\n',strInVideoFNameArray_compressed{i},proc(i));
        end
        fclose(fid);
    end
    
    % Extract file names from files-structure and remove white spaces
    NFiles = numel(files);
    strInVideoFNameArray = cell(NFiles,1);
    strInVideoFNameArray_compressed = cell(NFiles,1); % no white spaces
    for i=1:NFiles,
        ind = setdiff(1:length(files(i).name),strfind(files(i).name,' '));
        strInVideoFNameArray_compressed{i} = files(i).name(ind);
        strInVideoFNameArray{i} = files(i).name;
    end

    % Search for existing log file
    logfname = [basepath '*.log'];
    logfile = dir(logfname); 
    
    if isempty(logfile), logfile(1).bytes = 0; end
    if ~logfile.bytes,
        % If log file does not exists, create a new one
        % Set the 'completely processed index' to Zero
        fid = fopen([basepath logname '.log'],'w');
        for i=1:NFiles,
            fprintf(fid,'%40s %1g\n',strInVideoFNameArray_compressed{i},0);
        end
        fclose(fid);
    else
        % Read log file
        [strInVideoFNameArray_compressed, proc] = read_logfile([basepath logname '.log']);

        % If a file already finished exclude it from the todo list
        ind = find(~proc);
        if numel(ind),
            strInVideoFNameArray = strInVideoFNameArray(ind);
        else
            strInVideoFNameArray = [];
        end

        NFiles = numel(strInVideoFNameArray);
    end
    
    if NFiles,
        % Locate ROI/calibration file
        roi = 0;
        roifile = dir([basepath strInVideoFNameArray{1}(1:end-4) '*_1_roi.mat']);
        if isempty(roifile), roifile(1).bytes = 0; end
        if roifile(1).bytes,
            roi = 1;
        else
            % If there is none we need take the ROI/calibration infos
            % from another movie, assuming ROI/calibration conditions
            % within the selected files stay constant
            roifile = dir([basepath '*_1_roi.mat']);
            if isempty(roifile), roifile(1).bytes = 0; end
            if roifile(1).bytes,
                roi = 1;
                % Determine number of chambers (sourcefiles to copy)
                sourcefiles = dir([basepath roifile(1).name(1:end-9) '*_roi.mat']);
                % Copy ROI/calibration info to first movie in list
                for ii=1:numel(sourcefiles),
                    src = [basepath sourcefiles(ii).name];
                    tgt = [basepath strInVideoFNameArray{1}(1:end-4) '_' ...
                        num2str(ii) '_roi.mat'];
                    copyfile(src,tgt);
                end
            end
        end
        if ~roi,
            fprintf('\nNo calibrated data.\nPlease calibrate on movie first.\n');
            NFiles = 0;
        end
        Files.strInVideoPath = basepath(1:end-1);
        Files.strVideoFExt = strInVideoFNameArray{1}(end-2:end);
        if NFiles == 1, strInVideoFNameArray = cell2mat(strInVideoFNameArray); end
    else
        strInVideoFNameArray = 0;
        NFiles = 0;
        Files = 0;
    end
else
    strInVideoFNameArray = 0;
    NFiles = 0;
    Files = 0;
end


function [strInVideoFNameArray, proc] = read_logfile(logfile)

NFiles = 0;
fid = fopen(logfile,'r');
while ~feof(fid),
    Line = textscan(fid,'%40s %1s',1);
    if ~isempty(Line{1}), NFiles = NFiles + 1; end
end
fclose(fid);

strInVideoFNameArray = cell(NFiles,1); proc = zeros(NFiles,1);

NFiles = 0;
fid = fopen(logfile,'r');
while ~feof(fid),
    Line = textscan(fid,'%40s %1s',1);
    if ~isempty(Line{1}), 
        NFiles = NFiles + 1;
        strInVideoFNameArray{NFiles} = cell2mat(Line{1});
        proc(NFiles) = str2double(Line{2});
    end
end
fclose(fid);
