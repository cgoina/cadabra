%% Write_html
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
%% Creates webpages from given video clips and displays everything 
%% in the standard web browser of the system

% WEBPAGE CREATOR FOR MOVIE CLIPS OF DETECTED ACTIONS
% WEBPAGE IS CALLED BY LOCAL STANDARD WEB BROWSER
function write_html(path)
% PARAMETERS
imgs_per_row = 10; % define the number of clips per row in browser
clip_path = 'clips';
fname_main = 'clips_main';

% MOVIES ARE STORED BY ACTION AND BY GENOTYPE
% CREATE A MAIN WEBPAGE LINKING THE ACTION FOLDERS
% CREATE ONE WEBPAGE PER ACTION
% CREATE A WEBPAGE FOR EACH GENOTYPE PER ACTION
% LINK ALL MOVIES FOR EACH GENOTYPE

path0 = [path clip_path '/'];
scl = dir(path0);
if ~isempty(scl),
    % Create main web page
    fid0 = fopen([path0 fname_main '.html'],'w');
    print_header(fid0,'Clips of Actions');
    for icl=1:length(scl),
        % Linking of action folders
        path = [scl(icl).name '/'];
        in = strfind(path,'clips');
        if ~isempty(in) && isempty(strfind(path,'main')),
            fname = path(in+6:end-1);
            fprintf(fid0, ['<a href=\''clips_' fname '/' fname '.html\''>' fname '</a>\n<p>\n']);

            % Create action webpage
            fid1 = fopen([path0 path fname '.html'],'w');
            print_header(fid1,fname);
            s = dir([path0 path]);
            for i=1:length(s),
                % Linking of genotype folders
                if s(i).isdir,
                    f = dir([path0 path s(i).name '/*.gif']);
                    if ~isempty(f),
                        fprintf(fid1, ['<a href=\''' s(i).name '/' s(i).name '.html\''>' s(i).name '</a>\n<p>\n']);
                        
                        % Create genotype webpage
                        fid2 = fopen([path0 path s(i).name '/' s(i).name '.html'],'w');
                        print_header(fid2,s(i).name);
                        icnt = 0;
                        for ii=1:length(f),
                            % Linking of movies
                            if ~mod(icnt,imgs_per_row),
                                fprintf(fid2,'<br>\n');
                            end
                            fprintf(fid2, ['<img src=\''' f(ii).name '\'' />\n']);
                            icnt = icnt + 1;
                        end
                        fclose(fid2);
                    end
                end
            end
            fprintf(fid1, '</body>\n'); fprintf(fid1, '</html>');
            fclose(fid1);
        end
    end
    fclose(fid0);
    
    % CALL WEBBROWSER WITH MAIN WEBPAGE
    web([path0 fname_main '.html'],'-browser');
end


function print_header(fid,fname)
fprintf(fid, '<!DOCTYPE HTML PUBLIC \''-//W3C//DTD HTML 4.01 Transitional//EN\''>\n');
fprintf(fid, '<html lang=\''en\''>\n');
fprintf(fid, '<head>\n');
fprintf(fid, ['<title>' fname '</title>\n']);
fprintf(fid, '</head>\n');
fprintf(fid, '<body>\n');
fprintf(fid, [fname '\n']);
fprintf(fid, '<p><p>\n');
