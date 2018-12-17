%% Condens_bouts
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
%% Preprocessing of paths and file names provided by Analysis GUI
%% Feature file names follow a naming convention for the tested genotypes;
%% furthermore the GUI offers to take to take the data sub
%% folder name as genotype naming; otherwise the user may assign 
%% fly pairs to genotypes manually

function [sorted_names, fnames, mov_count, gen_ident, genotypes, paths] = filelist(path,bool_sub_dirs,bool_screen)
%%% bool_sub_dirs ... files in sub-folders?
%%% bool_screen   ... = 0: take sub-folder name as genotype, othewise from file-names

sorted_names = cell(1);
fnames = cell(1);
paths = cell(1);

if bool_screen,
    if bool_sub_dirs,
        subdirs = dir(path);
    else
        subdirs.name = '.';
    end    
    icnt = 0; lrcnt = 0;
    for i=1:length(subdirs),
        if (~strcmp(subdirs(i).name,'.') && ~strcmp(subdirs(i).name,'..')) || ~bool_sub_dirs,
            files = dir([path subdirs(i).name '/*.feat']);
            for j=1:length(files),
                icnt = icnt + 1;
                fname = files(j).name; ind = strfind(fname,'_'); 
                strt = ind(1)+1; strt2 = ind(2)-1;
                L = strfind(fname,'L'); R = strfind(fname,'R');
                S = strfind(fname,'ingle'); C = strfind(fname,'chamber');
                S1 = []; s1 = 0;
                if numel(L) && numel(R),
                    L = L(1); R = R(1);
                    if ~lrcnt,
                        fnames{icnt} = fname(L+1:R-1);
                    else
                        l = strfind(fname(R:end),'_'); l = R + l(1)-1;
                        fnames{icnt} = fname(R+1:l);
                    end
                    lrcnt = lrcnt + 1; if lrcnt > 1, lrcnt = 0; end
                elseif numel(S),
                    fnames{icnt} = fname(strt:S-2);
                elseif numel(S1),
                    fnames{icnt} = fname(strt:S1-2);
                elseif numel(C),
                    fnames{icnt} = fname(strt:end-14);
                else
                    fnames{icnt} = fname(strt:strt2);
                end
                if length(fnames{icnt}) > 1,
                    if strcmp(fnames{icnt}(end-1),'_'), fnames{icnt} = fnames{icnt}(1:end-2); end
                    if strcmp(fnames{icnt}(end),'_'),
                        fnames{icnt} = fnames{icnt}(1:end-1);
%                     elseif ~strcmp(fnames{icnt}(end),'m') && ~strcmp(fnames{icnt}(end),'-'),
%                         fnames{icnt} = fnames{icnt}(1:end-3);
                    end
                end
                if ~strcmp(fnames{icnt},'none'),
%                     if strfind(fnames{icnt},'w-'), fnames{icnt} = 'w-Exel'; end
                    if strfind(fnames{icnt},'CS'), fnames{icnt} = 'CS-Tully'; end
%                     if strfind(fnames{icnt},'Kir'), fnames{icnt} = fnames{icnt}(1:end-4); end
    %                 if strcmp(fnames{icnt}(end),'_'), fnames{icnt} = fnames{icnt}(1:end-1); end
                    paths{icnt} = [path subdirs(i).name '/'];
                    sorted_names{icnt} = [paths{icnt} fname];
                else
                    icnt = icnt - 1;
                end
            end
        end
    end
    gen_ident = unique(fnames); nident = length(gen_ident); 
    mov_count = zeros(nident,1); genotypes = zeros(length(fnames),1);
    for i=1:length(fnames),
        for j=1:nident,
            if strcmp(fnames{i},gen_ident{j}), 
                genotypes(i) = j;
                mov_count(j) = mov_count(j) + 1;
            end
        end
    end    
else
    if bool_sub_dirs,
        dirs = dir(path);
        ndirs = length(dirs);
    else
        ndirs = 1;
        dirs(4) = '.';
    end
    icnt = 0; mov_count = [];
    for idir=1:ndirs,
        if ~bool_sub_dirs || (bool_sub_dirs && ~strcmp(dirs(idir).name,'.') && ~strcmp(dirs(idir).name,'..')),
            if bool_sub_dirs,
                path1 = [path dirs(idir).name filesep];
            else
                path1 = path;
            end
            dir_struct = dir(fullfile(path1,'*.feat'));
            if numel(dir_struct),
                sorted_names1 = sortrows({dir_struct.name}');
                for i=1:length(sorted_names1),
                    if numel(sorted_names1{i}),
                        icnt = icnt + 1;
                        sorted_names{icnt} = [path1 sorted_names1{i}];
                        in = strfind(path1,'/'); ident = path1(in(end-1)+1:in(end)-1);
                        if strcmp(ident,'f02508'), ident = 'Sr'; end
                        if strcmp(ident,'f03827'), ident = 'CG9733'; end
                        if strcmp(ident,'f03446'), ident = 'CG14186'; end
                        if strcmp(ident,'e00517'), ident = 'Hr38'; end
                        if strcmp(ident,'d06417'), ident = 'CG15745'; end
                        fnames{icnt} = ident;
                        paths{icnt} = path1;
                    end
                end
                if length(sorted_names1), mov_count = [mov_count length(sorted_names1)]; end
            end
        end
    end   
    gen_ident = unique(fnames); genotypes = zeros(length(fnames),1);
    for i=1:length(gen_ident), genotypes(strcmp(fnames,gen_ident(i))) = i; end
end
end