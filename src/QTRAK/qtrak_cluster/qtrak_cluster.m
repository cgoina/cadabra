%% qtrak_cluster: cluster version of the qtrak
%
% Usage:  qtrak_cluster
%         qtrak_cluster fullfilename
%
%
% fullfilename is the full file name of the preprocessd.config which
% include the absolute path. If no fullfilename followed, a pop-up window
% will be used for user to locate the preprocessed.config.
%
% qtrak_cluster is derived from qtrak.m. In order to be adapted for the
% cluster computing, qtrak.m was modifed to two parts: qtrak_preprocess and
% qtrak_cluster. qtrak_cluster inludes all the GUI for user to input
% experimental parameters. qtrak_cluster only includes the computing part.
% qtrak_cluster need to load a configure file which is generated from the
% qtrak_preprocess. The configure file can be input as a function arguments
% or use a pop-up window to for user to locate it.
%
% * Original implementation by Heiko Dankert, modified by Jinyang Liu
% * Optimization and documentation by Edwin Soedarmadji

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
% * Original implementation by Heiko Dankert, modified by Jinyang Liu
% * Optimization and documentation by Edwin Soedarmadji
%
%%
% This is the main entry point of the entire tracking application. Most
% actions performed in this module is self-explanatory, and are explained
% in inline comments. The function makes heavy use of global variables, a
% feature that can be revised in the next iteration of this code.
function qtrak_cluster(inputFile)

%#function ProcessFrameCapWin
%#function PlotImage_Captioning
%#function histfit
%#function fhist gmmactiv gmmpost
%#function FFGrab seg_fullfly_inflexbox
%#function fhistc ellipsefit write_sequence
%#function gauss write_sequence_crop callback
%#function seg_fullfly_inconstbox segfly fotsu
%#function intersect_flies fly_bodymeas switch_flies
%#function fly_headtailwings mask2ellipses gmm gmminit consist
%#function kmeans dist2 gmmem gmmprob histfit
%#function auto_cal auto_cal_circ CircularHough_Grd fhistc openavi
%#function calibrate calibrate fhistfit GNU_message

% global params chamber scale;
% global slashstr Files NFiles strInVideoFNameArray;
% global DimX DimY DimZ Spacer;
% global intStartFrm intNFrms dt;
% global nframes_mean mean_image std_image frmindx;
global object object_1 object_2;
global currentFrame_O1 currentFrame_O2;

% Per Eric request, the moive might be moved to archieve, the inputfile can
% be either a moive file or a configuration file
if nargin
    [inputFilePath, inputFileName] = fileparts(inputFile);
    inputFilePath = [inputFilePath '/'];
    
else
    [inputFileNameWithExt, inputFilePath, inputFileIndex] = ...
        uigetfile('*.avi;*.wmv;*.config','choose a file to be processed' );
    inputFile = fullfile(inputFilePath, inputFileNameWithExt);
    [temp, inputFileName]=fileparts(inputFile);
    
    if ~inputFileIndex
        errordlg('A specified file should be chosen','File Error');
        return;
    end
end

tic;

%no matter input file is a moive file or a configuration file, we need the variable for configFile, featFile, and inputMoiveFile. 
configFileName = [inputFileName '.config'];
configFilePath = inputFilePath;
configFile=fullfile(configFilePath, configFileName);

featFileName = [inputFileName '_1.feat'];
featFile = fullfile(configFilePath, featFileName);

movieAVIName = [inputFileName, '.avi'];
movieWMVName = [inputFileName, '.wmv'];

movieAVI = fullfile(inputFilePath, movieAVIName);
movieWMV = fullfile(inputFilePath, movieWMVName);

if (exist(movieAVI, 'file') == 2)
    inputMovieFile = movieAVI;
elseif (exist(movieWMV, 'file') == 2)
    inputMovieFile = movieWMV;
end


% Check to see if updated config file exists
if (exist(configFile, 'file')==2)
    load(configFile, '-mat');
    if ~isfield(params,'analysis')
        %This code will update old config.params structure produced with revision 66 and earlier of qtrak_preprocess.
        params.analysis.radius0 = (chamber.ncols / 2) * scale.x;
        params.analysis.border_width = 2.0;                         % This is temporarily hard-coded, make sure to remove after all data is reanalyzed.
        params.analysis.radius = params.analysis.radius0 - params.analysis.border_width;
        params.analysis.tuningthreshold = 0.7;                       % Also temporarily hard-coded.
        params.analysis.correct_orient = 1;                          % temporarily hard-coded
        params.analysis.correct_positions = 1;                       % temporarily hard-coded
        params.analysis.max_frames = intNFrms;
        % save config file here
        save(configFile,'params', 'Files', 'NFiles',...
            'DimX', 'DimY', 'DimZ', 'Spacer', 'chamber',...
            'intStartFrm', 'intNFrms', 'dt', 'nframes_mean',...
            'frmindx', 'scale');
    end
else
    errordlg('A preprocessed.config is wanted! If you cannot find the file, please run qtrak_preprocess to generate it.','File Error');
    return
    
end


% Check to see if a feat file already exist.  If so, skip qtrak routine.
if ~(exist(featFile, 'file')==2)
    
    if ~(exist(inputMovieFile, 'file')==2)
        errordlg('Cannot find the feature file! Rerun qtrak_cluster from movie file.','File Error');
        return;
    end
    
    params.bool_plotcount = 0;
    ObjBuf = cell(params.nchambers,1);
    % Files.strInVideoFName = Files.strInVideoFName(1:end-4);
    % Files.strInFName = [Files.strInVideoPath slashstr Files.strInVideoFName '.' Files.strVideoFExt];
    
    FeatureFileName = cell(params.nchambers,1);
    ErrorFileName = cell(params.nchambers,1);
    Files.FeatureFID = cell(params.nchambers,1);
    Files.ErrorFID = cell(params.nchambers,1);
    
    for i=1:params.nchambers
        
        FeatureFileName{i} = [inputFilePath inputFileName '_' num2str(i) '.feat'];
        ErrorFileName{i} = [inputFilePath inputFileName '_' num2str(i) '.err'];
        
        Files.FeatureFID{i} = fopen(FeatureFileName{i},'w');
        if Files.FeatureFID{i} == -1,
            disp('ERROR: The Following Feature file could not be Created');
            disp(['         NAME: ', FeatureFileName{i}]);
            disp(['         PATH: ', Files.strInVideoPath slashstr]);
            return
        end
        
        fprintf( Files.FeatureFID{i}, '%14s', ...
            'frame'       ,'time [s]'    ,'fly1_dir'    ,'fly2_dir'    , ...
            'fly1_mvdir'  ,'fly2_mvdir'  ,'fly1_ori'    ,'fly2_ori'    , ...
            'dir_diff'    ,'mvdir_diff'  ,'1to2_mvdird' ,'2to1_mvdird' , ...
            'fly1_mean'   ,'fly2_mean'   ,'e_len_1 [mm]','e_len_2 [mm]', ...
            'e_ar1 [mm2]' ,'e_ar2 [mm2]' ,'vel1 [mm/s]' ,'vel2 [mm/s]' , ...
            'ac1 [mm/s^2]','ac2 [mm/s^2]','c1dist [mm]' ,'c2dist [mm]' , ...
            'dist [mm]'   ,'h1-t2 [mm]'  ,'h2-t1 [mm]'  ,'h-h [mm]'    , ...
            't-t [mm]'    ,'ddist [mm]'  ,'d h-h [mm]'  ,'d t-t [mm]'  , ...
            'Area1 [mm2]' ,'Area2 [mm2]' ,'Length1 [mm]','Length2 [mm]', ...
            'fly_closing' ,'fly_facing'  ,'fly_acting'  ,'fly_conn'    , ...
            'fly1_wingr'  ,'fly1_wingl'  ,'fly1_wing'   ,'fly2_wingr'  , ...
            'fly2_wingl'  ,'fly2_wing'   ,'chpos1 [mm]' ,'chpos2 [mm]' , ...
            'fly1 x [mm]' ,'fly1 y [mm]' ,'fly2 x [mm]' ,'fly2 y [mm]' , ...
            'phi1_r'      ,'phi1_l'      ,'fly1_r [mm]' ,'fly1_l [mm]' , ...
            'phi2_r'      ,'phi2_l'      ,'fly2_r [mm]' ,'fly2_l [mm]' );
        fprintf(Files.FeatureFID{i},'\n');
        
        Files.ErrorFID{i} = fopen(ErrorFileName{i},'a');
        if Files.ErrorFID{i} == -1
            disp('ERROR: The Following Error file could not be Created');
            disp(['         NAME: ', ErrorFileName{i}]);
            disp(['         PATH: ', Files.strInVideoPath slashstr]);
            return
        end
        
        % inintialize some global varibles
        ObjBuf{i} = struct;
        object{i}.dirdiff = 0;
        object{i}.mvdirdiff = 0;
        object{i}.distc = 0;
        object{i}.disth = 0;
        object{i}.distt = 0;
        object{i}.der_distc = 0;
        object{i}.der_disth = 0;
        object{i}.der_distt = 0;
        object_1{i}.disthto2 = 0;
        object_1{i}.pos_x = zeros(1,intNFrms);
        object_1{i}.pos_y = zeros(1,intNFrms);
        object_1{i}.headdir = zeros(1,intNFrms);
        object_1{i}.movedir = 0;
        object_1{i}.to2mvdirdiff = 0;
        object_1{i}.orient = 0;
        object_1{i}.mea = 0;
        object_1{i}.vel = 0;
        object_1{i}.acc = 0;
        object_1{i}.area = 0;
        object_1{i}.Area = 0;
        object_1{i}.length = 0;
        object_1{i}.r = 0;
        object_2{i} = object_1{i};
        
        currentFrame_O1{i}.xc = zeros(1,intNFrms);
        currentFrame_O1{i}.yc = zeros(1,intNFrms);
        currentFrame_O1{i}.head = zeros(1,intNFrms);
        currentFrame_O1{i}.xh = zeros(1,intNFrms);
        currentFrame_O1{i}.xt = zeros(1,intNFrms);
        currentFrame_O1{i}.yh = zeros(1,intNFrms);
        currentFrame_O1{i}.yt = zeros(1,intNFrms);
        currentFrame_O1{i}.phir = zeros(1,intNFrms);
        currentFrame_O1{i}.phil = zeros(1,intNFrms);
        currentFrame_O1{i}.r = zeros(1,intNFrms);
        currentFrame_O1{i}.l = zeros(1,intNFrms);
        
        currentFrame_O2{i}.xc = zeros(1,intNFrms);
        currentFrame_O2{i}.yc = zeros(1,intNFrms);
        currentFrame_O2{i}.head = zeros(1,intNFrms);
        currentFrame_O2{i}.xh = zeros(1,intNFrms);
        currentFrame_O2{i}.xt = zeros(1,intNFrms);
        currentFrame_O2{i}.yh = zeros(1,intNFrms);
        currentFrame_O2{i}.yt = zeros(1,intNFrms);
        currentFrame_O2{i}.phir = zeros(1,intNFrms);
        currentFrame_O2{i}.phil = zeros(1,intNFrms);
        currentFrame_O2{i}.r = zeros(1,intNFrms);
        currentFrame_O2{i}.l = zeros(1,intNFrms);
        
    end
    
    mearoidat = [inputFilePath inputFileName '_mearoi.mat'];
    load(mearoidat);
    
    if ispc
        mexDDGrab( 'setChambers', mea, roiCorners, params.bool_dot );
    else
        FFGrab( 'setChambers', mea, roiCorners, params.bool_dot );
    end
    
    try
        mmread(inputMovieFile, intStartFrm : intStartFrm+intNFrms, ...
            [], false, true, 'ProcessFrameCapWin', false );
    catch err
        if ~strcmp(err.message(end-14:end),'STOP PROCESSING')
            rethrow(err);
        end
    end
    
    %Save traking infomration
    trackinginfoFileName = [Files.strInVideoPath '/' Files.strInVideoFName '_trackinginfo.mat'];
    save(trackinginfoFileName, 'currentFrame_O1', 'currentFrame_O2', 'object_1', 'object_2');
    
    for i=1:params.nchambers
        fclose(Files.FeatureFID{i});
        fclose(Files.ErrorFID{i});
    end
    
end

% embedded analysis

params.analysis.nchambers = params.nchambers;
params.analysis.courtship = params.bool_court;


if params.bool_nfly
    params.analysis.oneobj = 0;
else
    params.analysis.oneobj = 1;
end

embedded_analysis(inputFilePath, inputFileName, params.analysis);

telapsed = toc;
fprintf('time elapsed: %02.0f:%02.0f mm:ss\n', floor(telapsed/60), mod(telapsed,60));
return
