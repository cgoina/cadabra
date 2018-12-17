mmread('./media/092208_CSMH_C1_S13N.wmv',1:30,[],false,false,'processFrame',false);

load('./media/092208_CSMH_C1_S13N_dat.mat');

try
    FFGrab( 'setChambers', mea, roiCorners, dot );
catch err
    err;
end

mmread('./media/092208_CSMH_C1_S13N.wmv',31:60,[],false,false,'processFrame',false);