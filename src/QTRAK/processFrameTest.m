function processFrameTest( data, images, width, height, frameNr, time )
    try
        
        PlotImageTest( data, images, width, height, frameNr );            

    catch err
        err.message
    end

end


function PlotImageTest( data, images, width, height, FrameNumber ) %#ok<INUSL>

    if (FrameNumber == 31),
            load('./media/092208_CSMH_C1_S13N_dat.mat');
            if ispc,
                mexDDGrab( 'setChambers', mea, roiCorners, 0 );
            else
                FFGrab( 'setChambers', mea, roiCorners, 0 );
            end
    end
  
end