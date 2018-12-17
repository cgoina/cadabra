/*
 * MATLAB Compiler: 4.9 (R2008b)
 * Date: Mon Jan 11 19:15:18 2010
 * Arguments: "-B" "macro_default" "-C" "-o" "qtrak_preprocess" "-W" "main"
 * "-d" "/home/liujin/SourceCode/compiledqtrak/qtrak_preprocess/src" "-T"
 * "link:exe" "-v" "-C"
 * "/home/liujin/SourceCode/compiledqtrak/qtrak_preprocess.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/auto_cal.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/auto_cal_circ.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/calibrate.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/callback.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/CircularHough_Grd.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/consist.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/dist2.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/ellipsefit.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/FFGrab.mexa64" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/fhist.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/fhistc.mexa64" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/fhistfit.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/fly_bodymeas.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/fly_headtailwings.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/fotsu.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/gauss.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/gmm.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/gmmactiv.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/gmmem.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/gmminit.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/gmmpost.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/gmmprob.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/GNU_message.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/histfit.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/intersect_flies.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/kmeans.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/mask2ellipses.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/mmread.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/openavi.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/PlotImage_Captioning.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/ProcessFrameCapWin.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/ProcessFrameMeanWin.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/seg_fullfly.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/seg_fullfly_inconstbox.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/seg_fullfly_inflexbox.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/segfly.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/switch_flies.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/write_sequence.m" "-a"
 * "/home/liujin/SourceCode/compiledqtrak/write_sequence_crop.m" 
 */

#include "mclmcrrt.h"

#ifdef __cplusplus
extern "C" {
#endif
const unsigned char __MCC_qtrak_preprocess_session_key[] = {
    '4', '6', 'C', 'D', '7', '4', 'F', '6', 'A', 'F', '5', '6', '9', 'B', 'F',
    'E', 'C', 'B', 'D', '3', '1', 'A', '2', '1', 'F', '7', '3', '6', '0', '5',
    'E', '6', 'A', '5', 'D', '6', 'C', '5', '7', '2', 'E', 'D', '7', '4', '9',
    '9', '7', 'A', 'E', '3', '9', '3', '8', '5', '5', 'E', 'A', '3', '9', '7',
    'C', '4', '0', '3', '0', '2', '9', '9', 'F', '8', 'B', 'E', 'F', '5', '9',
    '3', 'D', '6', '7', '8', 'B', '2', 'C', '8', 'E', '4', '8', '6', 'E', '1',
    '3', '5', '3', '1', '7', 'B', 'D', '3', 'C', '1', '1', 'A', 'B', 'D', 'C',
    '4', '4', '6', '8', '6', 'F', '3', '4', '7', '1', '8', 'D', 'A', 'B', '7',
    'B', 'C', '8', '7', 'F', 'F', '1', 'C', '4', '1', '0', 'B', 'C', '1', 'E',
    'F', '4', '1', '8', '0', 'A', '2', '1', 'A', 'C', 'F', '7', 'F', '2', '8',
    'F', '6', '8', '5', '8', '4', 'E', '0', 'D', 'F', '9', '1', '3', '8', 'E',
    'F', '4', '4', '0', '3', '6', '0', '1', '6', 'C', '1', '8', '7', 'B', 'F',
    'C', '3', '6', 'E', '5', '4', '7', '2', 'C', '0', '5', '2', '2', '0', 'E',
    '2', 'C', '9', '2', 'D', 'F', 'E', '3', '3', '8', '7', 'A', '0', '8', '9',
    '2', '2', '8', '9', '9', '0', '2', 'B', 'F', '5', '8', 'D', 'B', '8', '8',
    '7', '5', '6', 'B', 'C', '7', 'C', 'E', 'F', '1', '3', '0', '2', '0', 'E',
    '2', 'B', '1', '4', 'A', 'D', '2', 'E', '8', '1', 'C', 'C', 'B', '8', '0',
    '2', '\0'};

const unsigned char __MCC_qtrak_preprocess_public_key[] = {
    '3', '0', '8', '1', '9', 'D', '3', '0', '0', 'D', '0', '6', '0', '9', '2',
    'A', '8', '6', '4', '8', '8', '6', 'F', '7', '0', 'D', '0', '1', '0', '1',
    '0', '1', '0', '5', '0', '0', '0', '3', '8', '1', '8', 'B', '0', '0', '3',
    '0', '8', '1', '8', '7', '0', '2', '8', '1', '8', '1', '0', '0', 'C', '4',
    '9', 'C', 'A', 'C', '3', '4', 'E', 'D', '1', '3', 'A', '5', '2', '0', '6',
    '5', '8', 'F', '6', 'F', '8', 'E', '0', '1', '3', '8', 'C', '4', '3', '1',
    '5', 'B', '4', '3', '1', '5', '2', '7', '7', 'E', 'D', '3', 'F', '7', 'D',
    'A', 'E', '5', '3', '0', '9', '9', 'D', 'B', '0', '8', 'E', 'E', '5', '8',
    '9', 'F', '8', '0', '4', 'D', '4', 'B', '9', '8', '1', '3', '2', '6', 'A',
    '5', '2', 'C', 'C', 'E', '4', '3', '8', '2', 'E', '9', 'F', '2', 'B', '4',
    'D', '0', '8', '5', 'E', 'B', '9', '5', '0', 'C', '7', 'A', 'B', '1', '2',
    'E', 'D', 'E', '2', 'D', '4', '1', '2', '9', '7', '8', '2', '0', 'E', '6',
    '3', '7', '7', 'A', '5', 'F', 'E', 'B', '5', '6', '8', '9', 'D', '4', 'E',
    '6', '0', '3', '2', 'F', '6', '0', 'C', '4', '3', '0', '7', '4', 'A', '0',
    '4', 'C', '2', '6', 'A', 'B', '7', '2', 'F', '5', '4', 'B', '5', '1', 'B',
    'B', '4', '6', '0', '5', '7', '8', '7', '8', '5', 'B', '1', '9', '9', '0',
    '1', '4', '3', '1', '4', 'A', '6', '5', 'F', '0', '9', '0', 'B', '6', '1',
    'F', 'C', '2', '0', '1', '6', '9', '4', '5', '3', 'B', '5', '8', 'F', 'C',
    '8', 'B', 'A', '4', '3', 'E', '6', '7', '7', '6', 'E', 'B', '7', 'E', 'C',
    'D', '3', '1', '7', '8', 'B', '5', '6', 'A', 'B', '0', 'F', 'A', '0', '6',
    'D', 'D', '6', '4', '9', '6', '7', 'C', 'B', '1', '4', '9', 'E', '5', '0',
    '2', '0', '1', '1', '1', '\0'};

static const char * MCC_qtrak_preprocess_matlabpath_data[] = 
  { "qtrak_prepro/", "$TOOLBOXDEPLOYDIR/",
    "$TOOLBOXMATLABDIR/general/", "$TOOLBOXMATLABDIR/ops/",
    "$TOOLBOXMATLABDIR/lang/", "$TOOLBOXMATLABDIR/elmat/",
    "$TOOLBOXMATLABDIR/randfun/", "$TOOLBOXMATLABDIR/elfun/",
    "$TOOLBOXMATLABDIR/specfun/", "$TOOLBOXMATLABDIR/matfun/",
    "$TOOLBOXMATLABDIR/datafun/", "$TOOLBOXMATLABDIR/polyfun/",
    "$TOOLBOXMATLABDIR/funfun/", "$TOOLBOXMATLABDIR/sparfun/",
    "$TOOLBOXMATLABDIR/scribe/", "$TOOLBOXMATLABDIR/graph2d/",
    "$TOOLBOXMATLABDIR/graph3d/", "$TOOLBOXMATLABDIR/specgraph/",
    "$TOOLBOXMATLABDIR/graphics/", "$TOOLBOXMATLABDIR/uitools/",
    "$TOOLBOXMATLABDIR/strfun/", "$TOOLBOXMATLABDIR/imagesci/",
    "$TOOLBOXMATLABDIR/iofun/", "$TOOLBOXMATLABDIR/audiovideo/",
    "$TOOLBOXMATLABDIR/timefun/", "$TOOLBOXMATLABDIR/datatypes/",
    "$TOOLBOXMATLABDIR/verctrl/", "$TOOLBOXMATLABDIR/codetools/",
    "$TOOLBOXMATLABDIR/helptools/", "$TOOLBOXMATLABDIR/demos/",
    "$TOOLBOXMATLABDIR/timeseries/", "$TOOLBOXMATLABDIR/hds/",
    "$TOOLBOXMATLABDIR/guide/", "$TOOLBOXMATLABDIR/plottools/",
    "toolbox/local/", "toolbox/shared/dastudio/",
    "$TOOLBOXMATLABDIR/datamanager/", "toolbox/compiler/",
    "toolbox/images/colorspaces/", "toolbox/images/images/",
    "toolbox/images/imuitools/", "toolbox/images/iptformats/",
    "toolbox/images/iptutils/", "toolbox/shared/imageslib/",
    "toolbox/shared/spcuilib/", "toolbox/signal/signal/" };

static const char * MCC_qtrak_preprocess_classpath_data[] = 
  { "java/jar/toolbox/images.jar" };

static const char * MCC_qtrak_preprocess_libpath_data[] = 
  { "" };

static const char * MCC_qtrak_preprocess_app_opts_data[] = 
  { "" };

static const char * MCC_qtrak_preprocess_run_opts_data[] = 
  { "" };

static const char * MCC_qtrak_preprocess_warning_state_data[] = 
  { "off:MATLAB:dispatcher:nameConflict" };


mclComponentData __MCC_qtrak_preprocess_component_data = { 

  /* Public key data */
  __MCC_qtrak_preprocess_public_key,

  /* Component name */
  "qtrak_preprocess",

  /* Component Root */
  "",

  /* Application key data */
  __MCC_qtrak_preprocess_session_key,

  /* Component's MATLAB Path */
  MCC_qtrak_preprocess_matlabpath_data,

  /* Number of directories in the MATLAB Path */
  46,

  /* Component's Java class path */
  MCC_qtrak_preprocess_classpath_data,
  /* Number of directories in the Java class path */
  1,

  /* Component's load library path (for extra shared libraries) */
  MCC_qtrak_preprocess_libpath_data,
  /* Number of directories in the load library path */
  0,

  /* MCR instance-specific runtime options */
  MCC_qtrak_preprocess_app_opts_data,
  /* Number of MCR instance-specific runtime options */
  0,

  /* MCR global runtime options */
  MCC_qtrak_preprocess_run_opts_data,
  /* Number of MCR global runtime options */
  0,
  
  /* Component preferences directory */
  "qtrak_prepro_8274B94230CA79FC2FF2C7CA08AB2E32",

  /* MCR warning status data */
  MCC_qtrak_preprocess_warning_state_data,
  /* Number of MCR warning status modifiers */
  1,

  /* Path to component - evaluated at runtime */
  NULL

};

#ifdef __cplusplus
}
#endif


