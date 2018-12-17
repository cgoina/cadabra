/*
 * MATLAB Compiler: 4.8 (R2008a)
 * Date: Thu May 28 14:37:51 2009
 * Arguments: "-B" "macro_default" "-o" "qtrak" "-W" "main" "-d" "C:\Documents
 * and Settings\liuj\My Documents\Matlab\Qtrak\qtrak\src" "-T" "link:exe" "-v"
 * "C:\Documents and Settings\liuj\My Documents\Matlab\Qtrak\qtrak.m" "-a"
 * "C:\Documents and Settings\liuj\My Documents\Matlab\Qtrak\auto_cal.m" "-a"
 * "C:\Documents and Settings\liuj\My Documents\Matlab\Qtrak\auto_cal_circ.m"
 * "-a" "C:\Documents and Settings\liuj\My Documents\Matlab\Qtrak\calibrate.m"
 * "-a" "C:\Documents and Settings\liuj\My Documents\Matlab\Qtrak\callback.m"
 * "-a" "C:\Documents and Settings\liuj\My
 * Documents\Matlab\Qtrak\CircularHough_Grd.m" "-a" "C:\Documents and
 * Settings\liuj\My Documents\Matlab\Qtrak\consist.m" "-a" "C:\Documents and
 * Settings\liuj\My Documents\Matlab\Qtrak\dist2.m" "-a" "C:\Documents and
 * Settings\liuj\My Documents\Matlab\Qtrak\ellipsefit.m" "-a" "C:\Documents and
 * Settings\liuj\My Documents\Matlab\Qtrak\fhist.m" "-a" "C:\Documents and
 * Settings\liuj\My Documents\Matlab\Qtrak\fhistc.mexw32" "-a" "C:\Documents
 * and Settings\liuj\My Documents\Matlab\Qtrak\fhistfit.m" "-a" "C:\Documents
 * and Settings\liuj\My Documents\Matlab\Qtrak\fly_bodymeas.m" "-a"
 * "C:\Documents and Settings\liuj\My
 * Documents\Matlab\Qtrak\fly_headtailwings.m" "-a" "C:\Documents and
 * Settings\liuj\My Documents\Matlab\Qtrak\fotsu.m" "-a" "C:\Documents and
 * Settings\liuj\My Documents\Matlab\Qtrak\gauss.m" "-a" "C:\Documents and
 * Settings\liuj\My Documents\Matlab\Qtrak\gmm.m" "-a" "C:\Documents and
 * Settings\liuj\My Documents\Matlab\Qtrak\gmmactiv.m" "-a" "C:\Documents and
 * Settings\liuj\My Documents\Matlab\Qtrak\gmmem.m" "-a" "C:\Documents and
 * Settings\liuj\My Documents\Matlab\Qtrak\gmminit.m" "-a" "C:\Documents and
 * Settings\liuj\My Documents\Matlab\Qtrak\gmmpost.m" "-a" "C:\Documents and
 * Settings\liuj\My Documents\Matlab\Qtrak\gmmprob.m" "-a" "C:\Documents and
 * Settings\liuj\My Documents\Matlab\Qtrak\GNU_message.m" "-a" "C:\Documents
 * and Settings\liuj\My Documents\Matlab\Qtrak\histfit.m" "-a" "C:\Documents
 * and Settings\liuj\My Documents\Matlab\Qtrak\intersect_flies.m" "-a"
 * "C:\Documents and Settings\liuj\My Documents\Matlab\Qtrak\kmeans.m" "-a"
 * "C:\Documents and Settings\liuj\My Documents\Matlab\Qtrak\mask2ellipses.m"
 * "-a" "C:\Documents and Settings\liuj\My
 * Documents\Matlab\Qtrak\mexDDGrab.mexw32" "-a" "C:\Documents and
 * Settings\liuj\My Documents\Matlab\Qtrak\mmread.m" "-a" "C:\Documents and
 * Settings\liuj\My Documents\Matlab\Qtrak\openavi.m" "-a" "C:\Documents and
 * Settings\liuj\My Documents\Matlab\Qtrak\PlotImage_Captioning.m" "-a"
 * "C:\Documents and Settings\liuj\My Documents\Matlab\Qtrak\processFrame.m"
 * "-a" "C:\Documents and Settings\liuj\My
 * Documents\Matlab\Qtrak\ProcessFrameCapWin.m" "-a" "C:\Documents and
 * Settings\liuj\My Documents\Matlab\Qtrak\ProcessFrameMeanWin.m" "-a"
 * "C:\Documents and Settings\liuj\My
 * Documents\Matlab\Qtrak\processFrameTest.m" "-a" "C:\Documents and
 * Settings\liuj\My Documents\Matlab\Qtrak\readlog.m" "-a" "C:\Documents and
 * Settings\liuj\My Documents\Matlab\Qtrak\seg_fullfly.m" "-a" "C:\Documents
 * and Settings\liuj\My Documents\Matlab\Qtrak\seg_fullfly_inconstbox.m" "-a"
 * "C:\Documents and Settings\liuj\My
 * Documents\Matlab\Qtrak\seg_fullfly_inflexbox.m" "-a" "C:\Documents and
 * Settings\liuj\My Documents\Matlab\Qtrak\segfly.m" "-a" "C:\Documents and
 * Settings\liuj\My Documents\Matlab\Qtrak\switch_flies.m" "-a" "C:\Documents
 * and Settings\liuj\My Documents\Matlab\Qtrak\test.m" "-a" "C:\Documents and
 * Settings\liuj\My Documents\Matlab\Qtrak\write_sequence.m" "-a" "C:\Documents
 * and Settings\liuj\My Documents\Matlab\Qtrak\write_sequence_crop.m" 
 */

#include "mclmcrrt.h"

#ifdef __cplusplus
extern "C" {
#endif
const unsigned char __MCC_qtrak_session_key[] = {
    '3', 'E', 'A', '6', '9', 'F', '7', '1', '0', '7', '7', 'C', 'B', '2', 'A',
    '9', '4', '8', 'E', 'E', 'B', '5', 'C', 'B', 'F', 'B', '1', '1', '7', '0',
    'D', 'C', 'A', '8', 'A', '3', 'F', 'F', 'F', '5', '7', '3', '2', 'E', '2',
    '7', '8', '0', '2', '4', '3', 'A', '1', '3', '6', '1', 'D', 'D', '9', '7',
    'A', '6', '3', '3', '4', 'E', 'A', '0', '6', '9', '5', '7', '0', '7', 'B',
    '5', '0', 'D', 'A', 'D', 'E', '0', '5', 'B', 'C', 'A', 'A', 'A', '0', 'B',
    '1', '9', '1', '2', '1', 'E', 'F', '4', '4', '9', 'B', 'A', 'E', 'D', '7',
    '7', '0', '0', 'B', 'A', '1', 'C', '6', '9', 'C', '7', 'F', '1', '6', 'E',
    '2', '2', '4', '1', 'E', 'D', '8', '4', '6', 'C', 'C', 'B', 'F', '6', 'C',
    '6', 'C', 'F', '4', '9', '1', '7', '9', 'E', '9', 'A', 'C', '5', '1', '2',
    '4', '6', '2', '5', '7', '1', '0', 'F', '5', 'F', 'A', 'C', '1', '0', '9',
    'A', '5', '6', 'D', 'B', '6', '4', '6', '9', '7', '9', '8', '8', '6', '0',
    'B', '2', 'D', 'D', '6', 'D', '7', '7', '6', '5', '9', '2', '3', '5', '3',
    '3', '1', '2', 'A', '4', 'B', 'B', 'A', '6', 'A', '5', '4', '1', '3', 'C',
    'D', '0', '6', 'E', 'A', '5', '7', '2', 'A', '9', '3', '8', '7', 'B', '2',
    '7', '2', 'F', '1', '8', 'B', '2', '7', '6', '6', '6', '6', '7', 'E', 'E',
    'F', '3', '9', '4', '1', 'C', '1', '1', '5', '6', '8', '8', 'D', 'B', '8',
    'E', '\0'};

const unsigned char __MCC_qtrak_public_key[] = {
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

static const char * MCC_qtrak_matlabpath_data[] = 
  { "qtrak/", "toolbox/compiler/deploy/", "$TOOLBOXMATLABDIR/general/",
    "$TOOLBOXMATLABDIR/ops/", "$TOOLBOXMATLABDIR/lang/",
    "$TOOLBOXMATLABDIR/elmat/", "$TOOLBOXMATLABDIR/elfun/",
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
    "$TOOLBOXMATLABDIR/helptools/", "$TOOLBOXMATLABDIR/winfun/",
    "$TOOLBOXMATLABDIR/demos/", "$TOOLBOXMATLABDIR/timeseries/",
    "$TOOLBOXMATLABDIR/hds/", "$TOOLBOXMATLABDIR/guide/",
    "$TOOLBOXMATLABDIR/plottools/", "toolbox/local/",
    "toolbox/shared/dastudio/", "$TOOLBOXMATLABDIR/datamanager/",
    "toolbox/compiler/", "toolbox/images/images/",
    "toolbox/images/imuitools/", "toolbox/images/iptutils/",
    "toolbox/shared/imageslib/", "toolbox/images/medformats/",
    "toolbox/shared/spcuilib/", "toolbox/signal/signal/" };

static const char * MCC_qtrak_classpath_data[] = 
  { "java/jar/toolbox/images.jar" };

static const char * MCC_qtrak_libpath_data[] = 
  { "" };

static const char * MCC_qtrak_app_opts_data[] = 
  { "" };

static const char * MCC_qtrak_run_opts_data[] = 
  { "" };

static const char * MCC_qtrak_warning_state_data[] = 
  { "off:MATLAB:dispatcher:nameConflict" };


mclComponentData __MCC_qtrak_component_data = { 

  /* Public key data */
  __MCC_qtrak_public_key,

  /* Component name */
  "qtrak",

  /* Component Root */
  "",

  /* Application key data */
  __MCC_qtrak_session_key,

  /* Component's MATLAB Path */
  MCC_qtrak_matlabpath_data,

  /* Number of directories in the MATLAB Path */
  45,

  /* Component's Java class path */
  MCC_qtrak_classpath_data,
  /* Number of directories in the Java class path */
  1,

  /* Component's load library path (for extra shared libraries) */
  MCC_qtrak_libpath_data,
  /* Number of directories in the load library path */
  0,

  /* MCR instance-specific runtime options */
  MCC_qtrak_app_opts_data,
  /* Number of MCR instance-specific runtime options */
  0,

  /* MCR global runtime options */
  MCC_qtrak_run_opts_data,
  /* Number of MCR global runtime options */
  0,
  
  /* Component preferences directory */
  "qtrak_BBA0934415E7B46C795195F36C8C7F90",

  /* MCR warning status data */
  MCC_qtrak_warning_state_data,
  /* Number of MCR warning status modifiers */
  1,

  /* Path to component - evaluated at runtime */
  NULL

};

#ifdef __cplusplus
}
#endif


