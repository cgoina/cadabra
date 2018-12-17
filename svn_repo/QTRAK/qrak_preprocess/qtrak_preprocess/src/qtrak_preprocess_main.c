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

#include <stdio.h>
#include "mclmcrrt.h"
#ifdef __cplusplus
extern "C" {
#endif

extern mclComponentData __MCC_qtrak_preprocess_component_data;

#ifdef __cplusplus
}
#endif

static HMCRINSTANCE _mcr_inst = NULL;


#ifdef __cplusplus
extern "C" {
#endif

static int mclDefaultPrintHandler(const char *s)
{
  return mclWrite(1 /* stdout */, s, sizeof(char)*strlen(s));
}

#ifdef __cplusplus
} /* End extern "C" block */
#endif

#ifdef __cplusplus
extern "C" {
#endif

static int mclDefaultErrorHandler(const char *s)
{
  int written = 0;
  size_t len = 0;
  len = strlen(s);
  written = mclWrite(2 /* stderr */, s, sizeof(char)*len);
  if (len > 0 && s[ len-1 ] != '\n')
    written += mclWrite(2 /* stderr */, "\n", sizeof(char));
  return written;
}

#ifdef __cplusplus
} /* End extern "C" block */
#endif

/* This symbol is defined in shared libraries. Define it here
 * (to nothing) in case this isn't a shared library. 
 */
#ifndef LIB_qtrak_preprocess_C_API 
#define LIB_qtrak_preprocess_C_API /* No special import/export declaration */
#endif

LIB_qtrak_preprocess_C_API 
bool MW_CALL_CONV qtrak_preprocessInitializeWithHandlers(
    mclOutputHandlerFcn error_handler,
    mclOutputHandlerFcn print_handler
)
{
  if (_mcr_inst != NULL)
    return true;
  if (!mclmcrInitialize())
    return false;
  if (!mclInitializeComponentInstance(&_mcr_inst,
                                      &__MCC_qtrak_preprocess_component_data,
                                      true, NoObjectType, ExeTarget,
                                      error_handler, print_handler))
    return false;
  return true;
}

LIB_qtrak_preprocess_C_API 
bool MW_CALL_CONV qtrak_preprocessInitialize(void)
{
  return qtrak_preprocessInitializeWithHandlers(mclDefaultErrorHandler,
                                                mclDefaultPrintHandler);
}

LIB_qtrak_preprocess_C_API 
void MW_CALL_CONV qtrak_preprocessTerminate(void)
{
  if (_mcr_inst != NULL)
    mclTerminateInstance(&_mcr_inst);
}

int run_main(int argc, const char **argv)
{
  int _retval;
  /* Generate and populate the path_to_component. */
  char path_to_component[(PATH_MAX*2)+1];
  separatePathName(argv[0], path_to_component, (PATH_MAX*2)+1);
  __MCC_qtrak_preprocess_component_data.path_to_component = path_to_component; 
  if (!qtrak_preprocessInitialize()) {
    return -1;
  }
  argc = mclSetCmdLineUserData(mclGetID(_mcr_inst), argc, argv);
  _retval = mclMain(_mcr_inst, argc, argv, "qtrak_preprocess", 0);
  if (_retval == 0 /* no error */) mclWaitForFiguresToDie(NULL);
  qtrak_preprocessTerminate();
  mclTerminateApplication();
  return _retval;
}

int main(int argc, const char **argv)
{
  if (!mclInitializeApplication(
    __MCC_qtrak_preprocess_component_data.runtime_options,
    __MCC_qtrak_preprocess_component_data.runtime_option_count))
    return 0;
  
  return mclRunMain(run_main, argc, argv);
}
