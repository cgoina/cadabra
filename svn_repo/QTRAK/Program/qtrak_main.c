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

#include <stdio.h>
#include "mclmcrrt.h"
#ifdef __cplusplus
extern "C" {
#endif

extern mclComponentData __MCC_qtrak_component_data;

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
#ifndef LIB_qtrak_C_API 
#define LIB_qtrak_C_API /* No special import/export declaration */
#endif

LIB_qtrak_C_API 
bool MW_CALL_CONV qtrakInitializeWithHandlers(
    mclOutputHandlerFcn error_handler,
    mclOutputHandlerFcn print_handler
)
{
  if (_mcr_inst != NULL)
    return true;
  if (!mclmcrInitialize())
    return false;
  if (!mclInitializeComponentInstanceWithEmbeddedCTF(&_mcr_inst,
                                                     &__MCC_qtrak_component_data,
                                                     true, NoObjectType,
                                                     ExeTarget, error_handler,
                                                     print_handler, 793292, NULL))
    return false;
  return true;
}

LIB_qtrak_C_API 
bool MW_CALL_CONV qtrakInitialize(void)
{
  return qtrakInitializeWithHandlers(mclDefaultErrorHandler,
                                     mclDefaultPrintHandler);
}

LIB_qtrak_C_API 
void MW_CALL_CONV qtrakTerminate(void)
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
  __MCC_qtrak_component_data.path_to_component = path_to_component; 
  if (!qtrakInitialize()) {
    return -1;
  }
  _retval = mclMain(_mcr_inst, argc, argv, "qtrak", 0);
  if (_retval == 0 /* no error */) mclWaitForFiguresToDie(NULL);
  qtrakTerminate();
  mclTerminateApplication();
  return _retval;
}

int main(int argc, const char **argv)
{
  if (!mclInitializeApplication(
    __MCC_qtrak_component_data.runtime_options,
    __MCC_qtrak_component_data.runtime_option_count))
    return 0;
  
  return mclRunMain(run_main, argc, argv);
}
