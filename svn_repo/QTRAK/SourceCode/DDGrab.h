/***************************************************
Copyright (C) Micah Richert

This file is part of mmread.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Written by Micah Richert.
07/14/2005
Changed by Edwin Soedarmadji
10/15/2008
**************************************************/

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif

#if defined(_DEBUG) && defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK,__FILE__,__LINE__)
#endif

#include "dshow.h"
#include "qedit.h"
#include <assert.h>

template<class T> class vector
{
    public:
        vector()
        {
            datavec = NULL;
            datavecSize = 0;
            nr = 0;
            datavecSize = resize(128);
        }

        ~vector()
        {
            free(datavec);
        }

        T at(unsigned int i)
        {
            if (i >= nr) return NULL;
            return datavec[i];
        }

        inline int size()
        {
            return nr;
        }

        void assign(unsigned int i, T data)
        {
            if (i >= nr) return;
            datavec[i] = data;
        }

        void add(T data)
        {
            nr++;
            if (nr > datavecSize)
            {
                datavecSize = resize(datavecSize*2);
                if (nr > datavecSize)
                {
                    //we've failed, but can't return an error here...
                    nr--;
                    return;
                }
            }
            datavec[nr-1] = data;
        }

        void clear()
        {
            nr = 0;
            datavecSize = resize(128);
        }

    private:
        T* datavec;
        unsigned int nr;
        unsigned int datavecSize;

        int resize(int newsize)
        {
            void* tmp = realloc(datavec,newsize*sizeof(T));
            if (tmp)
            {
                datavec = (T*)tmp;
                return newsize;
            }
            return datavecSize;
        }
};

// since the Audio and Video CB vectors are public we need to make the CB interface public too
class CSampleGrabberCB : public ISampleGrabberCB
{
public:
    CSampleGrabberCB();
    virtual ~CSampleGrabberCB();

    vector<BYTE*> frames;
    vector<int> frameBytes;
    vector<int> frameNrs;
    vector<double> frameTimes;

    // use this to get data format information, ie. bit depth, sampling rate...
    BYTE *pbFormat;
    GUID subtype; // what is the encoding

    unsigned int frameNr;
    bool disabled;
    bool done;
    bool isAudio;

    int bytesPerWORD;
    double rate;
    double startTime, stopTime;

    // Fake out any COM ref counting
    //
    STDMETHODIMP_(ULONG) AddRef() { return 2; }
    STDMETHODIMP_(ULONG) Release() { return 1; }

    // Fake out any COM QI'ing
    //
    STDMETHODIMP QueryInterface(REFIID riid, void ** ppv);

    // We don't implement this one
    //
    STDMETHODIMP SampleCB( double SampleTime, IMediaSample * pSample ){ return 0; }

    // The sample grabber is calling us back on its deliver thread.
    // This is NOT the main app thread!
    //
    STDMETHODIMP BufferCB( double dblSampleTime, BYTE * pBuffer, long lBufferSize );
};

/* The bulk of DDGrab.h remains the same, except for the portion enclosed by the        */
/* conditional compilation block marked by the preprocessor directive MATLAB_MEX_FILE.  */  
/* For completeness, however, the entire class definition DDGrabber is included below.  */
/* The modified portion is listed on the next section.                                  */

// this is the main grabber class.  I think the interfaces and names are fairly self explanatory
class DDGrabber
{
public:
    vector<CSampleGrabberCB*> VideoCBs;
    vector<CSampleGrabberCB*> AudioCBs;
    DDGrabber();

    HRESULT buildGraph(char *filename);
    HRESULT doCapture();
    HRESULT getVideoInfo(unsigned int id, int *width, int *height, double *rate, int *nrFramesCaptured, int *nrFramesTotal, double *totalDuration);
    HRESULT getAudioInfo(unsigned int id, int *nrChannels, double *rate, int* bits, int *nrFramesCaptured, int *nrFramesTotal, GUID *subtype, double *totalDuration);
    void getCaptureInfo(int *nrVideo, int *nrAudio);
    // data must be freed by caller
    HRESULT getVideoFrame(unsigned int id, int frameNr, BYTE **data, int *nrBytes, double *time);
    // data must be freed by caller
    HRESULT getAudioFrame(unsigned int id, int frameNr, BYTE **data, int *nrBytes, double *time);
    void setFrames(int *frameNrs, int nrFrames);
    void setTime(double startTime, double stopTime);
    void setTrySeeking(int tryseek);
    void disableVideo();
    void disableAudio();
    void cleanUp(); // must be called at the end, in order to render anything afterward.

/* The first  modification is the introduction of the function setChambers. Because     */
/* DDGrabber now processes the pixels instead of just delivering video frames to the    */
/* Matlab routine,  it needs to know the boundary of each fly chamber in the image.     */
/* This is provided by the Matlab routine through setChambers.  For each chamber, a     */ 
/* structure is defined to store the chamber’s parameters.  This structure stores a     */
/* parameter indicating whether one of the flies in each chamber is tagged with a dot.  */

#ifdef MATLAB_MEX_FILE
    HRESULT SaveGraphFile(IGraphBuilder *pGraph, WCHAR *wszPath);	
    void setChambers(const mxArray *chambers, const mxArray *corners, bool dot);
    void setMatlabCommand(char *matlabCommand);
    void runMatlabCommand();
#endif
private:
    IGraphBuilder *pGraphBuilder;
    bool stopForced;	
    bool tryseeking;
    vector<int> frameNrs;
    double startTime, stopTime;

#ifdef MATLAB_MEX_FILE
#define CS 0
    /* Except for lastFrames and matlabCommand, other variables are     	*/
    /* introduced specifically for QTracker and the setChambers method  	*/

    vector<int> lastFrames;		            /* The frame buffer 		    */
    char *matlabCommand;			        /* Matlab callback routine	    */

    bool isDot; 				            /* Does one fly have a dot?  	*/
    int noofChambers;			            /* The number of fly chambers	*/
    struct t_chamber			            /* Structure for each chamber	*/
    {   
        int R0, R1, C0, C1;		            /* ROI corner points 		    */
        int roiDims[2];			            /* ROI matrix dimension	        */
        double *pimg, *pimg2, *pmean;	    /* Pointer to pixel values	    */
        mxArray *mxaimg, *mxaimg2;		    /* Pointer to Matlab arrays	    */
    } 
    *pChamber;				                /* A pointer instance		    */
    
#endif

    void MyFreeMediaType(AM_MEDIA_TYPE& mt);
    PIN_INFO getPinInfo(IPin* pin);
    IPin* getInputPin(IBaseFilter* filt);
    IPin* getOutputPin(IBaseFilter* filt);
    bool isRenderer(IBaseFilter* filt);
    IPin* connectedToInput(IBaseFilter* filt);
    GUID getMajorType(IBaseFilter* filt);
    HRESULT insertCapture(IGraphBuilder* pGraphBuilder, IBaseFilter* pRenderer, AM_MEDIA_TYPE* mt, CSampleGrabberCB** grabberCB);
    HRESULT insertVideoCapture(IGraphBuilder* pGraphBuilder, IBaseFilter* pRenderer);
    HRESULT insertAudioCapture(IGraphBuilder* pGraphBuilder, IBaseFilter* pRenderer);
    HRESULT changeToNull(IGraphBuilder* pGraphBuilder, IBaseFilter* pRenderer);
    HRESULT mangleGraph(IGraphBuilder* pGraphBuilder);
};