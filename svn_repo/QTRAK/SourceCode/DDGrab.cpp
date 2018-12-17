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

This is the main Grabber code.  It builds the Direct Show graphs
and then inserts SampleGrabbers right before the renderers.  It
then changes all of the renderers to NullRenderers.  This code
supports any number of audio or video streams.  Raw midi streams
are not supported -- I didn't think I should bother.

This code was intended to be used inside of a matlab interface,
but can be used as a generic grabber class for anyone who needs
one.

07/14/2005  Written by Micah Richert.
09/21/2005  Fixed a bug that when no frames are specified 
            it wouldn't capture all frames
10/15/2008  Modified by Edwin Soedarmadji specifically for 
            Caltech QTracker project
**************************************************/

#include "windows.h"
#include "atlbase.h"
#include "DDGrab.h"

// To detect memory leaks, we have to put code in each .cpp file
#if defined(_DEBUG) && defined(_MSC_VER)
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#define new new(_NORMAL_BLOCK,__FILE__,__LINE__)
#endif

#ifdef _DEBUG
#include <sys/timeb.h>
#include <time.h>
#endif

#include <math.h>
#include <assert.h>

GUID NULLGUID = {0};

CSampleGrabberCB::CSampleGrabberCB()
{
	pbFormat = NULL;
	frameNr = 0;
	disabled = false;
	done = false;
	bytesPerWORD = 0;
	rate = 0;
	startTime = 0;
	stopTime = 0;
	isAudio = false;
}

CSampleGrabberCB::~CSampleGrabberCB()
{
	if (pbFormat) free(pbFormat);

	_RPT1(_CRT_WARN,"Freeing %d frames\n",frames.size());

	for (int f=0;f<frames.size();f++) if (frames.at(f)) 
	{
		delete frames.at(f);
	}
	frames.clear();
	frameBytes.clear();
	frameNrs.clear();
	frameTimes.clear();

	Release();
}

// Fake out any COM QI'ing
//
STDMETHODIMP CSampleGrabberCB::QueryInterface(REFIID riid, void ** ppv)
{
	if (!ppv) return E_POINTER;

	if( riid == IID_ISampleGrabberCB || riid == IID_IUnknown )
	{
		*ppv = (void *) static_cast<ISampleGrabberCB*> ( this );
		return NOERROR;
	}

	return E_NOINTERFACE;
}

// The sample grabber is calling us back on its deliver thread.
// This is NOT the main app thread!

STDMETHODIMP CSampleGrabberCB::BufferCB( double dblSampleTime, BYTE * pBuffer, long lBufferSize )
{
#ifdef _DEBUG
   struct _timeb timebuffer;
   time_t time1;
   _ftime64_s( &timebuffer );
#endif

	if (disabled) return 0;

	if (!pBuffer) return E_POINTER;

	frameNr++;

	int nrWORDs = lBufferSize/bytesPerWORD;

	int offset=0, len=0;

    _RPT2(_CRT_WARN,"dblSampleTime: %f  frameNr: %d  realTime: %d:%d ", \
        (float)dblSampleTime,frameNr, timebuffer.time, timebuffer.millitm);

    if (frameNrs.size() > 0)
	{
		// frames are being specified, so we must be a videoCB...
		// check to see if the frame is in our list
		bool foundNr = false;
		int lastFrameNr = 0;
		for (int i=0;i<frameNrs.size();i++)
		{
			if (frameNrs.at(i) == frameNr) foundNr = true;
			if (frameNrs.at(i) > lastFrameNr) lastFrameNr = frameNrs.at(i);
		}

		_RPT1(_CRT_WARN,"lastFrameNr: %d\n",lastFrameNr);

		if (frameNr > lastFrameNr)
		{
			done = true;
			return 0;
		}

		if (foundNr) len = lBufferSize; // since this has to be a video frame, we just capture everything
	} else {
		// either no frames are specified (capture all), or we have time specified
		if (stopTime)
		{
			if (isAudio)
			{
				// time is being used...

				if (dblSampleTime>=startTime)
				{
					// if we've reached the start...
					offset = min(max(0,((int)((startTime-dblSampleTime)*rate))*bytesPerWORD),lBufferSize+1);
					len = min(lBufferSize,((int)((stopTime-dblSampleTime)*rate))*bytesPerWORD);
					// if we have gone past our stop time...

					_RPT2(_CRT_WARN," startTime-dblSampleTime: %d stopTime-dblSampleTime: %d\n",((int)((startTime-dblSampleTime)*rate))*bytesPerWORD,((int)((stopTime-dblSampleTime)*rate))*bytesPerWORD);
					
					_RPT4(_CRT_WARN,"startTime: %f stopTime: %f dblSampleTime: %f offset: %d",startTime,stopTime,dblSampleTime,offset);
					_RPT1(_CRT_WARN," len: %d \n",len);

					if (len < 0)
					{
						_RPT0(_CRT_WARN,"We should stop now!!! \n");
						done = true;
						return 0;
					}
				}
			} else {
				done = stopTime <= dblSampleTime;
				len = (startTime <= dblSampleTime)?lBufferSize:0;
 				
				if (done) return 0;
			}
		} else {
			// capture everything... video or audio
			len = lBufferSize;
			_RPT0(_CRT_WARN,"capture everything\n");
		}
	}

	if (len)
	{
		BYTE* tmp = (BYTE*)malloc(len);
		if (!tmp) return E_OUTOFMEMORY;

		memcpy(tmp,pBuffer+offset,len);

        frames.add(tmp);
		frameBytes.add(len);
		frameTimes.add(dblSampleTime);
	}

	return 0;
}

DDGrabber::DDGrabber()
{
	stopForced = false;	
	tryseeking = true;
#ifdef MATLAB_MEX_FILE
	isDot = false;
	noofChambers = 0;
#endif
}

void DDGrabber::cleanUp()
{
	IMediaControl* pMediaControl;
	IEnumFilters* filterList;
	IBaseFilter* filt;
	int i;

	for (i=0; i<VideoCBs.size(); i++) delete VideoCBs.at(i);
	for (i=0; i<AudioCBs.size(); i++) delete AudioCBs.at(i);

	VideoCBs.clear();
	AudioCBs.clear();

	// make sure the graph is stopped
	pGraphBuilder->QueryInterface(IID_IMediaControl, (void **)&pMediaControl);
	pMediaControl->Stop();
	pMediaControl->Release();

    _RPT0(_CRT_WARN,"Releasing... \n");
	if (SUCCEEDED(pGraphBuilder->EnumFilters(&filterList)))
	{
		while (filterList->Next(1, &filt, NULL) == S_OK)
		{
			#ifdef _DEBUG
				FILTER_INFO info;
				filt->QueryFilterInfo(&info);
				char str[100];
				WideCharToMultiByte( CP_ACP, 0, info.achName, -1, str, 100, NULL, NULL );
				_RPT1(_CRT_WARN,"Releasing: %s\n",str);
			#endif
			pGraphBuilder->RemoveFilter(filt);
			filterList->Reset();
			filt->Release();
		}
		filterList->Release();
	}
	pGraphBuilder->Release();
	pGraphBuilder = NULL;

#ifdef MATLAB_MEX_FILE
	if (matlabCommand) free(matlabCommand);
	matlabCommand = NULL;
#endif

   _CrtDumpMemoryLeaks();
}

void DDGrabber::MyFreeMediaType(AM_MEDIA_TYPE& mt)
{
	if (mt.cbFormat != 0)
	{
		CoTaskMemFree((PVOID)mt.pbFormat);
		mt.cbFormat = 0;
		mt.pbFormat = NULL;
	}
	if (mt.pUnk != NULL)
	{
		// Unnecessary because pUnk should not be used, but safest.
		mt.pUnk->Release();
		mt.pUnk = NULL;
	}
}

PIN_INFO DDGrabber::getPinInfo(IPin* pin)
{
	PIN_INFO	info = {0};

	if (pin)
	{
		if (SUCCEEDED(pin->QueryPinInfo(&info)))
		{
			info.pFilter->Release();
		}
	}

	return info;
}

IPin* DDGrabber::getInputPin(IBaseFilter* filt)
{
	IPin* pin = NULL;
	IEnumPins* pinList;
	ULONG tmp;

	if (!filt) return NULL;

	//get the input
	if (SUCCEEDED(filt->EnumPins(&pinList)))
	{
		pinList->Reset();
		while (pinList->Next(1, &pin, &tmp) == S_OK && getPinInfo(pin).dir != PINDIR_INPUT);
		pinList->Release();

		if (getPinInfo(pin).dir != PINDIR_INPUT) return NULL;
	}

	return pin;
}

IPin* DDGrabber::getOutputPin(IBaseFilter* filt)
{
	IPin* pin = NULL;
	IEnumPins* pinList;
	ULONG tmp;

	if (!filt) return NULL;

	//get the output
	if (SUCCEEDED(filt->EnumPins(&pinList)))
	{
		pinList->Reset();
		while (pinList->Next(1, &pin, &tmp) == S_OK && getPinInfo(pin).dir != PINDIR_OUTPUT);
		pinList->Release();

		if (getPinInfo(pin).dir != PINDIR_OUTPUT) return NULL;
	}

	return pin;
}

bool DDGrabber::isRenderer(IBaseFilter* filt)
{
	if (!filt) return false;

	IEnumPins*	pinList;
	int nrOutput = 0;
	int nrInput = 0;
	IPin*		pin = NULL;
	ULONG		tmp;

	if (FAILED(filt->EnumPins(&pinList))) return false;
	pinList->Reset();
	while (pinList->Next(1, &pin, &tmp) == S_OK)
	{
		if (getPinInfo(pin).dir == PINDIR_OUTPUT) nrOutput++;
		else nrInput++;
		pin->Release();
	}
	pinList->Release();

	#ifdef _DEBUG
		FILTER_INFO info;
		filt->QueryFilterInfo(&info);
		char str[100];
		WideCharToMultiByte( CP_ACP, 0, info.achName, -1, str, 100, NULL, NULL );
		_RPT0(_CRT_WARN,str);
		_RPT2(_CRT_WARN," %d %d\n", nrOutput, nrInput);
	#endif

	return nrOutput == 0 && nrInput == 1;  // the only filters that have no outputs are renderers
}

IPin* DDGrabber::connectedToInput(IBaseFilter* filt)
{
	if (!filt) return NULL;

	IPin* inPin;
	IPin* outPin;

	inPin = getInputPin(filt);
	if (!inPin) return NULL;

	if (FAILED(inPin->ConnectedTo(&outPin))) return NULL;
	return outPin;
}

GUID DDGrabber::getMajorType(IBaseFilter* filt)
{
	if (!filt) return NULLGUID;

	IPin* inPin;
	inPin = getInputPin(filt);
	if (!inPin) return NULLGUID;

	AM_MEDIA_TYPE mt = {0};
	if (FAILED(inPin->ConnectionMediaType(&mt))) return NULLGUID;

	GUID ret = mt.majortype;
	MyFreeMediaType(mt);

	return ret;
}

HRESULT DDGrabber::getVideoInfo(unsigned int id, int* width, int* height, double* rate, int* nrFramesCaptured, int* nrFramesTotal, double* totalDuration)
{
	if (!width || !height || !nrFramesCaptured || !nrFramesTotal) return E_POINTER;

	if (id >= VideoCBs.size()) return E_NOINTERFACE;
	CSampleGrabberCB* CB = VideoCBs.at(id);

	if (!CB) return E_POINTER;

	if (!CB->pbFormat)
	{
		*width = 0;
		*height = 0;
		*rate = 0;
	} else {
		VIDEOINFOHEADER * h = (VIDEOINFOHEADER*) CB->pbFormat;
		if (!h) return E_POINTER;
		*width  = h->bmiHeader.biWidth;
		*height = h->bmiHeader.biHeight;
		if (h->AvgTimePerFrame == 0) *rate = 0; // can't calc rate correctly...
		else *rate = 10000000.0/h->AvgTimePerFrame; // make it samples per second.
	}
	*nrFramesCaptured = CB->frames.size();
	*nrFramesTotal = CB->frameNr;

    *totalDuration = 0;
    
	// see if we can get more reliable nrFramesTotal and rate information
	IMediaSeeking* pMediaSeeking;
	pGraphBuilder->QueryInterface(IID_IMediaSeeking, (void **)&pMediaSeeking);
	if (pMediaSeeking)
	{
		LONGLONG duration = 0;
		LONGLONG durationus = 0;

		if (SUCCEEDED(pMediaSeeking->SetTimeFormat(&TIME_FORMAT_MEDIA_TIME))) pMediaSeeking->GetDuration(&durationus);
        *totalDuration = durationus/10000000.0;

		if (SUCCEEDED(pMediaSeeking->SetTimeFormat(&TIME_FORMAT_FRAME)) && SUCCEEDED(pMediaSeeking->GetDuration(&duration)))
		{
			if (stopForced) *nrFramesTotal = duration; // if we stopped early, calculate nrFramesTotal
			if (*rate == 0 && durationus) *rate = duration/(durationus/10000000.0);
		} else {
			 // if we stopped early, calculate nrFramesTotal
			if (stopForced && durationus && *rate) *nrFramesTotal = -(*rate)*(durationus/10000000.0);
		}
        pMediaSeeking->Release();
	}

	// some things don't work right if rate is 0
	if (*rate == 0) *rate = 1;

	return S_OK;
}

HRESULT DDGrabber::getAudioInfo(unsigned int id, int* nrChannels, double* rate, int* bits, int* nrFramesCaptured, int* nrFramesTotal, GUID* subtype, double* totalDuration)
{
	if (!nrChannels || !rate || !bits || !nrFramesCaptured || !nrFramesTotal) return E_POINTER;

	if (id >= AudioCBs.size()) return E_NOINTERFACE;
	CSampleGrabberCB* CB = AudioCBs.at(id);

	if (!CB) return E_POINTER;

	if (!CB->pbFormat)
	{
		*nrChannels = 0;
		*rate = 0;
		*bits = 0;
	} else {
		WAVEFORMATEX * h = (WAVEFORMATEX*) CB->pbFormat;
		if (!h) return E_POINTER;
		*nrChannels = h->nChannels;
		*rate = h->nSamplesPerSec;
		*bits = h->wBitsPerSample;

		*subtype = CB->subtype;
	}
	*nrFramesCaptured = CB->frames.size();
	*nrFramesTotal = CB->frameNr;

    *totalDuration = 0;
    
	IMediaSeeking* pMediaSeeking;
	pGraphBuilder->QueryInterface(IID_IMediaSeeking, (void **)&pMediaSeeking);
	if (pMediaSeeking)
	{
		LONGLONG durationus = 0;

		if (SUCCEEDED(pMediaSeeking->SetTimeFormat(&TIME_FORMAT_MEDIA_TIME))) pMediaSeeking->GetDuration(&durationus);
        *totalDuration = durationus/10000000.0;
        pMediaSeeking->Release();
	}
    
    return S_OK;
}

void DDGrabber::getCaptureInfo(int* nrVideo, int* nrAudio)
{
	if (!nrVideo || !nrAudio) return;

	*nrVideo = (VideoCBs.size()>0&&VideoCBs.at(0) && VideoCBs.at(0)->disabled)?0:VideoCBs.size();
	*nrAudio = (AudioCBs.size()>0&&AudioCBs.at(0) && AudioCBs.at(0)->disabled)?0:AudioCBs.size();
}

// data must be freed by caller
HRESULT DDGrabber::getVideoFrame(unsigned int id, int frameNr, BYTE** data, int* nrBytes, double* time)
{
	if (!data || !nrBytes) return E_POINTER;

	if (id >= VideoCBs.size()) return E_NOINTERFACE;
	CSampleGrabberCB* CB = VideoCBs.at(id);
	if (!CB) return E_POINTER;
	if (CB->frameNr == 0) return E_NOINTERFACE;
	if (frameNr < 0 || frameNr >= CB->frames.size()) return E_NOINTERFACE;

	BYTE* tmp = (BYTE*)CB->frames.at(frameNr);
	if (!tmp) return E_NOINTERFACE;

	*nrBytes = CB->frameBytes.at(frameNr);
	*time = CB->frameTimes.at(frameNr);

	*data = (BYTE*) malloc(*nrBytes);
	if (!*data) return E_OUTOFMEMORY;
	memcpy(*data,tmp,*nrBytes);

	return S_OK;
}

// data must be freed by caller
HRESULT DDGrabber::getAudioFrame(unsigned int id, int frameNr, BYTE** data, int* nrBytes, double* time)
{
	if (!data || !nrBytes) return E_POINTER;

	if (id >= AudioCBs.size()) return E_NOINTERFACE;
	CSampleGrabberCB* CB = AudioCBs.at(id);
	if (!CB) return E_POINTER;
	if (CB->frameNr == 0) return E_NOINTERFACE;
	if (frameNr < 0 || frameNr >= CB->frames.size()) return E_NOINTERFACE;

	BYTE* tmp = (BYTE*)CB->frames.at(frameNr);
	if (!tmp) return E_NOINTERFACE;

	*nrBytes = CB->frameBytes.at(frameNr);
	*time = CB->frameTimes.at(frameNr);

	*data = (BYTE*) malloc(*nrBytes);
	if (!*data) return E_OUTOFMEMORY;
	memcpy(*data,tmp,*nrBytes);

	return S_OK;
}

void DDGrabber::setFrames(int* frameNrs, int nrFrames)
{
	if (!frameNrs) return;
	
	this->frameNrs.clear();
    for (int i=0; i<nrFrames; i++) this->frameNrs.add(frameNrs[i]);

	for (int i=0; i < VideoCBs.size(); i++)
	{
		CSampleGrabberCB* CB = VideoCBs.at(i);
		if (CB)
		{
			CB->frames.clear();
			CB->frameNrs.clear();
			for (int i=0; i<nrFrames; i++) CB->frameNrs.add(frameNrs[i]);
			CB->frameNr = 0;
		}
	}
    
	// the meaning of frames doesn't make much sense for audio...
}


void DDGrabber::setTime(double startTime, double stopTime)
{
	this->startTime = startTime;
    this->stopTime = stopTime;
	frameNrs.clear();

    for (int i=0; i < VideoCBs.size(); i++)
	{
		CSampleGrabberCB* CB = VideoCBs.at(i);
		if (CB)
		{
			CB->frames.clear();
			CB->frameNrs.clear();
			CB->frameNr = 0;
			CB->startTime = startTime;
			CB->stopTime = stopTime;
		}
	}

	for (int i=0; i < AudioCBs.size(); i++)
	{
		CSampleGrabberCB* CB = AudioCBs.at(i);
		if (CB)
		{
			CB->frames.clear();
			CB->frameNrs.clear();
			CB->startTime = startTime;
			CB->stopTime = stopTime;
		}
	}
}

void DDGrabber::disableVideo()
{
	for (int i=0; i < VideoCBs.size(); i++)
	{
		if (VideoCBs.at(i)) VideoCBs.at(i)->disabled = true;
	}
}

void DDGrabber::disableAudio()
{
	for (int i=0; i < AudioCBs.size(); i++)
	{
		if (AudioCBs.at(i)) AudioCBs.at(i)->disabled = true;
	}
}

void DDGrabber::setTrySeeking(int tryseek)
{
    tryseeking = tryseek;
}

HRESULT DDGrabber::changeToNull(IGraphBuilder* pGraphBuilder, IBaseFilter* pRenderer)
{
	HRESULT hr;

	if (!pGraphBuilder || !pRenderer) return E_POINTER;

	IPin* outPin = connectedToInput(pRenderer);
	if (!outPin) return E_NOINTERFACE;

	// Add the Null Renderer filter to the graph.
	IBaseFilter *pNull;
	if (FAILED(hr = ::CoCreateInstance(CLSID_NullRenderer,NULL,CLSCTX_INPROC_SERVER,IID_IBaseFilter, (void**)&pNull))) return hr;
	if (FAILED(hr = pGraphBuilder->AddFilter(pNull, L"NullRender"))) return hr;
	if (FAILED(hr = outPin->Disconnect())) return hr;

	IPin* inPin = getInputPin(pNull);
	if (!inPin) return E_NOINTERFACE;

	#ifdef _DEBUG
		FILTER_INFO info;
		pRenderer->QueryFilterInfo(&info);
		char str[100];
		WideCharToMultiByte( CP_ACP, 0, info.achName, -1, str, 100, NULL, NULL );
		_RPT0(_CRT_WARN,str);
		_RPT0(_CRT_WARN," Removed\n");
	#endif

	if (FAILED(hr = pGraphBuilder->RemoveFilter(pRenderer))) return hr;
	pRenderer->Release();

	hr = pGraphBuilder->Connect(outPin,inPin);
	return hr;
}

HRESULT DDGrabber::insertCapture(IGraphBuilder* pGraphBuilder, IBaseFilter* pRenderer, AM_MEDIA_TYPE* mt, CSampleGrabberCB** grabberCB)
{
	HRESULT hr;

	if (!pGraphBuilder || !pRenderer || !mt || !grabberCB) return E_POINTER;

	IPin* rendererPin = getInputPin(pRenderer);
	IPin* upStreamPin = connectedToInput(pRenderer);
	IPin* grabberInPin = NULL;
	IPin* grabberOutPin = NULL;
	IBaseFilter* pGrabberBaseFilter = NULL;
	ISampleGrabber* pSampleGrabber = NULL;

	if (!upStreamPin || !rendererPin) return E_NOINTERFACE;

	// use the for loop so that we can break
	for (int i=0;i<1;i++)
	{
		_RPT0(_CRT_WARN,"Making Grabber\n");
		if (FAILED(hr = ::CoCreateInstance(CLSID_SampleGrabber,NULL,CLSCTX_INPROC_SERVER,IID_IBaseFilter, (LPVOID *)&pGrabberBaseFilter))) return hr;
		pGrabberBaseFilter->QueryInterface(IID_ISampleGrabber, (void**)&pSampleGrabber);
		if (pSampleGrabber == NULL) return E_NOINTERFACE;
		if (FAILED(hr = pGraphBuilder->AddFilter(pGrabberBaseFilter,L"Grabber"))) break;
		if (FAILED(hr = pSampleGrabber->SetMediaType(mt))) break;
		if (!(grabberInPin = getInputPin(pGrabberBaseFilter))) { hr = E_NOINTERFACE; break;}
		if (!(grabberOutPin = getOutputPin(pGrabberBaseFilter))) { hr = E_NOINTERFACE; break;}
		
		_RPT0(_CRT_WARN,"Disconnecting Pins\n");
		if (FAILED(hr = upStreamPin->Disconnect())) break;
		if (FAILED(hr = rendererPin->Disconnect())) break;

		_RPT0(_CRT_WARN,"Connecting Pins\n");
		if (FAILED(hr = pGraphBuilder->Connect(upStreamPin,grabberInPin))) break;

		_RPT0(_CRT_WARN,"Connecting Renderer Pins\n");
		if (FAILED(hr = pGraphBuilder->Connect(grabberOutPin,rendererPin))) break;

		if (FAILED(hr = pSampleGrabber->SetBufferSamples(FALSE))) break;

		*grabberCB = new CSampleGrabberCB();
		if (!*grabberCB) { hr = E_OUTOFMEMORY; break;}

		_RPT0(_CRT_WARN,"Setting Callback\n");
		if (FAILED(hr = pSampleGrabber->SetCallback( *grabberCB, 1 ))) break;
		if (FAILED(hr = pSampleGrabber->GetConnectedMediaType(mt))) break;

		// I don't know how long this pointer will stay valid... so copy it
		if ((*grabberCB)->pbFormat) free((*grabberCB)->pbFormat);
		(*grabberCB)->pbFormat = (BYTE*)malloc(mt->cbFormat);
		if (!(*grabberCB)->pbFormat) { hr = E_OUTOFMEMORY; break;}
		memcpy((*grabberCB)->pbFormat,mt->pbFormat,mt->cbFormat);
	}

	rendererPin->Release();
	upStreamPin->Release();
	if (grabberOutPin) grabberOutPin->Release();
	if (grabberInPin) grabberInPin->Release();
	if (pGrabberBaseFilter) pGrabberBaseFilter->Release();
	if (pSampleGrabber) pSampleGrabber->Release();

	if (FAILED(hr) && pGrabberBaseFilter)
	{
		_RPT0(_CRT_WARN,"Reconnecting original Pins\n");
		upStreamPin->Disconnect();
		if (grabberInPin) grabberInPin->Disconnect();
		if (grabberOutPin) grabberOutPin->Disconnect();
		rendererPin->Disconnect();
		pGraphBuilder->Connect(upStreamPin,rendererPin);

		_RPT0(_CRT_WARN,"Removing Grabber\n");
		pGraphBuilder->RemoveFilter(pGrabberBaseFilter);
	}

	return hr;
}

HRESULT DDGrabber::insertVideoCapture(IGraphBuilder* pGraphBuilder, IBaseFilter* pRenderer)
{
	AM_MEDIA_TYPE mt = {0};
	CSampleGrabberCB* grabberCB;

	if (!pGraphBuilder || !pRenderer) return E_POINTER;

	mt.majortype = MEDIATYPE_Video;
	mt.subtype = MEDIASUBTYPE_RGB24;
	mt.formattype = FORMAT_VideoInfo;

	_RPT0(_CRT_WARN,"Trying to add a VideoCapture.\n");
	HRESULT hr = insertCapture(pGraphBuilder, pRenderer, &mt, &grabberCB);
	if (SUCCEEDED(hr))
	{
		VideoCBs.add(grabberCB);
		_RPT0(_CRT_WARN,"Added a VideoCapture.\n");

		// if the subtype was changed on us... complain.
		if (mt.subtype != MEDIASUBTYPE_RGB24) return E_NOTIMPL;

		//set the rate and bytes per WORD info
		int width, height, dummy;
        double dummyd;
		int currentCB = VideoCBs.size()-1;

		getVideoInfo(currentCB, &width, &height, &(VideoCBs.at(currentCB)->rate), &dummy, &dummy, &dummyd);
		VideoCBs.at(currentCB)->bytesPerWORD = width*height*3;
	}

	MyFreeMediaType(mt);

	return hr;
}

HRESULT DDGrabber::insertAudioCapture(IGraphBuilder* pGraphBuilder, IBaseFilter* pRenderer)
{
	AM_MEDIA_TYPE mt = {0};
	CSampleGrabberCB* grabberCB;

	if (!pGraphBuilder || !pRenderer) return E_POINTER;

	mt.majortype = MEDIATYPE_Audio;
	mt.subtype = MEDIASUBTYPE_PCM;
	mt.formattype = FORMAT_WaveFormatEx;

	_RPT0(_CRT_WARN,"Trying to add an AudioCapture.\n");
	_RPT0(_CRT_WARN,"Trying MEDIASUBTYPE_PCM\n");
	HRESULT hr = insertCapture(pGraphBuilder, pRenderer, &mt, &grabberCB);
	if (FAILED(hr))
	{
		mt.subtype = MEDIASUBTYPE_IEEE_FLOAT;
		_RPT0(_CRT_WARN,"Trying MEDIASUBTYPE_IEEE_FLOAT\n");
		hr = insertCapture(pGraphBuilder, pRenderer, &mt, &grabberCB);
		if (FAILED(hr))
		{
			GUID null = {0};
			mt.subtype = null;
			_RPT0(_CRT_WARN,"Trying NULL\n");
			hr = insertCapture(pGraphBuilder, pRenderer, &mt, &grabberCB);
		}
	}
	if (SUCCEEDED(hr))
	{
		grabberCB->isAudio = true;
		grabberCB->subtype = mt.subtype;
		AudioCBs.add(grabberCB);
		_RPT0(_CRT_WARN,"Added an AudioCapture.\n");

		//set the rate and bytes per WORD info
		int nrChannels, bits, dummy;
        double dummyd;
		GUID subtype;
		int currentCB = AudioCBs.size()-1;

		getAudioInfo(currentCB, &nrChannels, &(AudioCBs.at(currentCB)->rate), &bits, &dummy, &dummy, &subtype, &dummyd);
		AudioCBs.at(currentCB)->bytesPerWORD = bits/8.0*nrChannels;

		// an unsupported bit depth, nrChannels combination; should only be for 4bit.
		if (AudioCBs.at(currentCB)->bytesPerWORD != bits/8.0*nrChannels) hr = E_NOTIMPL;
	}

	MyFreeMediaType(mt);

	return hr;
}

HRESULT DDGrabber::mangleGraph(IGraphBuilder* pGraphBuilder)
{
	if (!pGraphBuilder) return E_POINTER;

	IEnumFilters* filterList;
	IBaseFilter* filt;
	vector<IBaseFilter*> renderers; // if there is an audio and a video stream then we could have two renderers
	ULONG tmp;
	HRESULT hr;

	if (FAILED(hr = pGraphBuilder->EnumFilters(&filterList))) return hr;
	filterList->Reset();
	while (filterList->Next(1, &filt, &tmp) == S_OK)
	{
		if (isRenderer(filt)) 
		{
			renderers.add(filt);
		}
		else filt->Release();
	}
	filterList->Release();

	for (int i=0; i<renderers.size(); i++)
	{
		if (MEDIATYPE_Video == getMajorType(renderers.at(i)))
		{
			_RPT0(_CRT_WARN,"Inserting Video Capture...\n");
			if (FAILED(hr = insertVideoCapture(pGraphBuilder,renderers.at(i)))) return hr;
			_RPT0(_CRT_WARN,"Inserted Video Capture\n");
		} else if (MEDIATYPE_Audio == getMajorType(renderers.at(i))) {
			_RPT0(_CRT_WARN,"Inserting Audio Capture...\n");
			if (FAILED(hr = insertAudioCapture(pGraphBuilder,renderers.at(i)))) return hr;
			_RPT0(_CRT_WARN,"Inserted Audio Capture\n");
		} else {
			_RPT0(_CRT_WARN,"Renderer type not recognized.\n");
		}

		_RPT0(_CRT_WARN,"Changing to Null...\n");
		if (FAILED(hr = changeToNull(pGraphBuilder,(IBaseFilter*)renderers.at(i)))) return hr;
		_RPT0(_CRT_WARN,"Changed to Null\n");
	}

	return S_OK;
}

HRESULT DDGrabber::buildGraph(char* filename)
{
	if (!filename) return E_POINTER;

	WCHAR uniFilename[MAX_PATH];
	HRESULT hr;

	MultiByteToWideChar(CP_ACP, 0, filename, -1, uniFilename, MAX_PATH);

	// Create the graph builder
	pGraphBuilder = NULL;
	if (FAILED(hr = ::CoCreateInstance(CLSID_FilterGraph, NULL, CLSCTX_INPROC_SERVER, IID_IGraphBuilder, (void**)&pGraphBuilder))) return hr;

	// make sure everything is back to "normal"; cleanUp() also releases the current pGraphBuilder, so we need to make a new one
	cleanUp();
	if (FAILED(hr = ::CoCreateInstance(CLSID_FilterGraph, NULL, CLSCTX_INPROC_SERVER, IID_IGraphBuilder, (void**)&pGraphBuilder))) return hr;

    // Add some Null Renderers to the graph so that we don't have weird hardware dependences.
    IBaseFilter *pNull1, *pNull2;
	if (FAILED(hr = ::CoCreateInstance(CLSID_NullRenderer,NULL,CLSCTX_INPROC_SERVER,IID_IBaseFilter, (void**)&pNull1))) return hr;
	if (FAILED(hr = pGraphBuilder->AddFilter(pNull1, L"NullRender"))) return hr;
	if (FAILED(hr = ::CoCreateInstance(CLSID_NullRenderer,NULL,CLSCTX_INPROC_SERVER,IID_IBaseFilter, (void**)&pNull2))) return hr;
	if (FAILED(hr = pGraphBuilder->AddFilter(pNull2, L"NullRender"))) return hr;
  
    _RPT0(_CRT_WARN,"Rendering File...\n");
	if (FAILED(hr = pGraphBuilder->RenderFile(uniFilename,NULL))) return hr;
	_RPT0(_CRT_WARN,"Rendered File\n");

	IPin* nullPin, *tmpPin;
    nullPin = getInputPin(pNull1);
    if (nullPin->ConnectedTo(&tmpPin) == VFW_E_NOT_CONNECTED)
    {
        pGraphBuilder->RemoveFilter(pNull1);
        pNull1->Release();
    } else tmpPin->Release();
    nullPin->Release();
    nullPin = getInputPin(pNull2);
    if (nullPin->ConnectedTo(&tmpPin) == VFW_E_NOT_CONNECTED)
    {
        pGraphBuilder->RemoveFilter(pNull2);
        pNull2->Release();
    } else tmpPin->Release();
    nullPin->Release();

	_RPT0(_CRT_WARN,"Mangling graph...\n");
	// insert the Capture CBs, and change the Renderers to Null
	if (FAILED(hr = mangleGraph(pGraphBuilder))) return hr;
	_RPT0(_CRT_WARN,"Mangled graph...\n");

    return hr;
}

HRESULT DDGrabber::SaveGraphFile(IGraphBuilder *pGraph, WCHAR *wszPath) 
{
    const WCHAR wszStreamName[] = L"ActiveMovieGraph"; 
    HRESULT hr;
    
    IStorage *pStorage = NULL;
    hr = StgCreateDocfile(
        wszPath,
        STGM_CREATE | STGM_TRANSACTED | STGM_READWRITE | STGM_SHARE_EXCLUSIVE,
        0, &pStorage);
    if(FAILED(hr)) 
    {
        return hr;
    }

    IStream *pStream;
    hr = pStorage->CreateStream(
        wszStreamName,
        STGM_WRITE | STGM_CREATE | STGM_SHARE_EXCLUSIVE,
        0, 0, &pStream);
    if (FAILED(hr)) 
    {
        pStorage->Release();    
        return hr;
    }

    IPersistStream *pPersist = NULL;
    pGraph->QueryInterface(IID_IPersistStream, (void**)&pPersist);
    hr = pPersist->Save(pStream, TRUE);
    pStream->Release();
    pPersist->Release();
    if (SUCCEEDED(hr)) 
    {
        hr = pStorage->Commit(STGC_DEFAULT);
    }
    pStorage->Release();
    return hr;
}

/* This section contains code specifically written for QTracker. The routine    */
/* runMatlabCommand is radically modified. Previously, runMatlabCommand is only */
/* responsible for returning an array of bytes representing the red, green, and */
/* blue pixel values of the current video frame, leaving most of the processing */
/* on the Matlab side. In this modified version, runMatlabCommand converts the  */
/* bytes into double values and uses them to produce the img and img2 arrays    */
/* for each fly chamber.                                                        */

#ifdef MATLAB_MEX_FILE
#define ddDouble( x, y ) mxCreateNumericArray( x, y, mxDOUBLE_CLASS, mxREAL )

void DDGrabber::setMatlabCommand(char * matlabCommand)
{
	this->matlabCommand = matlabCommand;
}

void DDGrabber::setChambers( const mxArray *chambers, const mxArray *corners, bool dot )
{
    /* This routine is called to specify the mean value image and the           */
    /* coordinates of the regions of interest. Based on this information,       */
    /* memory areas are allocated. Micah's program structure does not lend      */
    /* itself for easy memory clean-up, which will be addressed in the next     */
    /* version of this grabber code, which will be done from scratch.           */
    
    double *p;

    isDot = dot;
    p = mxGetPr(corners);
    noofChambers = mxGetM(chambers); 
    pChamber = (struct t_chamber *)malloc( noofChambers * sizeof(struct t_chamber) );

    for( int j = 0; j < noofChambers; j++ ) 
    {
        mxArray *mxa = mxGetCell( chambers, j );
	    pChamber[j].roiDims[0] = mxGetM( mxa );
	    pChamber[j].roiDims[1] = mxGetN( mxa );
        size_t sz = pChamber[j].roiDims[0] * pChamber[j].roiDims[1];
        pChamber[j].pmean = (double *)malloc( sz*sizeof(double) );
        memcpy( pChamber[j].pmean, mxGetPr( mxa ), sz*sizeof(double) );

        pChamber[j].R0 = (int)p[ j*4 + 0 ];
        pChamber[j].R1 = (int)p[ j*4 + 1 ];
        pChamber[j].C0 = (int)p[ j*4 + 2 ];
        pChamber[j].C1 = (int)p[ j*4 + 3 ];
    }
}

void DDGrabber::runMatlabCommand()
{
	if (matlabCommand)
	{
		IMediaControl* pMediaControl;
		pGraphBuilder->QueryInterface(IID_IMediaControl, (void **)&pMediaControl);

        /* The pointers ps and pd refer to the source and destination pointers  */
        /* to the byte and double-valued image arrays, respectively. The two    */
        /* pointers are also used as linear indices to their respective arrays, */
        /* which adopt the video frame and Matlab coordinates, respectively     */
        /* The variables d, M, and A are used for the R, G, B, mean, and some   */
        /* temporary image value. Other variables are used as loop variables    */
        BYTE *ps;
        double *pd;
		int dims[] = {1,1,3};
        register double d[3], M, A;
        register int f, i, j, k, l, r, c, v, R, C, V;
        const char *fields[] = { "img", "img2", "imgs", "bin" };

		int width = 0, height = 0;
		int ExitCode;
		mxArray* plhs[] = {NULL};
		mxArray* prhs[6];
		OAFilterState fs;

        prhs[0] = NULL;
		prhs[1] = NULL;
		mexSetTrapFlag(0);

        /* Pause the render graph so that we can safely transmit to Matlab      */

		pMediaControl->GetState(100,&fs);
		if (fs != State_Stopped) pMediaControl->Pause(); 

		for (i = 0; i < VideoCBs.size(); i++)
		{
			CSampleGrabberCB* CB = VideoCBs.at(i);

			if (CB)
			{
                /* If the buffer is ready, obtain its dimension and prepare all */
                /* necessary data structures as arguments to the Matlab routine */

				if (CB->pbFormat)
				{
					VIDEOINFOHEADER * h = (VIDEOINFOHEADER*) CB->pbFormat;
					width  = h->bmiHeader.biWidth;
					height = h->bmiHeader.biHeight;
				}

				prhs[1] = mxCreateDoubleMatrix(1,1,mxREAL); 
				prhs[2] = mxCreateDoubleMatrix(1,1,mxREAL); 
				prhs[3] = mxCreateDoubleMatrix(1,1,mxREAL); 
				prhs[4] = mxCreateDoubleMatrix(1,1,mxREAL); 
				prhs[5] = mxCreateDoubleMatrix(1,1,mxREAL);
                mxGetPr(prhs[2])[0] = width;
                mxGetPr(prhs[3])[0] = height;

                /* Process all frames in the buffer that are ready to be read   */

				for (f = lastFrames.at(i); f < CB->frames.size(); f++)
				{
					if (CB->frames.at(f))
					{
                        /* This is the beginning of QTracker's custom code      */
						/* Allocate matrices to pass into the Matlab functions  */
						/* We are creating a structure array containing as      */
						/* elements as there are chambers. The fields Image in  */
						/* each element of the structure array are identical:   */
						/* it is the normalized double-precision RGB image.     */ 
						/* The fields img and img2 contain processed images     */
						/* whose dimensions match the respective ROIs           */
                        
						dims[0] = height;
                        dims[1] = width;
					    
                        /* Allocate the matrix for the data variable only once  */
                        /* If the chamber data has been specified, then create  */
                        /* Matlab structures for region of interest pixels      */
                        
                        if( prhs[0] == NULL ) 
                        {
						    prhs[0] = ddDouble( 3, dims ); 
                        }

                        if( noofChambers > 0 && mxGetM(prhs[1]) == 1 )
                        {
				            mxDestroyArray(prhs[1]);
                            prhs[1] = mxCreateStructMatrix( noofChambers, 1, 4, fields );
                            for( j = 0; j < noofChambers; j++ ) 
						    {
                                pChamber[j].mxaimg = ddDouble( 2, pChamber[j].roiDims ); 
                                pChamber[j].pimg = mxGetPr( pChamber[j].mxaimg );
                                pChamber[j].mxaimg2 = ddDouble( 2, pChamber[j].roiDims );
                                pChamber[j].pimg2 = mxGetPr( pChamber[j].mxaimg2 );
							    mxSetField( prhs[1], j, fields[0], pChamber[j].mxaimg );
							    mxSetField( prhs[1], j, fields[1], pChamber[j].mxaimg2 );
                            }
                        }						

					    /* Perform transformation from DirectX coordinate to   */
                        /* Matlab coordinate, and copy the normalized double   */
                        /* precision values into the 'data' variable           */
                        /* The variables R, C, and V are Matlab's row, column, */
                        /* and color coordinates for a 3-D matrix. r, c, and v */
                        /* are used to compute the linear index.               */ 
                        
                        ps = CB->frames.at(f); 
                        pd = mxGetPr(prhs[0]);
						for( r = height-1, R = height, k = 0; 
                             r >= 0; 
                             r--, R-- )
                        {
                            for( c = 0, C = 1; 
                                 c < height*width; 
                                 c += height, C++ )
                            {
                                for( v = 2*height*width, V = 0; 
                                     v >=0 ; 
                                     v -= height*width, V++, k++ )
                                {
                                    d[V] = (double)(ps[k])/255.0;
                                    pd[v + c + r] = d[V];
                                }

                                /* Check if the above pixel lies in any of the */
                                /* chambers, i.e., if R,C is between [R0, R1]  */
                                /* and [C0, C1]. If it is, then update         */
                                
                                for( j = 0; j < noofChambers; j++ ) 
						        {
                                    if( pChamber[j].R0 + CS <= R && pChamber[j].R1 - CS >= R &&
                                        pChamber[j].C0 + CS <= C && pChamber[j].C1 - CS >= C )
                                    {
                                        l = pChamber[j].roiDims[0]*(C-pChamber[j].C0) + (R-pChamber[j].R0);
                                        M = pChamber[j].pmean[l]; 
                                        A = 1.0 - (0.6*d[1] + 0.4*d[2])*M;
                                        pChamber[j].pimg2[l] = A + 0.8*d[1] - 1.0*d[2];
                                        pChamber[j].pimg[l] = isDot ? A + 0.6*d[1] - 1.0*d[2] : A ;

                                    }
                                }
                            }   
                        }

						mxGetPr(prhs[4])[0] = CB->frameNrs.size() == 0 ? f+1 : CB->frameNrs.at(f);
						mxGetPr(prhs[5])[0] = CB->frameTimes.size() == 0 ? f+1 : CB->frameTimes.at(f);

						/* Free the video buffer memory */
                        
						free((BYTE*)CB->frames.at(f));
						CB->frames.assign(f,NULL);
						CB->frameBytes.assign(f,0);

						_RPT3(_CRT_WARN,"mexCallMATLAB  %s f# %d #fs %d\n",matlabCommand,f,CB->frames.size());

						/* Execute matlab callback routine */
                        
						ExitCode = mexCallMATLAB( 0, plhs, 6, prhs, matlabCommand );

						_RPT1(_CRT_WARN,"ErrorCode = %d\n",ExitCode);
					}
				}
				lastFrames.assign(i,f);
				
				/* Manually free Matlab structures */
				mxDestroyArray(prhs[0]);
				mxDestroyArray(prhs[1]);
				mxDestroyArray(prhs[2]);
				mxDestroyArray(prhs[3]);
				mxDestroyArray(prhs[4]);
				mxDestroyArray(prhs[5]);
				prhs[0] = NULL;
				prhs[1] = NULL;
			}
		}

		pMediaControl->GetState(100,&fs);
		if (fs != State_Stopped) pMediaControl->Run(); // restart rendering...

		pMediaControl->Release();
	}
}
#undef ddDouble
#endif


HRESULT DDGrabber::doCapture()
{
	HRESULT hr;
	IMediaControl* pMediaControl;
	IMediaEvent* pMediaEventEx;
	IMediaFilter* pMediaFilter;
    IMediaSeeking* pMediaSeeking;

	pGraphBuilder->QueryInterface(IID_IMediaControl, (void **)&pMediaControl);
	pGraphBuilder->QueryInterface(IID_IMediaEvent, (void **)&pMediaEventEx);
	pGraphBuilder->QueryInterface(IID_IMediaFilter, (void **)&pMediaFilter);
	pGraphBuilder->QueryInterface(IID_IMediaSeeking, (void **)&pMediaSeeking);

	if (pMediaControl == NULL || pMediaEventEx == NULL || pMediaFilter == NULL || pMediaSeeking == NULL) return E_NOINTERFACE;

	_RPT1(_CRT_WARN,"tryseeking: %d\n",tryseeking);
    if (tryseeking)
    {
    	_RPT2(_CRT_WARN,"frameNrs.size(): %d  stopTime: %f\n",tryseeking,(float)stopTime);
        if (frameNrs.size() > 0)
        {
            LONGLONG firstFrame = frameNrs.at(0);
            _RPT0(_CRT_WARN,"Trying SetTimeFormat TIME_FORMAT_FRAME\n");
            if (SUCCEEDED(pMediaSeeking->SetTimeFormat(&TIME_FORMAT_FRAME)) && 
              SUCCEEDED(pMediaSeeking->SetPositions(&firstFrame, AM_SEEKING_AbsolutePositioning,NULL,AM_SEEKING_NoPositioning )))
            {
                _RPT0(_CRT_WARN,"SetTimeFormat SUCCEEDED\n");
                for (int i=0; i < VideoCBs.size(); i++)
                {
                    CSampleGrabberCB* CB = VideoCBs.at(i);
                    for (int j=0; j<frameNrs.size(); j++) CB->frameNrs.assign(j,frameNrs.at(j)-firstFrame+1);
                }
            } else {_RPT0(_CRT_WARN,"SetTimeFormat FAILED\n");}
        } else if (stopTime) {
            REFERENCE_TIME llstartTime = startTime*10000000;
            _RPT0(_CRT_WARN,"Trying SetTimeFormat TIME_FORMAT_MEDIA_TIME\n");
            if (SUCCEEDED(pMediaSeeking->SetTimeFormat(&TIME_FORMAT_MEDIA_TIME)))
            {
                _RPT0(_CRT_WARN,"SetTimeFormat SUCCEEDED\n");
                if (SUCCEEDED(pMediaSeeking->SetPositions(&llstartTime, AM_SEEKING_AbsolutePositioning,NULL,AM_SEEKING_NoPositioning ))) {_RPT0(_CRT_WARN,"SetPositions SUCCEEDED\n");}
            } else {_RPT0(_CRT_WARN,"SetTimeFormat FAILED\n");}
        }
    }
    
	//turn off the clock, so that it will run as fast as possible
	pMediaFilter->SetSyncSource(NULL);

	// Run the graph and wait for completion.
	if (FAILED(hr = pMediaControl->Run())) return hr;

	long evCode = 0;
	bool allDone = false;
	stopForced = false;

#ifdef MATLAB_MEX_FILE
	// initialize the lastFrames vector
	lastFrames.clear();
	for (int i=0; i < VideoCBs.size(); i++) lastFrames.add(0);
#endif

	while (pMediaEventEx->WaitForCompletion(1000, &evCode) == E_ABORT)
	{
		allDone = true;
		for (int i=0; i < VideoCBs.size(); i++)
		{
			if (VideoCBs.at(i) && !VideoCBs.at(i)->disabled && !VideoCBs.at(i)->done) allDone = false;
			_RPT2(_CRT_WARN,"Testing VideoCB: %d allDone: %d \n",i,allDone);
		}
		for (int i=0; i < AudioCBs.size(); i++)
		{
			if (AudioCBs.at(i) && !AudioCBs.at(i)->disabled && !AudioCBs.at(i)->done) allDone = false;
			_RPT2(_CRT_WARN,"Testing AudioCBs: %d allDone: %d \n",i,allDone);
		}

		if (allDone)
		{
			_RPT0(_CRT_WARN,"STOP!!! \n");
			// if all of the video streams are done, then force a stop
			if (FAILED(hr = pMediaControl->Stop())) return hr;
			stopForced = true;
		}

#ifdef MATLAB_MEX_FILE
		runMatlabCommand();
#endif
	}

	// make sure everything has really stopped before returning.
	if (FAILED(hr = pMediaControl->Stop())) return hr;

#ifdef MATLAB_MEX_FILE
	runMatlabCommand();
#endif

	pMediaEventEx->Release();
	pMediaControl->Release();
	pMediaFilter->Release();
	pMediaSeeking->Release();

	return S_OK;
}