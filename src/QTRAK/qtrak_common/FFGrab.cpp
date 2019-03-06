/***************************************************
This is the main Grabber code.  It uses AVbin and ffmpeg
to capture video and audio from video and audio files.
Because of this, mmread and supporting code is now
distributed under the LGPL.  See

The code supports any number of audio or video streams and
is a cross platform solution to replace DDGrab.cpp.

This code was intended to be used inside of a matlab interface,
but can be used as a generic grabber class for anyone who needs
one.

Copyright 2008 Micah Richert

This file is part of mmread.

mmread is free software; you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 3 of
the License, or (at your option) any later version.

mmread is distributed WITHOUT ANY WARRANTY.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public
License along with mmread.  If not, see <http://www.gnu.org/licenses/>.
**************************************************/

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#define FFprintf(...) mexPrintf(__VA_ARGS__)
#else
#define FFprintf(...) printf(__VA_ARGS__)
#endif

#include <stdlib.h>
#include <string.h>
#include <vector>
#include <map>

using namespace std;

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "FFGrab.h"
#include "avbin.h"

map<unsigned int,double> keyframes;
unsigned int startDecodingAt;


int Grabber::Grab(AVbinPacket *packet)
{
    if (done) return 0;

    frameNr++;
    packetNr++;
    if (DEBUG_LEVEL > 0) FFprintf("frameNr %d %d\n",frameNr,packetNr);
    int offset = 0, len = 0;
    double timestamp = (packet->timestamp - start_time) / 1000.0 / 1000.0;
    if (DEBUG_LEVEL > 0) FFprintf("time %lld %lld %lf\n",packet->timestamp, start_time, timestamp);

    // either no frames are specified (capture all), or we have time specified
    if (stopTime)
    {
        if (isAudio)
        {
            // time is being used...
            if (timestamp >= startTime)
            {
                // if we've reached the start...
                int bytesPerWord = info.audio.sample_bits * info.audio.channels / 8;
                offset = max(0, ((int)((timestamp - startTime) * rate)) * bytesPerWord);
                len = ((int)((stopTime - timestamp) * rate)) * bytesPerWord;
                // if we have gone past our stop time...
                done = len < 0;
            }
        } else {
            done = stopTime <= timestamp;
            len = (startTime <= timestamp) ? 0x7FFFFFFF : 0;
            if (DEBUG_LEVEL) FFprintf("startTime: %lf, stopTime: %lf, current: %lf, done: %d, len: %d\n",startTime,stopTime,timestamp,done,len);
        }
    } else {
        // capture everything... video or audio
        len = 0x7FFFFFFF;
    }

    if (isAudio)
    {
        return grabAudioPacket(packet, timestamp, len);
    }
    else
    {
        return grabVideoPacket(packet, timestamp, len);
    }
    return 0;
}

int Grabber::grabAudioPacket(AVbinPacket *packet, double from_timestamp, int capture_length)
{
    if (trySeeking && (len<=0 || done)) return 0;

    // !!!!! TODO
    return 0;
}

int Grabber::grabVideoPacket(AVbinPacket *packet, double from_timestamp, int capture_length)
{
    bool skip = false;
    if (frameNrs.size() > 0)
    {
        //frames are being specified
        // check to see if the frame is in our list
        bool foundNr = false;
        unsigned int lastFrameNr = 0;
        for (int i=0; i < frameNrs.size(); i++)
        {
            if (frameNrs.at(i) == frameNr) foundNr = true;
            if (frameNrs.at(i) > lastFrameNr) lastFrameNr = frameNrs.at(i);
        }
        done = frameNr > lastFrameNr;
        if (!foundNr) {
            if (DEBUG_LEVEL > 0) FFprintf("Skipping frame %d\n", frameNr);
            skip = true;
        }
    }
    if ((trySeeking && skip && packetNr < startDecodingAt && packetNr != 1) || done) return 0;

    // allocate buffer for the 
    int numBytes = av_image_get_buffer_size(AV_PIX_FMT_RGB24, 
                                            stream->codec_context->width,
                                            stream->codec_context->height,
                                            32);
    uint8_t* videobuf = (uint8_t *) av_malloc(numBytes * sizeof(uint8_t));

    if (videobuf == NULL) {
        FFprintf("Error allocating %d bytes for videobuf\n", numBytes);
        return 2;
    }
    if (DEBUG_LEVEL > 0) FFprintf("avbin_decode_video %d to %x\n", packet->size, videobuf);

    int nBytesRead = avbin_decode_video(stream, packet);
    if (nBytesRead <= 0)
    {
        FFprintf("avbin_decode_video FAILED!\n");
        // silently ignore decode errors
        frameNr--;
        free(videobuf);
        return 3;
    }
    if (stream->frame->key_frame)
    {
        keyframes[packetNr] = from_timestamp;
    }

    if (skip || capture_length == 0)
    {
        if (DEBUG_LEVEL > 0) FFprintf("Free videobuf %x\n", videobuf);
        av_free(videobuf);
        return 0;
    } else {
        size_t buflen = min(capture_length, numBytes);
        if (DEBUG_LEVEL > 0) FFprintf("Push videobuf %f: %d bytes to %x\n", from_timestamp, buflen, videobuf);
        frames.push_back(videobuf);
        frameBytes.push_back(buflen);
        frameTimes.push_back(from_timestamp);
    }
    return 0;
}

FFGrabber::FFGrabber()
{
    stopForced = false;
    tryseeking = true;
    file = NULL;
    filename = NULL;

#ifdef MATLAB_MEX_FILE
    isDot = false;
    noofChambers = 0;
#endif

    if (avbin_init() != AVBIN_RESULT_OK) FFprintf("avbin_init init failed!\n");
    avbin_set_log_level(DEBUG_LEVEL);
}

void FFGrabber::cleanUp()
{
    if (!file) return; // nothing to cleanup.

    for (streammap::iterator i = streams.begin(); i != streams.end(); i++)
    {
        avbin_close_stream(i->second->stream);
        delete i->second;
    }

    streams.clear();
    videos.clear();
    audios.clear();

    avbin_close_file(file);
    file = NULL;

#ifdef MATLAB_MEX_FILE
    if (matlabCommand) free(matlabCommand);
    matlabCommand = NULL;
#endif
}

int FFGrabber::getVideoInfo(unsigned int id, int* width, int* height, double* rate, int* nrFramesCaptured, int* nrFramesTotal, double* totalDuration)
{
    if (!width || !height || !nrFramesCaptured || !nrFramesTotal) return -1;

    if (id >= videos.size()) return -2;
    Grabber* CB = videos.at(id);

    if (!CB) return -1;

    *width  = CB->info.video.width;
    *height = CB->info.video.height;
    *rate = CB->rate;
    *nrFramesCaptured = CB->frames.size();
    *nrFramesTotal = CB->frameNr;

    *totalDuration = fileinfo.duration/1000.0/1000.0;
    if (stopForced) *nrFramesTotal = (int)(-(*rate)*(*totalDuration));

    return 0;
}

int FFGrabber::getAudioInfo(unsigned int id, int* nrChannels, double* rate, int* bits, int* nrFramesCaptured, int* nrFramesTotal, int* subtype, double* totalDuration)
{
    if (!nrChannels || !rate || !bits || !nrFramesCaptured || !nrFramesTotal) return -1;

    if (id >= audios.size()) return -2;
    Grabber* CB = audios.at(id);

    if (!CB) return -1;

    *nrChannels = CB->info.audio.channels;
    *rate = CB->info.audio.sample_rate;
    *bits = CB->info.audio.sample_bits;
    *subtype = CB->info.audio.sample_format;
    *nrFramesCaptured = CB->frames.size();
    *nrFramesTotal = CB->frameNr;

    *totalDuration = fileinfo.duration/1000.0/1000.0;

    return 0;
}

void FFGrabber::getCaptureInfo(int* nrVideo, int* nrAudio)
{
    if (!nrVideo || !nrAudio) return;

    *nrVideo = videos.size();
    *nrAudio = audios.size();
}

// data must be freed by caller
int FFGrabber::getVideoFrame(unsigned int id, unsigned int frameNr, uint8_t** data, unsigned int* nrBytes, double* time)
{
    if (DEBUG_LEVEL > 0) FFprintf("getting Video frame %d\n",frameNr);

    if (!data || !nrBytes) return -1;

    if (id >= videos.size()) return -2;
    Grabber* CB = videos[id];
    if (!CB) return -1;
    if (CB->frameNr == 0) return -2;
    if (frameNr < 0 || frameNr >= CB->frames.size()) return -2;

    uint8_t* tmp = CB->frames[frameNr];
    if (!tmp) return -2;

    *nrBytes = CB->frameBytes[frameNr];
    *time = CB->frameTimes[frameNr];

    *data = tmp;
    CB->frames[frameNr] = NULL;

    return 0;
}

// data must be freed by caller
int FFGrabber::getAudioFrame(unsigned int id, unsigned int frameNr, uint8_t** data, unsigned int* nrBytes, double* time)
{
    if (!data || !nrBytes) return -1;

    if (id >= audios.size()) return -2;
    Grabber* CB = audios[id];
    if (!CB) return -1;
    if (CB->frameNr == 0) return -2;
    if (frameNr < 0 || frameNr >= CB->frames.size()) return -2;

    uint8_t* tmp = CB->frames[frameNr];
    if (!tmp) return -2;

    *nrBytes = CB->frameBytes[frameNr];
    *time = CB->frameTimes[frameNr];

    *data = tmp;
    CB->frames[frameNr] = NULL;

    return 0;
}

void FFGrabber::setFrames(unsigned int* frameNrs, int nrFrames)
{
    if (!frameNrs) return;

    unsigned int minFrame=nrFrames>0?frameNrs[0]:0;

    this->frameNrs.clear();
    for (int i=0; i < nrFrames; i++) this->frameNrs.push_back(frameNrs[i]);

    for (int j=0; j < videos.size(); j++)
    {
        Grabber* CB = videos.at(j);
        if (CB)
        {
            CB->frames.clear();
            CB->frameNrs.clear();
            for (int i=0; i<nrFrames; i++)
            {
                CB->frameNrs.push_back(frameNrs[i]);
                minFrame=frameNrs[i]<minFrame?frameNrs[i]:minFrame;
            }
            CB->frameNr = 0;
            CB->packetNr = 0;
        }
    }

    if (tryseeking && nrFrames > 0)
    {
        startDecodingAt = 0;
        for (map<unsigned int,double>::const_iterator it=keyframes.begin(); it != keyframes.end(); it++)
        {
            if (it->first <= minFrame && it->first > startDecodingAt) startDecodingAt = it->first;
            if (DEBUG_LEVEL > 0) FFprintf("%d %d\n",it->first,startDecodingAt);
        }
    }

    // the meaning of frames doesn't make much sense for audio...
}


void FFGrabber::setTime(double startTime, double stopTime)
{
    this->startTime = startTime;
    this->stopTime = stopTime;
    frameNrs.clear();

    for (int i=0; i < videos.size(); i++)
    {
        Grabber* CB = videos.at(i);
        if (CB)
        {
            CB->frames.clear();
            CB->frameNrs.clear();
            CB->frameNr = 0;
            CB->packetNr = 0;
            CB->startTime = startTime;
            CB->stopTime = stopTime;
        }
    }

    for (int i=0; i < audios.size(); i++)
    {
        Grabber* CB = audios.at(i);
        if (CB)
        {
            CB->frames.clear();
            CB->frameNrs.clear();
            CB->startTime = startTime;
            CB->stopTime = stopTime;
        }
    }
}

#ifdef MATLAB_MEX_FILE


void FFGrabber::setMatlabCommand(char * matlabCommand)
{
    this->matlabCommand = matlabCommand;
}

void FFGrabber::setChambers( const mxArray *chambers, const mxArray *corners, bool dot )
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

void FFGrabber::runMatlabCommand(Grabber* G)
{
    if (matlabCommand)
    {
        /* The pointers ps and pd refer to the source and destination pointers  */
        /* to the byte and double-valued image arrays, respectively. The two    */
        /* pointers are also used as linear indices to their respective arrays, */
        /* which adopt the video frame and Matlab coordinates, respectively     */
        /* The variables d, M, and A are used for the R, G, B, mean, and some   */
        /* temporary image value. Other variables are used as loop variables    */
        uint8_t *ps;
        double *pd;
        mwSize dims[] = {1,1,3};
        double d[3], M, A;
        int f, i, j, k, l, r, c, v, R, C, V;
        const char *fields[] = { "img", "img2", "imgs", "bin" };

        int width=G->info.video.width, height = G->info.video.height;
        int ExitCode;
        mxArray* plhs[] = {NULL};

        // mexSetTrapFlag(0);

        if (G->frames.size() == 0) return;
        vector<uint8_t*>::iterator lastframe = --(G->frames.end());

        if (*lastframe == NULL) return;

        dims[0] = height;
        dims[1] = width;

        /* Allocate the matrix for the data variable only once  */
        /* If the chamber data has been specified, then create  */
        /* Matlab structures for region of interest pixels      */

        if( prhs[0] == NULL ) 
        {
            prhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
                                        
        
            prhs[1] = mxCreateDoubleMatrix(1,1,mxREAL); 
            prhs[2] = mxCreateDoubleMatrix(1,1,mxREAL); 
            prhs[3] = mxCreateDoubleMatrix(1,1,mxREAL); 
            prhs[4] = mxCreateDoubleMatrix(1,1,mxREAL); 
            prhs[5] = mxCreateDoubleMatrix(1,1,mxREAL);
            mxGetPr(prhs[2])[0] = width;
            mxGetPr(prhs[3])[0] = height;
        }

        /* This is the beginning of QTracker's custom code      */
        /* Allocate matrices to pass into the Matlab functions  */
        /* We are creating a structure array containing as      */
        /* elements as there are chambers. The fields Image in  */
        /* each element of the structure array are identical:   */
        /* it is the normalized double-precision RGB image.     */ 
        /* The fields img and img2 contain processed images     */
        /* whose dimensions match the respective ROIs           */



        if( noofChambers > 0 && mxGetM(prhs[1]) == 1 )
        {
            mxDestroyArray(prhs[1]);
            prhs[1] = mxCreateStructMatrix( noofChambers, 1, 4, fields );
            for( j = 0; j < noofChambers; j++ ) 
            {
                pChamber[j].mxaimg = mxCreateNumericArray(2, pChamber[j].roiDims, mxDOUBLE_CLASS, mxREAL); 
                pChamber[j].pimg = mxGetPr( pChamber[j].mxaimg );
                pChamber[j].mxaimg2 = mxCreateNumericArray(2, pChamber[j].roiDims, mxDOUBLE_CLASS, mxREAL);
                pChamber[j].pimg2 = mxGetPr( pChamber[j].mxaimg2 );
                mxSetField( prhs[1], j, fields[0], pChamber[j].mxaimg );
                mxSetField( prhs[1], j, fields[1], pChamber[j].mxaimg2 );
            }
        }

        /* Perform transformation from FFMPEG coordinate to   */
        /* Matlab coordinate, and copy the normalized double   */
        /* precision values into the 'data' variable           */
        /* The variables R, C, and V are Matlab's row, column, */
        /* and color coordinates for a 3-D matrix. r, c, and v */
        /* are used to compute the linear index.               */ 

        ps = *lastframe; 
        pd = mxGetPr(prhs[0]);

        for( r = 0, R = 1, k = 0; r <= height-1; r++, R++ )
        {
            for( c = 0, C = 1; c < height*width; c += height, C++ )
            {
                for( v = 2*height*width, V = 0;  v >=0 ;  v -= height*width, V++, k++ )
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

        mxGetPr(prhs[4])[0] = G->frameNrs.size()==0?G->frameTimes.size():G->frameNrs[G->frameTimes.size()-1];
        mxGetPr(prhs[5])[0] = G->frameTimes.back();

        /* Free the video buffer memory */
        free(*lastframe);
        *lastframe = NULL;


        /* Execute matlab callback routine */

        ExitCode = mexCallMATLAB( 0, plhs, 6, prhs, matlabCommand );
    }
}
#undef ddDouble

#endif

int FFGrabber::build(const char* filename, bool disableVideo, bool disableAudio, bool tryseeking)
{
    if (DEBUG_LEVEL > 0) FFprintf("avbin_open_filename %s %s\n", filename);
    file = avbin_open_filename(filename);
    if (!file) {
        if (DEBUG_LEVEL > 0) FFprintf("Error in avbin_open_filename %s\n", filename);
        return -4;
    }

    //detect if the file has changed
    struct stat fstat;
    stat(filename, &fstat);

    if (!this->filename || strcmp(this->filename,filename)!=0 || filestat.st_mtime != fstat.st_mtime)
    {
        free(this->filename);
        this->filename=strdup(filename);
        memcpy(&filestat,&fstat,sizeof(fstat));

        keyframes.clear();
        startDecodingAt = 0xFFFFFFFF;
    }

    fileinfo.structure_size = sizeof(fileinfo);

    if (avbin_file_info(file, &fileinfo)) {
        if (DEBUG_LEVEL > 0) FFprintf("ERROR avbin_file_info\n");
        return -1;
    }
    if (DEBUG_LEVEL > 0) FFprintf("avbin_file_info size: %ld\n", fileinfo.structure_size);

    for (int stream_index=0; stream_index<fileinfo.n_streams; stream_index++)
    {
        AVbinStreamInfo streaminfo;
        streaminfo.structure_size = sizeof(streaminfo);

        if (DEBUG_LEVEL > 0) FFprintf("avbin_stream_info\n");
        avbin_stream_info(file, stream_index, &streaminfo);

        if (DEBUG_LEVEL > 0) FFprintf("%lld\n",streaminfo,fileinfo.start_time);

        if (streaminfo.type == AVBIN_STREAM_TYPE_VIDEO && !disableVideo)
        {
            AVbinStream* tmp = avbin_open_stream(file, stream_index);
            if (tmp)
            {
                double rate = streaminfo.video.frame_rate_num/(0.00001+streaminfo.video.frame_rate_den);
                streams[stream_index] = new Grabber(false,
                                                    tmp,
                                                    tryseeking,
                                                    rate,
                                                    streaminfo,
                                                    fileinfo.start_time);
                videos.push_back(streams[stream_index]);
            }
            else {
                FFprintf("Could not open video stream\n");
            }
        }
        if (streaminfo.type == AVBIN_STREAM_TYPE_AUDIO && !disableAudio)
        {
            AVbinStream* tmp = avbin_open_stream(file, stream_index);
            if (tmp)
            {
                streams[stream_index] = new Grabber(true,
                                                    tmp,
                                                    tryseeking,
                                                    streaminfo.audio.sample_rate,
                                                    streaminfo,
                                                    fileinfo.start_time);
                audios.push_back(streams[stream_index]);
            } else {
                FFprintf("Could not open audio stream\n");
            }
        }
    }
    this->tryseeking = tryseeking;
    stopForced = false;

    return 0;
}

int FFGrabber::doCapture()
{
    AVbinPacket packet;
    packet.packet = av_packet_alloc();
    if (packet.packet == NULL) {
        // error allocating the packet
        return AVBIN_RESULT_ERROR;
    }
    packet.structure_size = sizeof(packet);
    streammap::iterator tmp;
    int needseek=1;
    bool allDone = false;
    while (!avbin_read_next_packet(file, &packet))
    {
        if (DEBUG_LEVEL > 0) FFprintf("AVbinPacket: (%ld: %d %ld)\n", packet.timestamp, packet.stream_index, packet.size);

        if ((tmp = streams.find(packet.stream_index)) != streams.end())
        {
            Grabber* G = tmp->second;
            G->Grab(&packet);

            if (G->done)
            {
                    allDone = true;
                    for (streammap::iterator i = streams.begin(); i != streams.end() && allDone; i++)
                    {
                            allDone = allDone && i->second->done;
                    }
            }

#ifdef MATLAB_MEX_FILE
            if (!G->isAudio) runMatlabCommand(G);
#endif
        }

        if (tryseeking && needseek)
        {
            if (stopTime && startTime > 0) {
                    if (DEBUG_LEVEL > 0) FFprintf("try seeking to %lf\n",startTime);
                    av_seek_frame(file->format_context, -1, (AVbinTimestamp)(startTime*1000*1000), AVSEEK_FLAG_BACKWARD);
            }
            needseek = 0;
        }

        if (allDone)
        {
            if (DEBUG_LEVEL > 0) FFprintf("stopForced\n");
            stopForced = true;
            break;
        }
    }

#ifdef MATLAB_MEX_FILE
    if (prhs[0])
    {
        mxDestroyArray(prhs[0]);
        if (prhs[1]) mxDestroyArray(prhs[1]);
        if (prhs[2]) mxDestroyArray(prhs[2]);
        if (prhs[3]) mxDestroyArray(prhs[3]);
        if (prhs[4]) mxDestroyArray(prhs[4]);
    }
    prhs[0] = NULL;
#endif

    return 0;
}

#ifdef MATLAB_MEX_FILE
FFGrabber FFG;

const char* message(int err)
{
    switch (err)
    {
        case 0: return "";
        case -1: return "Unable to initialize";
        case -2: return "Invalid interface";
        case -4: return "Unable to open file";
        case -5: return "AVbin version 8 or greater is required!";
        default: return "Unknown error";
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs < 1 || !mxIsChar(prhs[0])) mexErrMsgTxt("First parameter must be the command (a string)");

    char cmd[100];
    mxGetString(prhs[0],cmd,100);

    if (!strcmp("build",cmd)) {
    	if (nrhs < 5 || !mxIsChar(prhs[1])) mexErrMsgTxt("build: parameters must be the filename (as a string), disableVideo, disableAudio, trySeeking");
    	if (nlhs > 0) mexErrMsgTxt("build: there are no outputs");
    	int filenamelen = mxGetN(prhs[1])+1;
        char* filename = new char[filenamelen];
    	if (!filename) mexErrMsgTxt("build: out of memory");
    	mxGetString(prhs[1],filename,filenamelen);
    	const char* errmsg =  message(FFG.build(filename, mxGetScalar(prhs[2]), mxGetScalar(prhs[3]), mxGetScalar(prhs[4])));
    	delete[] filename;
    	if (strcmp("",errmsg)) mexErrMsgTxt(errmsg);
    } else if (!strcmp("doCapture",cmd)) {
        if (nlhs > 0) mexErrMsgTxt("doCapture: there are no outputs");
        const char* errmsg =  message(FFG.doCapture());
        if (strcmp("",errmsg)) mexErrMsgTxt(errmsg);
    } else if (!strcmp("getVideoInfo",cmd)) {
        if (nrhs < 2 || !mxIsNumeric(prhs[1])) mexErrMsgTxt("getVideoInfo: second parameter must be the video stream id (as a number)");
        if (nlhs > 6) mexErrMsgTxt("getVideoInfo: there are only 5 output values: width, height, rate, nrFramesCaptured, nrFramesTotal");

        unsigned int id = (unsigned int)mxGetScalar(prhs[1]);
        int width,height,nrFramesCaptured,nrFramesTotal;
        double rate, totalDuration;
        const char* errmsg =  message(FFG.getVideoInfo(id, &width, &height,&rate, &nrFramesCaptured, &nrFramesTotal, &totalDuration));

        if (strcmp("",errmsg)) mexErrMsgTxt(errmsg);

        if (nlhs >= 1) {plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[0])[0] = width; }
        if (nlhs >= 2) {plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[1])[0] = height; }
        if (nlhs >= 3) {plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[2])[0] = rate; }
        if (nlhs >= 4) {plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[3])[0] = nrFramesCaptured; }
        if (nlhs >= 5) {plhs[4] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[4])[0] = nrFramesTotal; }
        if (nlhs >= 6) {plhs[5] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[5])[0] = totalDuration; }
    } else if (!strcmp("getAudioInfo",cmd)) {
        if (nrhs < 2 || !mxIsNumeric(prhs[1])) mexErrMsgTxt("getAudioInfo: second parameter must be the audio stream id (as a number)");
        if (nlhs > 7) mexErrMsgTxt("getAudioInfo: there are only 6 output values: nrChannels, rate, bits, nrFramesCaptured, nrFramesTotal, subtype");

        unsigned int id = (unsigned int)mxGetScalar(prhs[1]);
        int nrChannels,bits,nrFramesCaptured,nrFramesTotal,subtype;
        double rate, totalDuration;
        const char* errmsg =  message(FFG.getAudioInfo(id, &nrChannels, &rate, &bits, &nrFramesCaptured, &nrFramesTotal, &subtype, &totalDuration));

        if (strcmp("",errmsg)) mexErrMsgTxt(errmsg);

        if (nlhs >= 1) {plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[0])[0] = nrChannels; }
        if (nlhs >= 2) {plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[1])[0] = rate; }
        if (nlhs >= 3) {plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[2])[0] = bits; }
        if (nlhs >= 4) {plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[3])[0] = nrFramesCaptured; }
        if (nlhs >= 5) {plhs[4] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[4])[0] = nrFramesTotal; }
        if (nlhs >= 6) {plhs[5] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[5])[0] = subtype==AVBIN_SAMPLE_FORMAT_FLOAT?1:0; }
        if (nlhs >= 7) {plhs[6] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[6])[0] = totalDuration; }
    } else if (!strcmp("getCaptureInfo",cmd)) {
        if (nlhs > 2) mexErrMsgTxt("getCaptureInfo: there are only 2 output values: nrVideo, nrAudio");

        int nrVideo, nrAudio;
        FFG.getCaptureInfo(&nrVideo, &nrAudio);

        if (nlhs >= 1) {plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[0])[0] = nrVideo; }
        if (nlhs >= 2) {plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[1])[0] = nrAudio; }
    } else if (!strcmp("getVideoFrame",cmd)) {
        if (nrhs < 3 || !mxIsNumeric(prhs[1]) || !mxIsNumeric(prhs[2])) mexErrMsgTxt("getVideoFrame: second parameter must be the audio stream id (as a number) and third parameter must be the frame number");
        if (nlhs > 2) mexErrMsgTxt("getVideoFrame: there are only 2 output value: data");

        unsigned int id = (unsigned int)mxGetScalar(prhs[1]);
        unsigned int frameNr = (unsigned int)mxGetScalar(prhs[2]);
        uint8_t* data;
        unsigned int nrBytes;
        double time;
        mwSize dims[2];
        dims[1]=1;
        const char* errmsg =  message(FFG.getVideoFrame(id, frameNr, &data, &nrBytes, &time));

        if (strcmp("",errmsg)) mexErrMsgTxt(errmsg);

        dims[0] = nrBytes;
        plhs[0] = mxCreateNumericArray(2, dims, mxUINT8_CLASS, mxREAL); // empty 2d matrix
        memcpy(mxGetPr(plhs[0]),data,nrBytes);
        free(data);
        if (nlhs >= 2) {plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[1])[0] = time; }
    } else if (!strcmp("getAudioFrame",cmd)) {
        if (nrhs < 3 || !mxIsNumeric(prhs[1]) || !mxIsNumeric(prhs[2])) mexErrMsgTxt("getAudioFrame: second parameter must be the audio stream id (as a number) and third parameter must be the frame number");
        if (nlhs > 2) mexErrMsgTxt("getAudioFrame: there are only 2 output value: data");

        unsigned int id = (unsigned int)mxGetScalar(prhs[1]);
        unsigned int frameNr = (unsigned int)mxGetScalar(prhs[2]);
        uint8_t* data;
        unsigned int nrBytes;
        double time;
        mwSize dims[2];
        dims[1]=1;
        mxClassID mxClass;
        const char* errmsg =  message(FFG.getAudioFrame(id, frameNr, &data, &nrBytes, &time));

        if (strcmp("",errmsg)) mexErrMsgTxt(errmsg);

        int nrChannels,bits,nrFramesCaptured,nrFramesTotal,subtype;
        double rate, totalDuration;
        FFG.getAudioInfo(id, &nrChannels, &rate, &bits, &nrFramesCaptured, &nrFramesTotal, &subtype, &totalDuration);

        switch (bits)
        {
            case 8:
            {
                dims[0] = nrBytes;
                mxClass = mxUINT8_CLASS;
                break;
            }
            case 16:
            {
                mxClass = mxINT16_CLASS;
                dims[0] = nrBytes/2;
                break;
            }
            case 24:
            {
                int* tmpdata = (int*)malloc(nrBytes/3*4);
                int i;

                //I don't know how 24bit float data is organized...
                for (i=0;i<nrBytes/3;i++)
                {
                        tmpdata[i] = (((0x80&data[i*3+2])?-1:0)&0xFF000000) | ((data[i*3+2]<<16)+(data[i*3+1]<<8)+data[i*3]);
                }

                free(data);
                data = (uint8_t*)tmpdata;

                mxClass = mxINT32_CLASS;
                dims[0] = nrBytes/3;
                nrBytes = nrBytes/3*4;
                break;
            }
            case 32:
            {
                mxClass = subtype==AVBIN_SAMPLE_FORMAT_S32?mxINT32_CLASS:subtype==AVBIN_SAMPLE_FORMAT_FLOAT?mxSINGLE_CLASS:mxUINT32_CLASS;
                dims[0] = nrBytes/4;
                break;
            }
            default:
            {
                dims[0] = nrBytes;
                mxClass = mxUINT8_CLASS;
                break;
            }
        }

        plhs[0] = mxCreateNumericArray(2, dims, mxClass, mxREAL); // empty 2d matrix
        memcpy(mxGetPr(plhs[0]),data,nrBytes);
        free(data);
        if (nlhs >= 2) {plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL); mxGetPr(plhs[1])[0] = time; }
    } else if (!strcmp("setFrames",cmd)) {
        if (nrhs < 2 || !mxIsDouble(prhs[1])) mexErrMsgTxt("setFrames: second parameter must be the frame numbers (as doubles)");
        if (nlhs > 0) mexErrMsgTxt("setFrames: has no outputs");
        int nrFrames = mxGetN(prhs[1]) * mxGetM(prhs[1]);
        unsigned int* frameNrs = new unsigned int[nrFrames];
        if (!frameNrs) mexErrMsgTxt("setFrames: out of memory");
        double* data = mxGetPr(prhs[1]);
        for (int i=0; i<nrFrames; i++) frameNrs[i] = (unsigned int)data[i];

        FFG.setFrames(frameNrs, nrFrames);

        delete[] frameNrs;

    // QTracker specific: specify mean image and ROI coordinates. 
    } else if (!strcmp("setChambers",cmd)) {
        if( nrhs < 4 || !mxIsCell(prhs[1]) || !mxIsDouble(prhs[2]) )
            mexErrMsgTxt("setChambers: Mean images with ROIs and corner points required!");
        FFG.setChambers( prhs[1], prhs[2], (bool)*(mxGetPr(prhs[3])) ); 

    } else if (!strcmp("setTime",cmd)) {
        if (nrhs < 3 || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2])) mexErrMsgTxt("setTime: start and stop time are required (as doubles)");
        if (nlhs > 0) mexErrMsgTxt("setTime: has no outputs");

        FFG.setTime(mxGetScalar(prhs[1]), mxGetScalar(prhs[2]));
    } else if (!strcmp("setMatlabCommand",cmd)) {
        if (nrhs < 2 || !mxIsChar(prhs[1])) mexErrMsgTxt("setMatlabCommand: the command must be passed as a string");
        if (nlhs > 0) mexErrMsgTxt("setMatlabCommand: has no outputs");

        char * matlabCommand = (char*)calloc(100,1);
        mxGetString(prhs[1],matlabCommand,100);

        if (strlen(matlabCommand)==0)
        {
                FFG.setMatlabCommand(NULL);
                free(matlabCommand);
        } else FFG.setMatlabCommand(matlabCommand);
    } else if (!strcmp("cleanUp",cmd)) {
        if (nlhs > 0) mexErrMsgTxt("cleanUp: there are no outputs");
        FFG.cleanUp();
    }
}
#endif

#ifdef TEST_FFGRAB
int main(int argc, char** argv)
{
    FFGrabber FFG;
    printf("%s\n",argv[1]);
    FFG.build(argv[1],false,false,true);
    FFG.doCapture();
    int nrVideo, nrAudio;
    FFG.getCaptureInfo(&nrVideo, &nrAudio);

    printf("there are %d video streams, and %d audio.\n",nrVideo,nrAudio);
}
#endif
