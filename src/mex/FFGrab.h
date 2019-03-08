#ifndef FFGRAB_H
#define FFGRAB_H

#ifndef DEBUG_LEVEL
#define DEBUG_LEVEL 0
#endif

#ifdef __cplusplus
extern "C" {
#endif

#include "avbin.h"


#ifdef __cplusplus
}
#endif
class Grabber
{
public:
    Grabber(bool isAudio,
            AVbinStream* stream,
            bool trySeeking,
            double rate,
            AVbinStreamInfo info,
            int64_t start_time)
    {
        this->stream = stream;
        frameNr = 0;
        packetNr = 0;
        done = false;
        this->rate = rate;
        startTime = 0;
        nrFramesTotal = isAudio ? 0 : info.video.nb_frames;
        stopTime = 0;
        this->isAudio = isAudio;
        this->info = info;
        this->trySeeking = trySeeking;
        this->start_time = start_time > 0 ? start_time : 0;
    };

    ~Grabber()
    {
        // clean up any remaining memory...
        for (vector<uint8_t*>::iterator i=frames.begin();i != frames.end(); i++) av_free(*i);
    }

    AVbinStream* stream;
    AVbinStreamInfo info;
    int64_t start_time;
    int64_t nrFramesTotal;

    vector<uint8_t*> frames;
    vector<unsigned int> frameBytes;
    vector<double> frameTimes;

    vector<unsigned int> frameNrs;

    unsigned int frameNr;
    unsigned int packetNr;
    bool done;
    bool isAudio;
    bool trySeeking;

    double rate;
    double startTime, stopTime;

    int Grab(AVbinPacket* packet);
private:
	int grabAudioPacket(AVbinPacket* packet, double from_timestamp, int capture_length);
	int grabVideoPacket(AVbinPacket* packet, double from_timestamp, int capture_length);
};

typedef map<int,Grabber*> streammap;

class FFGrabber
{
public:
    FFGrabber();

    int build(const char* filename, bool disableVideo, bool disableAudio, bool tryseeking);
    int doCapture();

    int getVideoInfo(unsigned int id, int* width, int* height, double* rate, int* nrFramesCaptured, int* nrFramesTotal, double* totalDuration);
    int getAudioInfo(unsigned int id, int* nrChannels, double* rate, int* bits, int* nrFramesCaptured, int* nrFramesTotal, int* subtype, double* totalDuration);
    void getCaptureInfo(int* nrVideo, int* nrAudio);
    // data must be freed by caller
    int getVideoFrame(unsigned int id, unsigned int frameNr, uint8_t** data, unsigned int* nrBytes, double* time);
    // data must be freed by caller
    int getAudioFrame(unsigned int id, unsigned int frameNr, uint8_t** data, unsigned int* nrBytes, double* time);
    void setFrames(unsigned int* frameNrs, int nrFrames);
    void setTime(double startTime, double stopTime);
    void disableVideo();
    void disableAudio();
    void cleanUp(); // must be called at the end, in order to render anything afterward.

#ifdef MATLAB_MEX_FILE
    void setChambers(const mxArray *chambers, const mxArray *corners, bool dot);
    void setMatlabCommand(char * matlabCommand);
    void runMatlabCommand(Grabber* G);
#endif
private:
    streammap streams;
    vector<Grabber*> videos;
    vector<Grabber*> audios;

    AVbinFile* file;
    AVbinFileInfo fileinfo;

    bool stopForced;
    bool tryseeking;
    vector<unsigned int> frameNrs;
    double startTime, stopTime;

    char* filename;
    struct stat filestat;


#ifdef MATLAB_MEX_FILE
#define CS 0

    char* matlabCommand;
    mxArray* prhs[6];

    bool isDot;                            /* Does one fly have a dot?         */
    int noofChambers;                      /* The number of fly chambers       */
    struct t_chamber                       /* Structure for each chamber       */
    {   
        int R0, R1, C0, C1;                /* ROI corner points                */
        mwSize roiDims[2];                 /* ROI matrix dimension             */
        double *pimg, *pimg2, *pmean;      /* Pointer to pixel values          */
        mxArray *mxaimg, *mxaimg2;         /* Pointer to Matlab arrays         */
    } 
    *pChamber;                             /* A pointer instance               */

#endif
};

#endif // end FFGRAB_H
