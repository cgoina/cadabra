#ifndef AVBIN_H
#define AVBIN_H

#ifdef __cplusplus
extern "C" {
	#include <libavcodec/avcodec.h>
	#include <libavformat/avformat.h>
	#include <libavutil/avutil.h>
	#include <libavutil/dict.h>
	#include <libavutil/imgutils.h>
	#include <libavutil/mathematics.h>
	#include <libswscale/swscale.h>
}
#endif

/**
 * Point in time, or a time range; given in microseconds.
 */
typedef int64_t AVbinTimestamp;

/**
 * Error-checked function result.
 */
typedef enum _AVbinResult {
    AVBIN_RESULT_ERROR = -1,
    AVBIN_RESULT_OK = 0
} AVbinResult;

/**
 * Type of a stream; currently only video and audio streams are supported.
 */
typedef enum _AVbinStreamType {
    AVBIN_STREAM_TYPE_UNKNOWN = 0,
    AVBIN_STREAM_TYPE_VIDEO = 1,
    AVBIN_STREAM_TYPE_AUDIO = 2
} AVbinStreamType;

/**
 * The sample format for audio data.
 */
typedef enum _AVbinSampleFormat {
    AVBIN_SAMPLE_FORMAT_UNKNOWN = -1,
    /** Unsigned byte */
    AVBIN_SAMPLE_FORMAT_U8 = 0,
    /** Signed 16-bit integer */
    AVBIN_SAMPLE_FORMAT_S16 = 1,
    /** Signed 24-bit integer
     *  AVBIN_SAMPLE_FORMAT_S24 removed upstream.  Removed here in AVbin 11
     */
    AVBIN_SAMPLE_FORMAT_UNUSED = 2,
    /** Signed 32-bit integer */
    AVBIN_SAMPLE_FORMAT_S32 = 3,
    /** 32-bit IEEE floating-point */
    AVBIN_SAMPLE_FORMAT_FLOAT = 4
} AVbinSampleFormat;

/**
 * Initialization Options
 */
typedef struct _AVbinOptions {
    /**
     * Size of this structure, in bytes.  This must be filled in by the
     * application before passing to AVbin.
     */
    size_t structure_size;

    /**
     * Number of threads to attempt to use.  Using the recommended thread_count of
     * 0 means try to detect the number of CPU cores and set threads to
     * (num cores + 1).  A thread_count of 1 or a negative number means single
     * threaded.  Any other number will result in an attempt to set that many threads.
     */
    int32_t thread_count;
} AVbinOptions;

/**
 * File details.  The info struct is filled in by avbin_get_file_info.
 */
typedef struct _AVbinFileInfo {
    /**
     * Size of this structure, in bytes.  This must be filled in by the
     * application before passing to AVbin.
     */
    size_t structure_size;

    /**
     * Number of streams contained in the file.
     */
    int32_t n_streams;

    /**
     * Starting time of all streams.
     */
    AVbinTimestamp start_time;

    /**
     * Duration of the file.  Does not include the time given in start_time.
     */
    AVbinTimestamp duration;

    /**
     * @name Metadata fields
     *
     * File metadata.
     *
     * Strings are NUL-terminated and may be omitted (the first character
     * NUL) if the file does not contain appropriate information.  The
     * encoding of the strings is unspecified.
     */
    /*@{*/
    char title[512];
    char author[512];
    char copyright[512];
    char comment[512];
    char album[512];
    int32_t year;
    int32_t track;
    char genre[32];
    /*@}*/
} AVbinFileInfo;

typedef struct _AVbinFile {
    AVFormatContext *format_context;
} AVbinFile;

/**
 * Stream details.
 *
 * A stream is a single audio track or video.  Most audio files contain one
 * audio stream.  Most video files contain one audio stream and one video
 * stream.  More than one audio stream may indicate the presence of multiple
 * languages which can be selected (however at this time AVbin does not
 * provide language information).
 */
typedef struct _AVbinStreamInfo {
    /**
     * Size of this structure, in bytes.  This must be filled in by the
     * application before passing to AVbin.
     */
    size_t structure_size;

    /**
     * The type of stream; either audio or video.
     */
    AVbinStreamType type;

    union {
        struct {
            /**
             * Width of the video image, in pixels.  This is the width
             * of actual video data, and is not necessarily the size the
             * video is to be displayed at (see sample_aspect_num).
             */
            uint32_t width;

            /**
             * Height of the video image, in pixels.
             */
            uint32_t height;

            /**
             * Aspect-ratio of each pixel.  The aspect is given by dividing
             * sample_aspect_num by sample_aspect_den.
             */
            uint32_t sample_aspect_num;
            uint32_t sample_aspect_den;

            /** Frame rate, in frames per second.  The frame rate is given by
             * dividing frame_rate_num by frame_rate_den.
             *
             * @version Version 8.  Requires frame_rate feature.
             *
             * REMOVED IN VERSION 11 - see note on avbin_have_feature()
             */
            uint32_t frame_rate_num;
            uint32_t frame_rate_den;
        } video;

        struct {
            /**
             * Data type of audio samples.
             */
            AVbinSampleFormat sample_format;

            /**
             * Number of samples per second, in Hz.
             */
            uint32_t sample_rate;

            /**
             * Number of bits per sample; typically 8 or 16.
             */
            uint32_t sample_bits;

            /**
             * Number of interleaved audio channels.  Typically 1 for
             * monoaural, 2 for stereo.  Higher channel numbers are used for
             * surround sound, however AVbin does not currently provide a way
             * to access the arrangement of these channels.
             */
            uint32_t channels;
        } audio;
    };
} AVbinStreamInfo;

/**
 * A single packet of stream data.
 *
 * The structure size must be initialised before passing to avbin_read.  The
 * data will point to a block of memory allocated by AVbin -- you must not
 * free it.  The data will be valid until the next time you call avbin_read,
 * or until the file is closed.
 */
typedef struct _AVbinPacket {
    /**
     * Size of this structure, in bytes.  This must be filled in by the
     * application before passing to AVbin.
     */
    size_t structure_size;

    /**
     * The time at which this packet is to be played.  This can be used
     * to synchronise audio and video data.
     */
    AVbinTimestamp timestamp;

    /**
     * The stream this packet contains data for.
     */
    uint32_t stream_index;

    AVPacket *packet;
    uint8_t *data;
    size_t size;
} AVbinPacket;

typedef struct _AVbinStream {
    int type;
    AVFormatContext *format_context;
    AVCodecContext *orig_codec_context;
    AVCodecContext *codec_context;
    AVFrame *frame;
} AVbinStream;


AVbinResult avbin_init();
AVbinResult avbin_init_options(AVbinOptions * options_ptr);
void avbin_set_log_level(int level);
AVbinFile* avbin_open_filename(const char* filename);
AVbinFile* avbin_open_filename_with_format(const char* filename, char* format);
void avbin_close_file(AVbinFile* file);
AVbinResult avbin_file_info(AVbinFile *file, AVbinFileInfo *info);
AVbinResult avbin_stream_info(AVbinFile* file, int32_t stream_index, AVbinStreamInfo* info);
AVbinStream* avbin_open_stream(AVbinFile *file, int32_t stream_index);
void avbin_close_stream(AVbinStream *stream);
int32_t avbin_read_next_packet(AVbinFile* file, AVbinPacket* packet);
int32_t avbin_decode_audio(AVbinStream* stream, AVbinPacket* packet);
int32_t avbin_decode_video(AVbinStream* stream, AVbinPacket* packet);

#endif // end AVBIN_H
