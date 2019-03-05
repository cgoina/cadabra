#include "avbin.h"

static int32_t avbin_thread_count = 1;

AVbinResult avbin_init()
{
    return avbin_init_options(NULL);
}

AVbinResult avbin_init_options(AVbinOptions* options_ptr)
{
    if (options_ptr == NULL)
    {
        options_ptr = (AVbinOptions*)malloc(sizeof(options_ptr));
        if (options_ptr == NULL) return AVBIN_RESULT_ERROR;

        // Set defaults...
        options_ptr->structure_size = sizeof(AVbinOptions);
        options_ptr->thread_count = 1; // single thread
    }

    // What version did we get?
    AVbinOptions * options = NULL;
    if (options_ptr->structure_size == sizeof(AVbinOptions))
    {
        options = options_ptr;
    } else {
        return AVBIN_RESULT_ERROR;
    }

    // Stupid choices deserve single-threading
    if (options->thread_count < 0) options->thread_count = 1;

    avbin_thread_count = options->thread_count;

    av_register_all(); // register all compiled muxers, demuxers and protocols.
    avcodec_register_all(); // register all codecs

    return AVBIN_RESULT_OK;
}

void avbin_set_log_level(int level)
{
    av_log_set_level(level);
}

AVbinFile* avbin_open_filename(const char* filename)
{
  return avbin_open_filename_with_format(filename, NULL);
}

AVbinFile* avbin_open_filename_with_format(const char* filename, char* format)
{
    AVbinFile* file = (AVbinFile*)malloc(sizeof *file);
    AVInputFormat *avformat = NULL;
    if (format) avformat = av_find_input_format(format);

    file->format_context = NULL;    // Zero-initialize
    if (avformat_open_input(&file->format_context, filename, avformat, NULL) != 0)
      goto error;

    if (avformat_find_stream_info(file->format_context, NULL) < 0)
      goto error;

    return file;

error:
    free(file);
    return NULL;
}

void avbin_close_file(AVbinFile* file)
{
    avformat_close_input(&file->format_context);
    free(file);
}

AVbinResult avbin_file_info(AVbinFile *file, AVbinFileInfo *info)
{
    if (info->structure_size < sizeof *info)
        return AVBIN_RESULT_ERROR;

    info->n_streams = file->format_context->nb_streams;
    info->start_time = file->format_context->start_time;
    info->duration = file->format_context->duration;

    // Zero-initialize fields first
    memset(info->title, 0, sizeof(info->title));
    memset(info->author, 0, sizeof(info->author));
    memset(info->copyright, 0, sizeof(info->copyright));
    memset(info->comment, 0, sizeof(info->comment));
    memset(info->album, 0, sizeof(info->album));
    memset(info->genre, 0, sizeof(info->genre));
    info->year = 0;
    info->track = 0;

    AVDictionaryEntry* entry;
    if ((entry = av_dict_get(file->format_context->metadata, "title", NULL, 0)) != NULL)  {
      strncpy(info->title, entry->value, sizeof(info->title));
    }

    if (((entry = av_dict_get(file->format_context->metadata, "artist", NULL, 0)) != NULL) ||
       (entry = av_dict_get(file->format_context->metadata, "album_artist", NULL, 0)) != NULL) {
      strncpy(info->author, entry->value, sizeof(info->author));
    }
    if ((entry = av_dict_get(file->format_context->metadata, "copyright", NULL, 0)) != NULL)  {
      strncpy(info->copyright, entry->value, sizeof(info->copyright));
    }
    if ((entry = av_dict_get(file->format_context->metadata, "comment", NULL, 0)) != NULL)  {
      strncpy(info->comment, entry->value, sizeof(info->comment));
    }
    if ((entry = av_dict_get(file->format_context->metadata, "album", NULL, 0)) != NULL)  {
      strncpy(info->album, entry->value, sizeof(info->album));
    }
    if ((entry = av_dict_get(file->format_context->metadata, "date", NULL, 0)) != NULL)  {
      info->year = atoi(entry->value);
    }
    if ((entry = av_dict_get(file->format_context->metadata, "track", NULL, 0)) != NULL)  {
      info->track = atoi(entry->value);
    }
    if ((entry = av_dict_get(file->format_context->metadata, "genre", NULL, 0)) != NULL)  {
      strncpy(info->genre, entry->value, sizeof(info->genre));
    }

    return AVBIN_RESULT_OK;
}

AVbinResult avbin_stream_info(AVbinFile* file, int32_t stream_index, AVbinStreamInfo* info)
{
    AVCodec* pCodec = avcodec_find_decoder(file->format_context->streams[stream_index]->codecpar->codec_id);
    if (pCodec == NULL) {
    	// unsupported codec
    	return AVBIN_RESULT_ERROR;
    }
    AVCodecContext* context = avcodec_alloc_context3(pCodec); // [7]

    /* Error if not large enough for version 1 */
    if (info->structure_size < sizeof *info)
        return AVBIN_RESULT_ERROR;

    switch (context->codec_type)
    {
        case AVMEDIA_TYPE_VIDEO:
            info->type = AVBIN_STREAM_TYPE_VIDEO;
            info->video.width = context->width;
            info->video.height = context->height;
            info->video.sample_aspect_num = context->sample_aspect_ratio.num;
            info->video.sample_aspect_den = context->sample_aspect_ratio.den;
            break;
        case AVMEDIA_TYPE_AUDIO:
            info->type = AVBIN_STREAM_TYPE_AUDIO;
            info->audio.sample_rate = context->sample_rate;
            info->audio.channels = context->channels;
            switch (context->sample_fmt)
            {
                case AV_SAMPLE_FMT_U8:
                    info->audio.sample_format = AVBIN_SAMPLE_FORMAT_U8;
                    info->audio.sample_bits = 8;
                    break;
                case AV_SAMPLE_FMT_S16:
                    info->audio.sample_format = AVBIN_SAMPLE_FORMAT_S16;
                    info->audio.sample_bits = 16;
                    break;
                case AV_SAMPLE_FMT_S32:
                    info->audio.sample_format = AVBIN_SAMPLE_FORMAT_S32;
                    info->audio.sample_bits = 32;
                    break;
                case AV_SAMPLE_FMT_FLT:
                    info->audio.sample_format = AVBIN_SAMPLE_FORMAT_FLOAT;
                    info->audio.sample_bits = 32;
                    break;
                default:
                  // Unknown sample format
                  info->audio.sample_format = AVBIN_SAMPLE_FORMAT_UNKNOWN;
                  info->audio.sample_bits = (uint32_t)0;
                  break;

                // TODO: support planar formats
            }
            break;

        default:
            info->type = AVBIN_STREAM_TYPE_UNKNOWN;
            break;
    }

 	avcodec_free_context(&context);

    return AVBIN_RESULT_OK;
}

AVbinStream *avbin_open_stream(AVbinFile *file, int32_t stream_index)
{
    if (stream_index < 0 || stream_index >= file->format_context->nb_streams)
        return NULL;

    // Find the decoder for the stream
    AVCodec* codec = avcodec_find_decoder(file->format_context->streams[stream_index]->codecpar->codec_id);
    if (codec == NULL)
    {
        // codec not found
        return NULL;
    }

    /* The Libav api example does this (see libav/libavcodec-api-example.c).
     * The only explanation is "we do not send complete frames".  I tried
     * adding it, and there seemed to be no effect either way.  I'm going to
     * leave it here commented out just in case we find the need to enable it
     * in the future.
     */
    AVCodecContext* orig_codec_context = avcodec_alloc_context3(codec);
    int ret = avcodec_parameters_to_context(orig_codec_context, file->format_context->streams[stream_index]->codecpar);

    /**
     * Note that we must not use the AVCodecContext from the video stream
     * directly! So we have to use avcodec_copy_context() to copy the
     * context to a new location (after allocating memory for it, of
     * course).
     */

    // Copy context
    AVCodecContext* codec_context = avcodec_alloc_context3(codec); // [7]
    ret = avcodec_parameters_to_context(codec_context, file->format_context->streams[stream_index]->codecpar);
    if (ret != 0)
    {
        // error copying context
         return NULL;
    }
     
    if (avcodec_open2(codec_context, codec, NULL) < 0)
    {
        // could not open codec
        return NULL;
    }

    if (avbin_thread_count != 1)
        codec_context->thread_count = avbin_thread_count;

    AVFrame* frame = av_frame_alloc();
    if (frame == NULL)
    {
      printf("Error allocating the frame");
      return NULL;
    }

    AVbinStream* stream = (AVbinStream*)malloc(sizeof *stream);
    if (stream == NULL) {
      printf("Error allocating the stream");
      return NULL;
    }
    stream->format_context = file->format_context;
    stream->codec_context = codec_context;
    stream->orig_codec_context = orig_codec_context;
    stream->type = codec_context->codec_type;
    stream->frame = frame;
    return stream;
}

void avbin_close_stream(AVbinStream *stream)
{
    if (stream->frame)
        av_frame_free(&stream->frame);
    avcodec_close(stream->codec_context);
    free(stream);
}

int32_t avbin_read_next_packet(AVbinFile* file, AVbinPacket* packet)
{
    if (packet == NULL)
        return AVBIN_RESULT_ERROR;
    if (av_read_frame(file->format_context, packet->packet) < 0)
        return AVBIN_RESULT_ERROR;

    packet->timestamp = av_rescale_q(packet->packet->dts,
        AV_TIME_BASE_Q,
        file->format_context->streams[packet->packet->stream_index]->time_base);
    packet->stream_index = packet->packet->stream_index;
    packet->data = packet->packet->data;
    packet->size = packet->packet->size;

    return AVBIN_RESULT_OK;
}

int avbin_decode_audio(AVbinStream* stream)
{
    int bytes_used = 0;
    if (stream->type != AVMEDIA_TYPE_AUDIO)
        return AVBIN_RESULT_ERROR;

    // FIXME !!!!!
    return bytes_used;
}
