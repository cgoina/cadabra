#include "avbin.h"

#define Logprintf(...) printf(__VA_ARGS__)

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
    // Close the video file
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

void avbin_dump(AVbinFile* file, const char* filename)
{
    av_dump_format(file->format_context, 0, filename, 0);
}

AVbinResult avbin_stream_info(AVbinFile* file, int32_t stream_index, AVbinStreamInfo* info)
{
    AVCodec* pCodec = avcodec_find_decoder(file->format_context->streams[stream_index]->codecpar->codec_id);
    if (pCodec == NULL) {
       Logprintf("No codec found for %d\n", file->format_context->streams[stream_index]->codecpar->codec_id);
    	// unsupported codec
    	return AVBIN_RESULT_ERROR;
    }
    Logprintf("Found codec %s for %d\n", pCodec->name, file->format_context->streams[stream_index]->codecpar->codec_id);
    AVCodecParameters* codec_params = file->format_context->streams[stream_index]->codecpar;

    /* Error if not large enough for version 1 */
    if (info->structure_size < sizeof *info)
        return AVBIN_RESULT_ERROR;

    switch (codec_params->codec_type)
    {
        case AVMEDIA_TYPE_VIDEO:
            info->type = AVBIN_STREAM_TYPE_VIDEO;
            info->video.width = codec_params->width;
            info->video.height = codec_params->height;
            info->video.nb_frames = file->format_context->streams[stream_index]->nb_frames;
            AVRational sample_aspect_ratio; 
            if ((codec_params->sample_aspect_ratio.num != 0 && codec_params->sample_aspect_ratio.den != 0) ||
                codec_params->width == 0 || codec_params->height == 0)
            {
                sample_aspect_ratio = codec_params->sample_aspect_ratio;
            }
            else
            {
                sample_aspect_ratio = av_d2q((double)codec_params->height / (double)codec_params->width, 20);
            }
            info->video.sample_aspect_num = sample_aspect_ratio.num;
            info->video.sample_aspect_den = sample_aspect_ratio.den;
            break;
        case AVMEDIA_TYPE_AUDIO:
            info->type = AVBIN_STREAM_TYPE_AUDIO;
            info->audio.sample_rate = codec_params->sample_rate;
            info->audio.channels = codec_params->channels;
            switch (codec_params->format)
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
        Logprintf("No codec found for: %d\n", file->format_context->streams[stream_index]->codecpar->codec_id);
        // codec not found
        return NULL;
    }
    Logprintf("Found codec %s\n", codec->name);

    /**
     * Note that we must not use the AVCodecContext from the video stream
     * directly! So we have to use avcodec_copy_context() to copy the
     * context to a new location (after allocating memory for it, of
     * course).
     */

    // Copy context
    AVCodecContext* codec_context = avcodec_alloc_context3(codec); // [7]
    int ret = avcodec_parameters_to_context(codec_context, file->format_context->streams[stream_index]->codecpar);
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
    stream->type = codec_context->codec_type;
    stream->frame = frame;
    // initialize SWS context for software scaling
    stream->sws_ctx = sws_getContext(
        codec_context->width,
        codec_context->height,
        codec_context->pix_fmt,
        codec_context->width,
        codec_context->height,
        AV_PIX_FMT_RGB24,   // sws_scale destination color scheme
        SWS_BILINEAR,
        NULL,
        NULL,
        NULL
    );

    return stream;
}

void avbin_close_stream(AVbinStream *stream)
{
    if (stream->frame)
    {
        // Free the frame
        av_frame_free(&stream->frame);
        av_free(stream->frame);
    }

    // Close the codec
    avcodec_close(stream->codec_context);

    free(stream);
}

int32_t avbin_read_next_packet(AVbinFile* file, AVbinPacket* packet)
{
    if (packet == NULL)
        return AVBIN_RESULT_ERROR;
    if (av_read_frame(file->format_context, packet->packet) < 0)
        return AVBIN_RESULT_ERROR;

    packet->timestamp = av_rescale_q(
        packet->packet->dts,
        file->format_context->streams[packet->packet->stream_index]->time_base,
        AV_TIME_BASE_Q
    );
    packet->stream_index = packet->packet->stream_index;
    return AVBIN_RESULT_OK;
}

int32_t avbin_decode_audio(AVbinStream* stream, AVbinPacket* packet)
{
    int bytes_used = 0;
    if (stream->type != AVMEDIA_TYPE_AUDIO)
        return AVBIN_RESULT_ERROR;

    // FIXME !!!!!
    return bytes_used;
}

int32_t avbin_decode_video_frame(AVbinStream *stream, AVbinPacket* packet, 
                                 uint8_t* output_buffer, int output_size)
{
    AVFrame* pFrameRGB = av_frame_alloc();
    if (pFrameRGB == NULL)
    {
        // Could not allocate frame
        Logprintf("avbin_decode_video_frame: error allocating the RGB frame\n");
        return -1;
    }
    int nbytes = av_image_fill_arrays(
        pFrameRGB->data,
        pFrameRGB->linesize,
        output_buffer,
        AV_PIX_FMT_RGB24,
        stream->codec_context->width,
        stream->codec_context->height,
        32
    );
    Logprintf("avbin_decode_video_frame: av_image_fill_arrays linesize=%d, width=%d, height=%d nbytes=%d, output size=%d\n",
              pFrameRGB->linesize[0], stream->codec_context->width, stream->codec_context->height, nbytes, output_size);
    if (nbytes < 0)
    {
        return -1;
    }
	int32_t ret = avcodec_send_packet(stream->codec_context, packet->packet);
    if (ret < 0)
    {
        // could not send packet for decoding
        Logprintf("avcodec_send_packet could not send packet for decoding: %d\n", ret);
        av_frame_free(&pFrameRGB);
        av_free(pFrameRGB);
        av_packet_unref(packet->packet);
        return -1;
    } else {
        Logprintf("avcodec_send_packet OK!\n");       
    }
    while (ret >= 0)
    {
        ret = avcodec_receive_frame(stream->codec_context, stream->frame);
        Logprintf("avcodec_receive_frame result=%d!\n", ret);
        if (ret == AVERROR(EAGAIN) || ret == AVERROR_EOF)
        {
            // EOF
            break;
        }
        else if (ret < 0)
        {
            // other error
            break;
        } else {
            // successful decoding => convert the image from its native format to RGB
            Logprintf("avcodec_receive_frame successfully decoded packet result=%d!\n", ret);
            sws_scale(stream->sws_ctx,
                      (uint8_t const * const *)stream->frame->data,
                      stream->frame->linesize,
                      0,
                      stream->codec_context->height,
                      pFrameRGB->data,
                      pFrameRGB->linesize
            );

            // align the line size to be <frame_width> * <bytes_per_pixels> 
            for (int y = 0; y < stream->frame->height; y++)
            {
                memcpy(output_buffer + y * stream->frame->width * 3, // in this case we have 3 bytes per pixel
                       output_buffer + y * pFrameRGB->linesize[0],
                       stream->frame->width * 3);
            }
            Logprintf("av_image_copy_to_buffer -> %d!\n", ret);
            break;
        }
    }
    av_frame_free(&pFrameRGB);
    av_free(pFrameRGB);
    av_packet_unref(packet->packet);
    return ret;
}
