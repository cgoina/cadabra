function [video, audio] = mmread(filename, frames, time, disableVideo, disableAudio, matlabCommand)
% function [video, audio] = mmread(filename, frames, time, disableVideo, disableAudio)
% mmread reads virtually any media file.  If Windows Media Play can play
% it, so should mmread.  It uses the Window's DirectX infrastructure to
% render the media, so other OSs are out of luck.
%
% INPUT
% filename      input file to read (mpg, avi, wmv, asf, wav, mp3, gif, jpg, ...)
% frames        specifies which video frames to capture, default [] for all or
%               to specify time
% time          [startTime stopTime], default [] for all
% disableVideo  disables ALL video capturing, to save memory or time
% disableAudio  disables ALL audio capturing, to save memory or time
% matlabCommand Do not return the video structure, but call the function
%               specified by matlabCommand.  The function definition must
%               match that of processFrame.m.  See processFrame.m for more
%               information.
%
% OUTPUT
% video is a struct with the following fields:
%   width           width of the video frames
%   height          height of the video frames
%   rate            the frame rate of the video, if it can't be determined
%                   it will be 1.
%   nrFramesTotal   the total number of frames in the movie regardless of
%                   how many were captured.  Unfortunately, this can not
%                   always be determined.  If it is negative then it
%                   is an estimate based upon the duration and rate
%                   (normally accurate to within .1%).   It can be 0,
%                   in which case it could not be determined at all.  If it
%                   is a possitive number then it should always be accurate.
%   frames          a struct array with the following fields:
%       cdata       [height X width X 3] uint8 matricies
%       colormap    always empty
%   times           the corresponding time stamps for the frames (in milliseconds)
%
% audio is a struct with the following fields:
%   nrChannels      the number of channels in the audio stream (1 or 2)
%   rate            sampling rate of the audio, ex. 44100.  If it can't be
%                   determined then it will be 1.
%   bits            bit depth of the samples (8 or 16)
%   data            the real data of the whole audio stream.  This can be
%                   played using wavplay.  If time ranges are specified,
%                   the length of the data may not correspond to the total
%                   time.  This normally happens with movies.  The issue is
%                   that the start of the audio stream is generally counted
%                   at the END of the first frame.  So, time is shifted by
%                   1/framerate.
%   nrFramesTotal   Audio comes in packets or frames when captured, the
%                   division of the audio into frames may or may not make
%                   sense.
%   frames          cell array of uint8s.  Probably not of great use.
%   times           the corresponding time stamps for the frames (in milliseconds)
%
% If there is no video or audio stream the corresponding structure will be
% empty.
%
% Specifying frames does not effect audio capturing.  If you want only a
% subsection of the audio use the 3rd parameter "time".  Specifying time
% effects both audio and video.  Time is specified in seconds (subsecond
% resolution is supported with fractional numbers ex. 1.125), starting at 0.
% Time is defined as startTime (inclusive) to stopTime (exclusive), or
% using set notation [startTime stopTime).
%
% If there are multiple video or audio streams, then the structure will be
% of length > 1.  For example: audio(1).data and audio(2).data.
%
% Images work, however the frames must be specified.  For some reason
% DirectShow doesn't ever stop when "playing" an image.  So to deal with
% this, I added support so that the processing stops once the last
% specified frame is captured instead of waiting until the media completes.
%
% EXAMPLES
% [video, audio] = mmread('chimes.wav'); % read whole wav file
% wavplay(audio.data,audio.rate);
%
% video = mmread('mymovie.mpg'); % read whole movie
% movie(video.frames);
%
% video = mmread('mymovie.mpg',1:10); %get only the first 10 frames
%
% video = mmread('mymovie.mpg',[],[0 3.5]); %read the first 3.5 seconds of the video
%
% [video, audio] = mmread('chimes.wav',[],[0 0.25]); %read the first 0.25 seconds of the wav
% [video, audio] = mmread('chimes.wav',[],[0.25 0.5]); %read 0.25 to 0.5 seconds of the wav, there is no overlap with the previous example.
%
% video = mmread('mymovie.mpg',[],[],false,true); %read all frames, disable audio
%
% Written by Micah Richert

if nargin < 6
    matlabCommand = '';
    if nargin < 5
        disableAudio = false;
        if nargin < 4
            disableVideo = false;
            if nargin < 3
                time = [];
                if nargin < 2
                    frames = [];
                end
            end
        end
    end
end

try
    mexDDGrab('buildGraph',filename);
    if (isempty(time))
        mexDDGrab('setFrames',frames);
    else
        if (numel(time) ~= 2)
            error('time must be a vector of length 2: [startTime stopTime]');
        end
        mexDDGrab('setTime',time(1),time(2));
    end
    if (disableVideo)
        mexDDGrab('disableVideo');
    end;
    if (disableAudio | nargout < 2)
        mexDDGrab('disableAudio');
    end;
    mexDDGrab('setMatlabCommand',matlabCommand);
        
    mexDDGrab('doCapture');

    [nrVideoStreams, nrAudioStreams] = mexDDGrab('getCaptureInfo');

    video = struct('width',{},'height',{},'nrFramesTotal',{},'frames',{});
    audio = struct('nrChannels',{},'rate',{},'bits',{},'nrFramesTotal',{},'data',{},'frames',{});

    warned = false;

    % we can only get the video frames if we don't process a matlabCommand
    if strcmp(matlabCommand,'')
        % loop through getting all of the video data from each stream
        for i=1:nrVideoStreams
            [width, height, rate, nrFramesCaptured, nrFramesTotal] = mexDDGrab('getVideoInfo',i-1);
            video(i).width = width;
            video(i).height = height;
            video(i).rate = rate;
            video(i).nrFramesTotal = nrFramesTotal;
            video(i).frames = struct('cdata',cell(1,nrFramesCaptured),'colormap',cell(1,nrFramesCaptured));

            if (nrFramesTotal > 0 && any(frames > nrFramesTotal))
                warning('mmread:general',['Frame(s) ' num2str(frames(frames>nrFramesTotal)) ' exceed the number of frames in the movie.']);
            end

            scanline = ceil(width*3/4)*4; % the scanline size must be a multiple of 4.

            for f=1:nrFramesCaptured
                [data, time] = mexDDGrab('getVideoFrame',i-1,f-1);

                if (numel(data) ~= scanline*height)
                    if (numel(data) > 3*width*height)
                        if (~warned)
                            warning('mmread:general','dimensions do not match data size. Guessing badly...');
                            warned = true;
                        end
                        scanline = width*3;
                        data = data(1:3*width*height);
                    else
                        if (f == 1)
                            error('dimensions do not match data size. Too little data.');
                        else
                            warning(['dimensions do not match data size. Too little data for ' num2str(f) 'th frame.']);
                            continue;
                        end
                    end
                end

                % if there is any extra scanline data, remove it
                data = reshape(data,scanline,height);
                data = data(1:3*width,:);

                % the data ordering is wrong for matlab images, so permute it
                tmp = permute(reshape(data, 3, width, height),[3 2 1]);
                % the images are also upside down and colors were backwards.
                video(i).frames(f).cdata = tmp(end:-1:1,:,3:-1:1);
                video(i).times(f) = time;
            end

            % if frames are specified then make sure that the order is the same
            if (~isempty(frames) && nrFramesCaptured > 0)
                [uniqueFrames, dummy, frameOrder] = unique(frames);
                if (length(uniqueFrames) > nrFramesCaptured)
                    warning('mmread:general','Not all frames specified were captured.  Returning what was captured, but order may be different than specified.');
                    remainingFrames = frames(frames<=uniqueFrames(nrFramesCaptured));
                    [dummy, dummy, frameOrder] = unique(remainingFrames);
                end

                video(i).frames = video(i).frames(frameOrder);
                video(i).times = video(i).times(frameOrder);
            end
        end
    end

    % loop through getting all of the audio data from each stream
    for i=1:nrAudioStreams
        [nrChannels, rate, bits, nrFramesCaptured, nrFramesTotal] = mexDDGrab('getAudioInfo',i-1);
        audio(i).nrChannels = nrChannels;
        audio(i).rate = rate;
        audio(i).bits = bits;
        audio(i).nrFramesTotal = nrFramesTotal;
        audio(i).frames = cell(1,nrFramesCaptured);
        for f=1:nrFramesCaptured
            [data, time] = mexDDGrab('getAudioFrame',i-1,f-1);
            audio(i).frames{f} = data;
            audio(i).times(f) = time;
        end
        % combine the data across frames
        d = double(cat(1,audio(i).frames{:}));
        % reshape and rescale the data so that it is nrChannels x Samples
        % and -1.0 to 1.0.  This should be the same output as wavread.
        audio(i).data = reshape(d/2^(bits-1),nrChannels,length(d)/nrChannels)';
    end

    mexDDGrab('cleanUp');
catch
    err = lasterror;
    mexDDGrab('cleanUp');
    rethrow(err);
end
