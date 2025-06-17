import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
from matplotlib.animation import FuncAnimation
import moviepy

# add src dir to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

"""
Rescales the framerate of a video clip to matech the desired total length of the video
Takes the current nuber of frames to recalculate the correct framerate and reencode

Returns: A new `moviepy.VideoFileClip` object with the desired fps
"""
# def rescale_framerate(video_clip, desired_clip_length):
#     clip = moviepy.VideoFileClip("/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/data/16-6/N10000-alignsimple-t0.5_smalldots.mp4")


old_filepath = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/data/16-6/N10000-alignsimple-t0.5_smalldots.mp4"
new_filepath = "/Users/blakedanziger/Documents/Grad/MSc Theoretical Physics/Dissertation/Dev/work/data/16-6/N10000-alignsimple-t0.5_smalldots_realtime.mp4"
clip = moviepy.VideoFileClip(old_filepath)

new_duration = 0.5

new_clip = clip.with_duration(new_duration)
new_fps = clip.n_frames / new_duration
new_clip_1 = clip.with_fps(new_fps, change_duration=True)


print(clip.n_frames, new_clip.n_frames, new_clip_1.n_frames)
print(clip.fps, new_clip.fps, new_clip_1.fps)
print(clip.duration, new_clip.duration, new_clip_1.duration)

new_clip_1.write_videofile(new_filepath)