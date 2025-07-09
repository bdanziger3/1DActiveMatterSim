import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
from matplotlib.animation import FuncAnimation
import moviepy

# add src dir to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

def write_videofile_realtime(video_clip:moviepy.VideoFileClip, desired_clip_length:float, output_filename:str=None):
    """
    Rescales the framerate of a video clip to matech the desired total length of the video
    Takes the current nuber of frames to recalculate the correct framerate and reencode

    Saves the video to the specified filepath if`output_filename` is proivded.

    Returns: A new `moviepy.VideoFileClip` object with the desired fps
    """
    new_fps = video_clip.n_frames / desired_clip_length
    new_clip = video_clip.with_fps(new_fps, change_duration=True)

    if output_filename is not None:
        new_clip.write_videofile(output_filename)

    return new_clip


def resize_videofile(video_clip:moviepy.VideoFileClip, desired_clip_length:float, desired_clip_fps:float, output_filename:str=None):
    """
    Rescales the framerate of a video clip to matech the desired total length of the video
    Takes the current nuber of frames to recalculate the correct framerate and reencode

    Saves the video to the specified filepath if`output_filename` is proivded.

    Returns: A new `moviepy.VideoFileClip` object with the desired fps
    """
    new_fps = video_clip.n_frames / desired_clip_length
    new_clip = video_clip.with_fps(new_fps, change_duration=True)
    new_clip_resized = new_clip.with_fps(desired_clip_fps, change_duration=False)

    if output_filename is not None:
        new_clip_resized.write_videofile(output_filename)

    return new_clip_resized