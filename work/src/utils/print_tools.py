# Blake Danziger
# MSc Theoretical Physics
# Dissertaion
# 
# Library of functions useful for printing when running processes in Python

import numpy as np
import time

def format_time(time_seconds):
    if time_seconds < 60:
        return f"{int(np.round(time_seconds))} seconds"
    elif time_seconds < 3600:
        time_minutes = time_seconds // 60
        return f"{int(np.round(time_minutes))} minutes {int(np.round(time_seconds - (60 * time_minutes)))} seconds"
    else:
        time_hours = time_seconds // 3600
        time_minutes = (time_seconds  - (3600 * time_hours)) // 60
        return f"{int(np.round(time_hours))} hours {int(np.round(time_minutes))} minutes {int(np.round(time_seconds - (3600 * time_hours) - (60 * time_minutes)))} seconds"

class ProgressBar:
    def __init__(self, total_steps:int, process_msg:str, finished_msg:str):
        self.total_steps = total_steps
        self.process_msg = process_msg
        self.finsihed_msg = finished_msg

        self.current_step = 0
        self.clock = time.time()
        self.step_time = None
        self.last_time = None

        self.sparse_measured = False

    def step(self):
        self.current_step += 1

    def reset(self):
        self.current_step = 0

    def measure_time(self):
        new_time = time.time()
        self.step_time = new_time - self.clock
        self.clock = new_time

    def sparse_progress(self, current_step:int, update_time_s:float = 0.5):
        # initialise sparse progress
        if current_step < 10:
            self.measure_time()

        # measure first step
        elif not self.sparse_measured:
            if time.time() - self.clock >= update_time_s:
                self.measure_time()
                self.sparse_measured = True

        # measure time to do previous step and print
        else:
            self.sparse_measured = False
            self.measure_time()

            estimate_seconds = self.step_time * (self.total_steps - current_step)
            time_msg = f". Approx. time remaining: {format_time(estimate_seconds)}"

            print(f"\033[K{self.process_msg}: {np.round(100 * current_step / self.total_steps, 2)} %{time_msg}", end="\r")
            if self.total_steps == current_step: 
                print(f"\033[K{self.finsihed_msg}")
        


        # print the results


    def print_progress(self, current_step, print_time=False):
        self.measure_time()
        if print_time:
            estimate_seconds = self.step_time * (self.total_steps - current_step)
            
            time_msg = f". Approx. time remaining: {format_time(estimate_seconds)}"
        else:
            time_msg = ""

        pct_steps = self.total_steps // 10
        if current_step % pct_steps == 0:
            print(f"\033[K{self.process_msg}: {np.round(100 * current_step / self.total_steps, 2)} %{time_msg}", end="\r")
        if self.total_steps == current_step: 
            print(f"\033[K{self.finished_msg}")




def print_progress(current_step: int, total_steps: int, process_msg:str, finished_msg:str):
    pct_steps = total_steps // 100
    if current_step % pct_steps == 0:
        print(f"\033[K{process_msg}: {np.round(100 * current_step / total_steps, 2)} %", end="\r")
    if total_steps == current_step: 
        print(f"\033[K{finished_msg}")