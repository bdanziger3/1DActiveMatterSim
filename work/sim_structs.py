import numpy as np
import json

"""
Simulation parameters on which to run particles
"""
class SimulationParameters:
    def __init__(self, num_particles:int, total_time:float, dt:float, v0:float, flip_rate:float, box_width:float=1, interaction:str="none", start_time:float=0):
        self._num_particles:int = int(num_particles)
        self._total_time = total_time
        self._dt = dt
        self._v0 = v0
        self._flip_rate = flip_rate
        self._box_width = box_width
        self.interaction = interaction
        self._start_time = start_time
        # self._starting_positions = starting_positions
        # self._random_starts = random_starts
        # self._random_spins = random_spins
        # self._spin_interaction = spin_interaction


    @classmethod
    def construct_from_dict(params_dict:dict):
        # Construct from Dictionary
        return SimulationParameters(params_dict["numparticles"], params_dict["totaltime"], params_dict["dt"], params_dict["v0"], params_dict["fliprate"], params_dict["boxwidth"], params_dict["starttime"])
    
    @staticmethod
    def construct_from_list(params_list:list):
        return SimulationParameters(*params_list)
    
    def get_times(self) -> np.ndarray[float]:
        end_time = self._start_time + self._total_time
        return np.arange(self._start_time, end_time, self._dt)

    def csv_serialize(self) -> str:
        return f"{self._num_particles},{self._total_time},{self._dt},{self._v0},{self._flip_rate},{self._box_width},{self._start_time}"

    def json_serialize(self) -> str:
        # paramsdict::Dict = Dict("numparticles" => simparams.numparticles, "totaltime" => simparams.totaltime, "dt" => simparams.dt, "v0" => simparams.v0, "fliprate" => simparams.fliprate, "boxwidth" => simparams.boxwidth, "starttime" => simparams.starttime)
        # json_str = JSON.json(asdict(ssimparams))
        # return json_str
        # return JSON.json(asdict(simparams), 2)
        pass

    def get_ntimes(self) -> int:
        return int(np.floor(self._total_time / self._dt) + 1)
    
    def as_array(self) -> list:
        return [self._num_particles, self._total_time, self._dt, self._v0, self._flip_rate, self._box_width, self._start_time]

    def __repr__(self):
        return "num_particles: {}\ntotal_time: {}\ndt: {}\nv0: {}\nflip_rate: {}\nbox_width: {}\nstart_time: {}".format(*self.as_array()) 


class SimulationData():
    def __init__(self, sim_params, times, positions, wrapped_positions, spins):
        self._sim_params: SimulationParameters = sim_params
        self.times = times
        self.positions = positions
        self.wrapped_positions = wrapped_positions
        self.spins = spins