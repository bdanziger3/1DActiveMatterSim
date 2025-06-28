import sys, os
from enum import Enum
import numpy as np
import time

# add src dir to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from simulation.sim_structs import SimulationData, SimulationParameters
from utils.print_tools import ProgressBar

DATE_START_INDEX = len("Saved at ")

class DataFileType(Enum):
    SEPERATE_FILES = 1
    SEQUENTIAL_TXT = 2
    ROWWISE_TXT = 3
    JSON = 4

def parse_julia_datetime(datetime_str:str) -> time.struct_time:
    t = time.strptime(datetime_str.strip(), r"%Y-%m-%dT%H:%M:%S")
    return t

def string_list_to_number_list(str_list, seperator=",", out_type=np.float64):
    if type(out_type) == type(np.float64):
        return [out_type(x) for x in str_list.split(seperator)]
    else:
        # different types for each param
        if len(out_type) >= 1:
            out_list = []
            for t, x in zip(out_type, str_list.split(seperator)):
                if t == bool:
                    out_list.append(x.lower() == "true")
                else:
                    out_list.append(t(x))

            return out_list
            # return [t(x) for t, x in zip(type_list, str_list.split(seperator))]



def loadsim(input_fn:str, file_type:DataFileType=DataFileType.ROWWISE_TXT) -> SimulationData:
    # open and read file
    fn = open(input_fn, "r")
    line = fn.readline()

    if file_type == DataFileType.ROWWISE_TXT:

        # Check that data file header is correct        
        assert (("row wise txt" in line.lower()) or ("rowwise txt" in line.lower())), "Rowwisetxt data file does not have correct first line. File may be corrupted."

        line = fn.readline()

        # if the current line is a date, read it and skip to next line. Else do nothing
        sim_dates:list[time.struct_time] = []
        datestartindex = len("Saved at ")
        try:
            sim_dates.append(parse_julia_datetime(line[datestartindex:]))
            line = fn.readline()
        except ValueError as e:
            sim_dates.append(parse_julia_datetime("0001-01-01T00:00:00"))

        assert "simulation parameters" in line.lower(), "ROWWISE_TXT data file does not have correct header structure. File may be corrupted."
       
        # read next line to get number of position data points
        sim_param_info = fn.readline()
        sim_param_info_list = string_list_to_number_list(sim_param_info, out_type=SimulationParameters.expected_types)
        sim_params:SimulationParameters = SimulationParameters.construct_from_list(sim_param_info_list)
        ntimes = sim_params.get_nsaves()
            
        # now read through the data in the order: [[positions], [spins]]
        line = fn.readline()

        # positions
        assert "particle states" in line.lower(), "ROWWISE_TXT  file does not have correct `Particle States` header. File may be corrupted."
        # initialize matrix and read in data from file
        pb = ProgressBar(ntimes, "Loading positions...", "Done")
        pos_matrix = np.zeros([ntimes, sim_params._num_particles], dtype=np.float64)
        spins_matrix = np.zeros([ntimes, sim_params._num_particles], dtype=np.int8)

        for i in range(ntimes):  
            pb.sparse_progress(i)              
            data_line = fn.readline()
            data_list = string_list_to_number_list(data_line)
            pos_matrix[i, :] = data_list[:sim_params._num_particles]        
            spins_matrix[i, :] = data_list[sim_params._num_particles:]        


    
    elif file_type == DataFileType.SEQUENTIAL_TXT:
        # Check that data file header is correct
        if "sequential txt" in line.lower():
            # skip line if first line is header for Sequential Text file
            line = fn.readline()

        assert "simulation parameters" in line.lower(), "SEQUENTIAL_TXT data file does not have correct header structure. File may be corrupted."
       
        # read next line to get number of position data points
        sim_param_info = fn.readline()
        sim_param_info_list = string_list_to_number_list(sim_param_info, out_type=SimulationParameters.expected_types)
        sim_params:SimulationParameters = SimulationParameters.construct_from_list(sim_param_info_list)
        ntimes = sim_params.get_nsaves()
            
        # now read through the data in the order: [positions, spins]
        line = fn.readline()

        # positions
        assert "positions" in line.lower(), "SEQUENTIAL_TXT data file does not have correct `Positions` header. File may be corrupted."
        # initialize matrix and read in data from file
        pb = ProgressBar(ntimes, "Loading positions...", "Done")
        pos_matrix = np.zeros([ntimes, sim_params._num_particles], dtype=np.float64)
    
        for i in range(ntimes):  
            pb.sparse_progress(i)              
            data_line = fn.readline()
            pos_matrix[i,:] = string_list_to_number_list(data_line)        

        # spins
        # initialize matrix and read in data from file
        line = fn.readline()
        assert "spins" in line.lower(), "SEQUENTUAL_TXT data file does not have correct `Spins` header line. File may be corrupted."
        
        pb = ProgressBar(ntimes, "Loading spins...", "Done")
        spins_matrix = np.zeros([ntimes, sim_params._num_particles], dtype=np.int8)
        for i in range(ntimes):
            pb.sparse_progress(i)
            data_line = fn.readline()
            spins_matrix[i,:] = string_list_to_number_list(data_line, out_type=np.int8)
        
    # close file
    fn.close()
        
    # store matrix data as `SimulationData` object and return
    sim_data = SimulationData(sim_params, sim_params.get_times, pos_matrix, spins_matrix)
    return sim_data


# def get_sim_segements(input_fn:str) -> SimulationData:
#     # open and read file
#     fn = open(input_fn, "r")
#     line = fn.readline()

#     # Check that data file header is correct        
#     assert (("row wise txt" in line.lower()) or ("rowwise txt" in line.lower())), "Rowwisetxt data file does not have correct first line. File may be corrupted."

#     line = fn.readline()

#     end_of_file = False

#     sim_dates:list[time.struct_time] = []
#     line_counter = 0

#     while not end_of_file:

#         # if the current line is a date, read it and skip to next line. Else do nothing
#         try:
#             sim_dates.append(parse_julia_datetime(line[DATE_START_INDEX:]))
#             line = fn.readline()
#         except ValueError as e:
#             sim_dates.append(parse_julia_datetime("0001-01-01T00:00:00"))

#         assert "simulation parameters" in line.lower(), "ROWWISE_TXT data file does not have correct header structure. File may be corrupted."
       
#         # read next line to get number of position data points
#         sim_param_info = fn.readline()
#         sim_param_info_list = string_list_to_number_list(sim_param_info, out_type=SimulationParameters.expected_types)
#         sim_params:SimulationParameters = SimulationParameters.construct_from_list(sim_param_info_list)
#         ntimes = sim_params.get_nsaves()
            
#         # now skip the particle state data
#         line = fn.readline()
#         linecounter += 1
#         assert "particle states" in line.lower(), "ROWWISE_TXT  file does not have correct `Particle States` header. File may be corrupted."
        
#         @assert contains(lowercase(line), "particle states") "Rowwisetxt data file does not have correct `Particle States` header. File may be corrupted."
        
#         # now read through the data in the order: [[positions], [spins]]

#         # positions
#         # initialize matrix and read in data from file
#         pb = ProgressBar(ntimes, "Loading positions...", "Done")
#         pos_matrix = np.zeros([ntimes, sim_params._num_particles], dtype=np.float64)
#         spins_matrix = np.zeros([ntimes, sim_params._num_particles], dtype=np.int8)

#         for i in range(ntimes):  
#             pb.sparse_progress(i)              
#             data_line = fn.readline()
#             data_list = string_list_to_number_list(data_line)
#             pos_matrix[i, :] = data_list[:sim_params._num_particles]        
#             spins_matrix[i, :] = data_list[sim_params._num_particles:]        

    # # close file
    # fn.close()
        
    # # store matrix data as `SimulationData` object and return
    # sim_data = SimulationData(sim_params, sim_params.get_times, pos_matrix, spins_matrix)
    # return sim_data


def loadsim_n_lines(input_fn:str, start_line:int, n_lines:int, file_type:DataFileType=DataFileType.ROWWISE_TXT, change_simparams=False) -> SimulationData:
    """
    Reads only the `n_lines` lines of data in a data file starting at line `start_line`.
    `start_line` is indexed at 0

    use `start_line=-1` to read the last `n_lines` lines of the file

    use `change_simparams=True` to change the `SimulationParams` object to match the length of the number of lines
    """
    # open and read file
    fn = open(input_fn, "r")
    line = fn.readline()

    if file_type == DataFileType.ROWWISE_TXT:

        # Check that data file header is correct        
        assert (("row wise txt" in line.lower()) or ("rowwise txt" in line.lower())), "Rowwisetxt data file does not have correct first line. File may be corrupted."

        line = fn.readline()

        # if the current line is a date, read it and skip to next line. Else do nothing
        sim_dates = []
        datestartindex = len("Saved at ")
        try:
            sim_dates.append(parse_julia_datetime(line[datestartindex:]))
            line = fn.readline()
        except ValueError:
            sim_dates.append(parse_julia_datetime("0001-01-01T00:00:00"))


        assert "simulation parameters" in line.lower(), "ROWWISE_TXT data file does not have correct header structure. File may be corrupted."
       
        # read next line to get number of position data points
        sim_param_info = fn.readline()
        sim_param_info_list = string_list_to_number_list(sim_param_info, out_type=SimulationParameters.expected_types)
        print(sim_param_info)
        print(sim_param_info_list)
        sim_params:SimulationParameters = SimulationParameters.construct_from_list(sim_param_info_list)
        ntimes = sim_params.get_nsaves()
            
        line = fn.readline()
        assert "particle states" in line.lower(), "ROWWISE_TXT  file does not have correct `Particle States` header. File may be corrupted."
        
        # now read through the requested lines of data in the order: [[positions], [spins]]
        pb = ProgressBar(ntimes, "Loading positions...", "Done")

        # initialize matrix and read in data from file
        pos_matrix = np.zeros([n_lines, sim_params._num_particles], dtype=np.float64)
        spins_matrix = np.zeros([n_lines, sim_params._num_particles], dtype=np.int8)
        
        # skip to start line
        if start_line == -1:
            lines_to_skip = ntimes - n_lines
        else:
            lines_to_skip = start_line

        for _ in range(lines_to_skip):
            fn.readline()

        # read `n_lines` lines
        for i in range(n_lines):  
            pb.sparse_progress(i)              
            data_line = fn.readline()
            data_list = string_list_to_number_list(data_line)
            pos_matrix[i, :] = data_list[:sim_params._num_particles]        
            spins_matrix[i, :] = data_list[sim_params._num_particles:]        

    # close file
    fn.close()

    # Change simparams if toggled
    if change_simparams:
        sim_params._total_time = sim_params._snapshot_dt * n_lines
        
    # store matrix data as `SimulationData` object and return
    sim_data = SimulationData(sim_params, sim_params.get_times, pos_matrix, spins_matrix)
    return sim_data






    
def write_dlm(file_name:str, data:list, delimiter:str):
    f = open(file_name, "w")

    array_dim = np.size(data)

    for r in range(array_dim[0]):
        for c in range(array_dim[1]):
            f.write(",".join([str(x) for x in data[r]])+"\n")
    f.close()

def save_sim(sim_data:SimulationData, output_fn:str, file_type:DataFileType=DataFileType.SEQUENTIAL_TXT, clear_file:bool=False):
    """
    Saves simulation data from a `SimulationData` object to a data file.

    Can specify the output file type with `filetype`.

    By default appends the data to the file at `outputfilename`.
    Set `clearfile` to `true` to clear the file before writing instead.

    """
    file_str = f"{output_fn}_simdata.txt"
    fn = open(file_str, "w" if clear_file else "a")
    if file_type == DataFileType.SEQUENTIAL_TXT:
        fn.write("Simulation Parameters\n")
        fn.write(",".join([str(x) for x in sim_data._sim_params.as_array()])+"\n")
        fn.write("Positions\n")
        fn.write(",".join([str(x) for x in sim_data.positions])+"\n")
        fn.write("Spins\n")
        fn.write(",".join([str(x) for x in sim_data.spins])+"\n")
    elif file_type == DataFileType.SEPERATE_FILES:
        sim_datafile_prefix = output_fn
        fn.write(",".join([str(x) for x in sim_data.spins])+"\n")
        write_dlm(f"{sim_datafile_prefix}_x.txt", sim_data.positions, ",")
        write_dlm(f"{sim_datafile_prefix}_S.txt", sim_data.spins, ",")
        write_dlm(f"{sim_datafile_prefix}_t.txt", sim_data.times, ",")
    
    fn.close()
