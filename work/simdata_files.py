from enum import Enum
from sim_structs import SimulationData, SimulationParameters
import numpy as np
from utils.print_tools import ProgressBar

class DataFileType(Enum):
    SEPERATE_FILES = 1
    SEQUENTIAL_TEXT = 2
    JSON = 3
    SEQUENTIAL_TEXT_ABS  = 4

def string_list_to_number_list(str_list, seperator=",", out_type=np.float64):
    if type(out_type) == type(np.float64):
        return [out_type(x) for x in str_list.split(seperator)]
    else:
        # different types for each param
        if len(out_type) >= 1:
            type_list = out_type
            return [t(x) for t, x in zip(type_list, str_list.split(seperator))]



def loadsim(input_fn:str, file_type:DataFileType=DataFileType.SEQUENTIAL_TEXT) -> SimulationData:
    if file_type == DataFileType.SEQUENTIAL_TEXT or file_type == DataFileType.SEQUENTIAL_TEXT_ABS:
        # read file
        fn = open(input_fn, "r")

        line = fn.readline()
        # Check that data file header is correct
        assert "simulation parameters" in line.lower(), "SEQUENTIAL_TEXT data file does not have correct first line. File may be corrupted."
       
        # read next line to get number of position data points
        sim_param_info = fn.readline()
        expected_types = [np.float64, np.float64, np.float64, np.float64, np.float64, np.float64, str, np.float64]
        sim_param_info_list = string_list_to_number_list(sim_param_info, out_type=expected_types)
        # sim_param_info_list = [float(x) for x in sim_param_info.split(",")]
        sim_params:SimulationParameters = SimulationParameters.construct_from_list(sim_param_info_list)
        print(sim_params)
        print(sim_params._num_particles)
        ntimes = sim_params.get_ntimes()
            
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
        print(i)

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
    
def write_dlm(file_name:str, data:list, delimiter:str):
    f = open(file_name, "w")

    array_dim = np.size(data)

    for r in range(array_dim[0]):
        for c in range(array_dim[1]):
            f.write(",".join([str(x) for x in data[r]])+"\n")
    f.close()

"""
Saves simulation data from a `SimulationData` object to a data file.

Can specify the output file type with `filetype`.

By default appends the data to the file at `outputfilename`.
Set `clearfile` to `true` to clear the file before writing instead.

"""
def save_sim(sim_data:SimulationData, output_fn:str, file_type:DataFileType=DataFileType.SEQUENTIAL_TEXT, clear_file:bool=False):
    file_str = f"{output_fn}_simdata.txt"
    fn = open(file_str, "w" if clear_file else "a")
    if file_type == DataFileType.SEQUENTIAL_TEXT:
        fn.write("Simulation Parameters\n")
        fn.write(",".join([str(x) for x in sim_data._sim_params.as_array()])+"\n")
        fn.write("Positions\n")
        fn.write(",".join([str(x) for x in sim_data.positions])+"\n")
        fn.write("Wrapped Positions\n")
        fn.write(",".join([str(x) for x in sim_data.wrapped_positions])+"\n")
        fn.write("Spins\n")
        fn.write(",".join([str(x) for x in sim_data.spins])+"\n")
    elif file_type == DataFileType.SEPERATE_FILES:
        sim_datafile_prefix = output_fn
        fn.write(",".join([str(x) for x in sim_data.spins])+"\n")
        write_dlm(f"{sim_datafile_prefix}_x.txt", sim_data.positions, ",")
        write_dlm(f"{sim_datafile_prefix}_xwrap.txt", sim_data.wrapped_positions, ",")
        write_dlm(f"{sim_datafile_prefix}_S.txt", sim_data.spins, ",")
        write_dlm(f"{sim_datafile_prefix}_t.txt", sim_data.times, ",")
    
    fn.close()






# sim_d = loadsim("./basic_N10_t1000.0_align_T1_sim_simdata.txt", DataFileType.SEQUENTIAL_TEXT)

# print()
# print("Positions")
# print(sim_d.positions)
# print()
# print("times")
# print(sim_d.times)
# print()
# print("wrapped")
# print(sim_d.wrapped_positions)
# print()
# print("spins")
# print(sim_d.spins)
# print()
# print("sim params")
# print(sim_d._sim_params)