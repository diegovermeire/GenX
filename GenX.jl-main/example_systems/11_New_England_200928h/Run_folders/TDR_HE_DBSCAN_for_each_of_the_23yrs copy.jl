using Base.Filesystem
include(raw"C:\Users\Diego\GenX\GenX.jl-main\example_systems\11_New_England_200928h\Compress_TDR_HO_HE.jl")

function list_subfolders(parent_dir::String)
    subfolders = []
    
    for (root, dirs, _) in walkdir(parent_dir)
        append!(subfolders, joinpath.(root, dirs))  # Get full path of each subfolder
    end
    
    return subfolders
end

# Example usage
parent_directory = raw"C:\Users\Diego\GenX\GenX.jl-main\example_systems\11_New_England_200928h_system_yr"
subfolders = list_subfolders(parent_directory)
epsilon_values_36 = [0.001 , 0.02559659, 0.04404403, 0.08093892,  0.10553551, 0.12398295, 0.1424304, 0.16087784, 0.19162358, 0.23466761, 0.27771164, 0.31460653, 
0.38224715, 0.57974331, 1.0, 0.20281818, 0.23309091, 0.27345455,  0.29363636, 0.32390909, 0.334, 0.36427273,  0.39454545, 0.40463636, 0.43490909, 0.47527273, 0.50554545,  0.53581818, 0.56609091, 0.59636364,
0.62663636, 0.65690909, 0.68718182, 0.71745455,  0.74772727, 0.76790909]

for folder in subfolders
    rep_period = basename(folder)
    save_folder = raw"C:\Users\Diego\GenX\GenX.jl-main\example_systems\11_New_England_200928h\HE_TDR_DBSCAN_InputData"
    save_path = joinpath(save_folder, "11_New_England_200928h_system_yr", string(rep_period, "_36eps"))

    HO_TDR_to_HE_DBSCAN(epsilon_values_36, folder, save_path)
end