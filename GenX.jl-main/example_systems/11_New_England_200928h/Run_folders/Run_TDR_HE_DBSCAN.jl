using Base.Filesystem
include(raw"C:\Users\Diego\GenX\GenX.jl-main\example_systems\11_New_England_200928h\Compress_TDR_HO_HE.jl")

epsilon_values_2 = [ 1.0, 0.334]
path_TDR = raw"C:\Users\Diego\GenX\GenX.jl-main\example_systems\11_New_England_200928h\HO_TDR_InputData\TDR_results_HO_number_repr_period=2time_rep_period=8736"
path_TDR_HE_DBSCAN = raw"C:\Users\Diego\GenX\GenX.jl-main\example_systems\11_New_England_200928h\HE_TDR_DBSCAN_InputData"
save_path = joinpath(path_TDR_HE_DBSCAN, "TDR_results_HO_number_repr_period=2time_rep_period=8736")

HO_TDR_to_HE_DBSCAN(epsilon_values_2, path_TDR, save_path)
