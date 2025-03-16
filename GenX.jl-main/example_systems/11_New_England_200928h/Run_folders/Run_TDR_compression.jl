# Set environment variable to disable precompilation
ENV["GENX_PRECOMPILE"] = "false"

# Import necessary packages
using GenX
using Gurobi
using JuMP
using DataFrames
using CSV
using FilePathsBase
using YAML

# Include custom scripts
# include(raw"/home/gridsan/dvermeire/GenX.jl-main/Work directrory/TDR_HO_HE.jl")
# include(raw"/home/gridsan/dvermeire/GenX.jl-main/Work directrory/TDR_HO_HO.jl")
# include(raw"/home/gridsan/dvermeire/GenX.jl-main/Work directrory/Analyse_HO_HE.jl")
include(raw"C:\Users\Diego\GenX\GenX.jl-main\Work directrory\create_resources_policies_enforce_cap.jl")

# Define unique result folders and base folder
unique_results_folder_GF = "22_02_Greenfield_yr_TDR_CO2=2_MC=0_UC=0"
results_folder_BF_writing = "22_02_Borwnfield_yr_TDR_CO2=2_MC=0_UC=0"
base_folder = "resource_enforce_TDR_V4"
CO2 = 2
mincap = 0
commit = 0 

# Define your Timesteps and number of representative periods
# Timesteps_per_period = [96, 168, 672, 4368, 8736]
# number_rep_period = Dict(
#     "96" => [6, 8, 40, 80, 182],
#     "168" => [8, 16, 32, 48, 64, 96, 112],
#     "672" => [2, 4, 8, 12, 16, 24, 28],
#     "4368" => [1, 2, 3, 4],
#     "8736" => [1, 2],
# )
Timesteps_per_period = [8736]

number_rep_period = Dict(
    "8736" => [1, 2, 4, 8, 16, 23]
)
# Timesteps_per_period = [48]
# number_rep_period = Dict(
#     "48" => [12]
# )

# Generate all combinations of Timesteps_per_period and number_rep_period
combinations = [(timestep, num_rep) for timestep in Timesteps_per_period for num_rep in number_rep_period[string(timestep)]]
comb_df = DataFrame(Timestep = [x[1] for x in combinations], NumberRepPeriod = [x[2] for x in combinations])

# # Parse worker ID and number of workers from ARGS
# worker_id = parse(Int, ARGS[1]) 
# num_workers = parse(Int, ARGS[2]) 

# # Distribute the combinations among workers
# assigned_rows = comb_df[worker_id + 1:num_workers:end, :]

# Iterate over the assigned rows for each worker
for row in eachrow(comb_df)
    println("row", row)
    time_rep_period = row.Timestep
    number_repr_period = row.NumberRepPeriod
    # println("row", comb_df)
    # time_rep_period = comb_df.Timestep
    # number_repr_period = comb_df.NumberRepPeriod
    
    ## Define variables
    path = dirname(dirname(@__FILE__))
    genx_settings = GenX.get_settings_path(path, "genx_settings.yml") # Settings YAML file path
    writeoutput_settings = GenX.get_settings_path(path, "output_settings.yml") # Write-output settings YAML file path
    mysetup = GenX.configure_settings(genx_settings, writeoutput_settings) # mysetup dictionary stores settings and GenX-specific parameters
    case= path
    optimizer = Gurobi.Optimizer
    settings_path = GenX.get_settings_path(case)
    outputs_folder = ["", ""]
    TDRpath_homogenous = ""

    # Change the TDR_settings
    myTDRsetup = YAML.load(open(joinpath(settings_path, "time_domain_reduction_settings.yml")))
    GenX.update_deprecated_tdr_inputs!(myTDRsetup)
    mysetup["TimestepsPerRepPeriod"] = time_rep_period
    mysetup["HeterogenousTimesteps"] = 0
    mysetup["TimeDomainReduction"] = 1
    mysetup["UCommit"] = commit
    mysetup["CO2Cap"] = CO2
    mysetup["MinCapReq"] = mincap
    mysetup["SystemFolder"] = "system"

    ### Cluster time series inputs if necessary and if specified by the user
    if (mysetup["TimeDomainReduction"] == 1) && (mysetup["HeterogenousTimesteps"] == 0)
        mysetup["TimeDomainReductionFolder"] = string("TDR_results_HO", "_number_repr_period=", number_repr_period, "time_rep_period=", time_rep_period)
        println("#########################################################")
        println("TimeDomainReductionFolder", mysetup["TimeDomainReductionFolder"])
        println("#########################################################")
        TDRpath = joinpath(case, mysetup["TimeDomainReductionFolder"])
        TDRpath_homogenous = TDRpath
        system_path = joinpath(case, mysetup["SystemFolder"])
        GenX.prevent_doubled_timedomainreduction(system_path)

        #Verify whether TRD files exists at TRDpath
        if !GenX.time_domain_reduced_files_exist(TDRpath)
            println("Clustering Time Series Data (Grouped)...")
            myTDRsetup = YAML.load(open(joinpath(settings_path,
            "time_domain_reduction_settings.yml")))
            myTDRsetup["MinPeriods"] = number_repr_period
            myTDRsetup["MaxPeriods"] = number_repr_period
            myTDRsetup["UseExtremePeriods"] = 0
            println("Modified the MaxPeriods and MinPeriods variables of myTDRsetup")
            random = true
            GenX.cluster_inputs(case, settings_path, mysetup,  -99,false; random = true, myTDRsetup_ = myTDRsetup)
            # TDR_HO_to_HO(TDRpath)
        else
            println("Time Series Data Already Clustered.!")
        end
    end
end