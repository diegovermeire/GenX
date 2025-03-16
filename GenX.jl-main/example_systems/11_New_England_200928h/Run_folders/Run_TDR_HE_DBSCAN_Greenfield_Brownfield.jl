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
include(raw"C:\Users\Diego\GenX\GenX.jl-main\example_systems\11_New_England_200928h\Compress_TDR_HO_HE.jl")
include(raw"C:\Users\Diego\GenX\GenX.jl-main\example_systems\11_New_England_200928h\Compress_TDR_HO_HE_per_rep_period.jl")

"""
Run file Inputs
"""
# Files
unique_results_folder_GF = "28_02_TDR_DBSCAN_Greenfield_CO2=2_MC=0_UC=0_LDS=1_V4"
results_folder_BF_writing = "28_02_TDR_DBSCAN_Borwnfield_CO2=2_MC=0_UC=0_LDS=1_V4"
base_folder = "resource_enforce_TDR__DBSCAN_V6"
TDR_folder = "TDR_DBSCAN_V5"

# Scenario settings
CO2 = 2
mincap = 0
commit = 0 

# # Data compression settings
# Timesteps_per_period = [8736]
# number_rep_period = Dict(
#     "8736" => [2]
# )
# epsilon_values = [0.001, 0.45, 1.0]

# Data compression settings
Timesteps_per_period = 8736
number_rep_period =[4]
eps_per_number_rep_period = Dict(
    "4" => [0.15]
)
# number_rep_period =[1, 2, 4, 16, 23]
# eps_per_number_rep_period = Dict(
#     "1" => [0.001, 0.15, 0.334000, 0.364273, 0.382247, 0.394545, 0.404636, 0.434909, 0.475273, 0.6, 0.7, 1.0],
#     "2" => [0.001, 0.15, 0.334000, 0.364273, 0.382247, 0.394545, 0.404636, 0.434909, 0.475273, 0.6, 0.7, 1.0],
#     "4" => [0.001, 0.15, 0.334000, 0.364273, 0.382247, 0.394545, 0.404636, 0.434909, 0.475273, 0.6, 0.7, 1.0],
#     "16" => [0.001, 0.15, 0.334000, 0.364273, 0.382247, 0.404636, 0.475273, 0.6, 0.7, 1.0],
#     "23" => [0.001, 0.364273, 0.475273, 1.0]
# )

"""
Model definition 
"""
# Generate all combinations of Timesteps_per_period and number_rep_period
combinations = [(num_rep, eps) for num_rep in number_rep_period for eps in eps_per_number_rep_period[string(num_rep)]]
comb_df = DataFrame(number_rep_period = [x[1] for x in combinations], Epsilon = [x[2] for x in combinations])


# # Generate all combinations of Timesteps_per_period and number_rep_period
# combinations = [(timestep, num_rep) for timestep in Timesteps_per_period for num_rep in number_rep_period[string(timestep)]]
# comb_df = DataFrame(Timestep = [x[1] for x in combinations], NumberRepPeriod = [x[2] for x in combinations])

# # Parse worker ID and number of workers from ARGS
# worker_id = parse(Int, ARGS[1]) 
# num_workers = parse(Int, ARGS[2]) 

# # Distribute the combinations among workers
# assigned_rows = comb_df[worker_id + 1:num_workers:end, :]


for row in eachrow(comb_df)
    println("row", row)
    time_rep_period = Timesteps_per_period
    number_repr_period = row.number_rep_period
    eps = row.Epsilon
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
    mysetup["HeterogenousTimesteps"] = 1
    mysetup["TimeDomainReduction"] = 1
    if eps == 0.001
        mysetup["HeterogenousTimesteps"] = 0
        if number_repr_period == 23
            mysetup["TimeDomainReduction"] = 0
        end
    end
    mysetup["UCommit"] = commit
    mysetup["CO2Cap"] = CO2
    mysetup["MinCapReq"] = mincap
    mysetup["SystemFolder"] = "system"
    ############ Modified ############
    

    ### Cluster time series inputs if necessary and if specified by the user
    if (mysetup["TimeDomainReduction"] == 1) && (mysetup["HeterogenousTimesteps"] == 1)
        mysetup["TimeDomainReductionFolder"] = string(TDR_folder, "_qt_rep_period=", number_repr_period, "_length_rep_period=", time_rep_period)
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
        else
            println("Time Series Data Already Clustered.!")
        end
        path_HO_TDR_HE_DBSCAN = string(TDRpath, "_eps")
        HO_TDR_to_HE_DBSCAN_seperate([eps], TDRpath, path_HO_TDR_HE_DBSCAN)
        mysetup["TimeDomainReductionFolder"] = joinpath(string(mysetup["TimeDomainReductionFolder"],"_eps"), string("system_epsilon_", eps))
    end
    
    ### Cluster time series inputs if necessary and if specified by the user
    if (mysetup["TimeDomainReduction"] == 1) && (mysetup["HeterogenousTimesteps"] == 0)
        mysetup["TimeDomainReductionFolder"] = joinpath(TDR_folder, string("_qt_rep_period=", number_repr_period, "_length_rep_period=", time_rep_period), string("system_epsilon_", eps))
        
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
        else
            println("Time Series Data Already Clustered.!")
        end
    end

    ## Configure solver
    println("Configuring Solver")
    OPTIMIZER = configure_solver(settings_path, optimizer)

    ## Load inputs
    println("Loading Inputs")
    myinputs = GenX.load_inputs(mysetup, case)

    ##Generate the model
    println("Generating the Optimization Model")
    time_elapsed = @elapsed EP = GenX.generate_model(mysetup, myinputs, OPTIMIZER)
    println("Time elapsed for model building is")
    println(time_elapsed)

    println("Solving Model")
    EP, solve_time = GenX.solve_model(EP, mysetup)
    myinputs["solve_time"] = GenX.solve_time # Store the model solve time in myinputs

    outputs_path_GF = ""
    try
        println("Generating outputs for number_repr_period $number_rep_period and number of timeseps per per representative period $time_rep_period")

        ############ Needs to be changed ############
        run_results_folder = string(GenX.get_default_output_folder(case), "_", unique_results_folder_GF)
        try
            if !isdir(run_results_folder)
                mkpath(run_results_folder)
            end
        catch e
            if e isa IOError && e.code == Base.UV_EEXIST
                println("Directory already exists: $run_results_folder")
            else
                println("Directory already exists: $run_results_folder")
            end
        end
        outputs_path = joinpath(run_results_folder, string( string("_number_rep_period=", number_repr_period), string("_time_rep_period=",time_rep_period)), string("eps_", eps))
        println("Outputs path", outputs_path)
        #############################################

        outputs_path = replace(outputs_path, '\\' => '/')

        if mysetup["OverwriteResults"] == 1
            # Overwrite existing results if dir exists
            # This is the default behaviour when there is no flag, to avoid breaking existing code
            if !(isdir(path))
                mkpath(path)
            end
        else
            # Find closest unused ouput directory name and create it
            outputs_path = choose_output_dir(outputs_path)
        end

        elapsed_time = @elapsed outputs_path_GF = GenX.write_outputs(EP,
            outputs_path,
            mysetup,
            myinputs)
        println("Time elapsed for writing is")
        println(elapsed_time)

        # Save results
        accuracy_run_time = NamedTuple[] 
        push!(accuracy_run_time, (CO2, mincap, commit, number_repr_period, time_rep_period, eps, total_cost= value(EP[:eObj]), Run_time=solve_time, Time_Indexes=myinputs["T"], Build_model_time=time_elapsed, write_outputs_time=elapsed_time, termination_Status=termination_status, has_primal=JuMP.has_values(EP)))

        # Convert rows into a DataFrame
        data = DataFrame(accuracy_run_time)
        directory = outputs_path
        base_name = "epsilon_values"
        extension = ".csv"
        version = 0

        # Generate a unique file name
        file_name = ""
        while true
            version_suffix = version > 0 ? "_V$version" : ""
            file_name = joinpath(directory, "$base_name$version_suffix$extension")
            if !isfile(file_name)
                break
            end
            version += 1
        end

        # Save the DataFrame to the unique file
        CSV.write(file_name, data, bufsize=50_000_000)
        CSV.write(joinpath(outputs_path,"mysetup.csv"), mysetup, bufsize=50_000_000)
        CSV.write(joinpath(outputs_path,"myinputs.csv"), myinputs, bufsize=50_000_000)

    catch e 
        println("Writing Outputs")
        termination_status = JuMP.termination_status(EP)
        println("Termination status: ", termination_status)
        println("Has primal values: ", JuMP.has_values(EP))
        ############ Needs to be changed ############
        run_results_folder = string(GenX.get_default_output_folder(case), "_", unique_results_folder_GF)
        if !isdir(run_results_folder)
            mkpath(run_results_folder)
        end
        outputs_path = joinpath(run_results_folder, string( string("_number_rep_period=", number_repr_period), string("_time_rep_period=",time_rep_period)), string("eps_", eps))
        println("Outputs path", outputs_path)
        #############################################
        
        if !(isdir(outputs_path))
            mkpath(outputs_path)
        end

        accuracy_run_time = NamedTuple[]
        push!(accuracy_run_time, (CO2, mincap, commit, number_repr_period, time_rep_period, eps, total_cost= value(EP[:eObj]), Run_time=0, Time_Indexes=myinputs["T"], termination_Status=termination_status, has_primal=JuMP.has_values(EP)))

        # Convert rows into a DataFrame
        data = DataFrame(accuracy_run_time)
        directory = outputs_path
        base_name = "epsilon_values"
        extension = ".csv"
        version = 0

        # Generate a unique file name
        file_name = ""
        while true
            version_suffix = version > 0 ? "_V$version" : ""
            file_name = joinpath(directory, "$base_name$version_suffix$extension")
            if !isfile(file_name)
                break
            end
            version += 1
        end

        # Save the DataFrame to the unique file
        CSV.write(file_name, data, bufsize=50_000_000)
        CSV.write(joinpath(outputs_path,"mysetup.csv"), mysetup, bufsize=50_000_000)
        CSV.write(joinpath(outputs_path,"myinputs.csv"), myinputs, bufsize=50_000_000)

    end

    ##############################################################################
    ################################# BROWNFIELD ################################# 
    ##############################################################################

    outputs_path_read_BF = outputs_path_GF

    # unique_results_folder_reading_cap =  "TDR_28_01_11.03am_MinCAP=1_LDS=1"

    # run_results_folder = string(GenX.get_default_output_folder(case), "_", unique_results_folder_reading_cap)
    # outputs_path_read = joinpath(run_results_folder, string( string("_number_rep_period=", number_repr_period), string("_time_rep_period=",time_rep_period), string("_ExtremePeriod=", extreme_period)))
    # println("outputs_path_read", outputs_path_read)

    res_folder = create_resource_policies_enforce_cap(dirname(dirname(@__FILE__)), outputs_path_read_BF, base_folder)
    dir_res = extract_subpath(res_folder, base_folder)

    mysetup["SystemFolder"] = "system"
    mysetup["TimeDomainReduction"] = 0
    mysetup["HeterogenousTimesteps"] = 0
    mysetup["ResourcesFolder"] = dir_res
    mysetup["PoliciesFolder"] = "policies"

    ## Configure solver
    println("Configuring Solver")
    OPTIMIZER = configure_solver(settings_path, optimizer)

    ## Load inputs
    println("Loading Inputs")
    myinputs = GenX.load_inputs(mysetup, case)

    ##Generate the model
    println("Generating the Optimization Model")
    time_elapsed = @elapsed EP = GenX.generate_model(mysetup, myinputs, OPTIMIZER)
    println("Time elapsed for model building is")
    println(time_elapsed)

    println("Solving Model")
    EP, solve_time = GenX.solve_model(EP, mysetup)
    myinputs["solve_time"] = GenX.solve_time # Store the model solve time in myinputs

    try
        println("Generating outputs for number_repr_period $number_rep_period and number of timeseps per per representative period $time_rep_period")

        ############ Needs to be changed ############
        run_results_folder = string(GenX.get_default_output_folder(case), "_", results_folder_BF_writing)
        try
            if !isdir(run_results_folder)
                mkpath(run_results_folder)
            end
        catch e
            if e isa IOError && e.code == Base.UV_EEXIST
                println("Directory already exists: $run_results_folder")
            else
                println("Directory already exists: $run_results_folder")
            end
        end
        outputs_path = joinpath(run_results_folder, string( string("_number_rep_period=", number_repr_period), string("_time_rep_period=",time_rep_period)), string("eps_", eps))
        println("Outputs path", outputs_path)
        #############################################

        outputs_path = replace(outputs_path, '\\' => '/')

        if mysetup["OverwriteResults"] == 1
            # Overwrite existing results if dir exists
            # This is the default behaviour when there is no flag, to avoid breaking existing code
            if !(isdir(path))
                mkpath(path)
            end
        else
            # Find closest unused ouput directory name and create it
            outputs_path = choose_output_dir(outputs_path)
        end

        elapsed_time = @elapsed outputs_path = GenX.write_outputs(EP,
            outputs_path,
            mysetup,
            myinputs)
        println("Time elapsed for writing is")
        println(elapsed_time)

        directory = outputs_path_GF
        base_name = "epsilon_values.csv"
        filename = joinpath(directory, base_name)

        println("filename", filename)

        accuracy_run_time = CSV.read(filename, DataFrame)
        println("accuracy_run_time", accuracy_run_time)
        new_row = DataFrame(
            :CO2 => [CO2],
            :mincap => [mincap],
            :commit => [commit],
            :number_repr_period => [number_repr_period], 
            :time_rep_period => [time_rep_period],  
            :eps => [eps], 
            :total_cost => [value(EP[:eObj])], 
            :Run_time => [solve_time], 
            :Time_Indexes => [myinputs["T"]], 
            :Build_model_time => [time_elapsed], 
            :write_outputs_time => [elapsed_time], 
            :termination_Status => [string(termination_status)], 
            :has_primal => [JuMP.has_values(EP)]
        )
        append!(accuracy_run_time, new_row)
        # Save the DataFrame to the unique file
        data = DataFrame(accuracy_run_time)
        CSV.write(filename, data, bufsize=50_000_000)
        CSV.write(joinpath(outputs_path,"mysetup.csv"), mysetup, bufsize=50_000_000)
        CSV.write(joinpath(outputs_path,"myinputs.csv"), myinputs, bufsize=50_000_000)
    catch e 
        println("Writing Outputs")
        termination_status = JuMP.termination_status(EP)
        println("Termination status: ", termination_status)
        println("Has primal values: ", JuMP.has_values(EP))
        ############ Needs to be changed ############
        run_results_folder = string(GenX.get_default_output_folder(case), "_", results_folder_BF_writing)
        if !isdir(run_results_folder)
            mkpath(run_results_folder)
        end
        outputs_path = joinpath(run_results_folder, string( string("_number_rep_period=", number_repr_period), string("_time_rep_period=",time_rep_period)), string("eps_", eps))
        println("Outputs path", outputs_path)
        #############################################
        
        if !(isdir(outputs_path))
            mkpath(outputs_path)
        end

        directory = outputs_path_GF
        base_name = "epsilon_values.csv"
        filename = joinpath(directory, base_name)

        println("filename", filename)

        accuracy_run_time = CSV.read(filename, DataFrame)
        println("accuracy_run_time", accuracy_run_time)
        new_row = DataFrame(
            :CO2 => [CO2],
            :mincap => [mincap],
            :commit => [commit],
            :number_repr_period => [number_repr_period], 
            :time_rep_period => [time_rep_period],  
            :eps => [eps], 
            :total_cost => [0], 
            :Run_time => [0], 
            :Time_Indexes => [myinputs["T"]], 
            :Build_model_time => [time_elapsed], 
            :write_outputs_time => [0], 
            :termination_Status => [string(termination_status)], 
            :has_primal => [JuMP.has_values(EP)]
        )
        append!(accuracy_run_time, new_row)
        # Save the DataFrame to the unique file
        data = DataFrame(accuracy_run_time)
        CSV.write(filename, data, bufsize=50_000_000)
        CSV.write(joinpath(outputs_path,"mysetup.csv"), mysetup, bufsize=50_000_000)
        CSV.write(joinpath(outputs_path,"myinputs.csv"), myinputs, bufsize=50_000_000)
    end
end
