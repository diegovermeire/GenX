##
ENV["GENX_PRECOMPILE"] = "false"
using GenX
using Gurobi
using JuMP
using DataFrames       
using CSV              
using FilePathsBase
using YAML
include(raw"C:\Users\Diego\GenX\GenX.jl-main\Work directrory\TDR_HO_HE.jl")
include(raw"C:\Users\Diego\GenX\GenX.jl-main\Work directrory\TDR_HO_HO.jl")
include(raw"C:\Users\Diego\GenX\GenX.jl-main\Work directrory\Analyse_HO_HE.jl")

# Timesteps_per_period = collect(range(24, stop=1008, length=25))
# Timesteps_per_period = map(x -> round(Int, x), Timesteps_per_period)
unique_results_folder = "TDR_no_sequential"

Timesteps_per_period = [48]
number_rep_period = Dict(
    "48" => [8, 40, 80, 140, 160, 182],
    # "168" => [4, 8],
    # "672" => [2]
    )

# number_rep_period = Dict(
#     "48" => [4, 8, 40, 80, 182],
#     "168" => [4, 8, 12, 24, 36, 48],
#     "672" => [2, 4, 8, 12]
#     )


# # Collect accuracy vs run_time results
# accuracy_run_time = Dict(
#     time_rep_period => Dict(
#         number_repr_period => Dict("total_cost" => 0.0, "Run_time" => 0.0)
#         for number_repr_period in number_rep_period[string(time_rep_period)]
#     )
#     for time_rep_period in Timesteps_per_period
# )


for time_rep_period in Timesteps_per_period
    time_rep_period
    ## Define variables
    path = dirname(@__FILE__)
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

    for number_repr_period in number_rep_period[string(time_rep_period)]
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
                println("Modified the MaxPeriods and MinPeriods variables of myTDRsetup")
                random = true
                GenX.cluster_inputs(case, settings_path, mysetup,  -99,false; random = true, myTDRsetup_ = myTDRsetup)
                # TDR_HO_to_HO(TDRpath)
            else
                println("Time Series Data Already Clustered.!")
            end
        end

        ### Cluster time series inputs if necessary and if specified by the user
        if (mysetup["TimeDomainReduction"] == 1) && (mysetup["HeterogenousTimesteps"] == 1)        
            mysetup["TimeDomainReductionFolder"] = "TDR_results_HE"
            TDRpath = joinpath(case, mysetup["TimeDomainReductionFolder"])
            system_path = joinpath(case, mysetup["SystemFolder"])
            GenX.prevent_doubled_timedomainreduction(system_path)

            #Verify whether TRD files exists at TRDpath
            if !GenX.time_domain_reduced_files_exist(TDRpath)
                println("Clustering Time Series Data (Grouped)...")
                destination_path = replace(TDRpath_homogenous, "HO" => "HE")
                # Check if the source folder exists
                if isdir(TDRpath_homogenous)
                    try
                        # Copy the folder
                        cp(TDRpath_homogenous, destination_path; force=true)
                        println("Folder successfully copied and renamed from '$TDRpath' to '$destination_path'.")
                    catch e
                        println("An error occurred: ", e)
                    end
                else
                    println("The folder at '$TDRpath' does not exist.")
                end
                # TDR_HO_to_HE(TDRpath)
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

        try
            println("Generating outputs for number_repr_period $number_rep_period and number of timeseps per per representative period $time_rep_period")

            ############ Needs to be changed ############
            run_results_folder = string(GenX.get_default_output_folder(case), "_", unique_results_folder)
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
            outputs_path = joinpath(run_results_folder, string( string("_number_rep_period=", number_repr_period), string("_time_rep_period=",time_rep_period)))
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

            # Save results
            accuracy_run_time = NamedTuple[] 
            push!(accuracy_run_time, (number_repr_period, time_rep_period, total_cost= value(EP[:eObj]), Run_time=solve_time, Time_Indexes=myinputs["T"], Build_model_time=time_elapsed, write_outputs_time=elapsed_time, termination_Status=termination_status, has_primal=JuMP.has_values(EP)))

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
            CSV.write(file_name, data)

        catch e 
            println("Writing Outputs")
            termination_status = JuMP.termination_status(EP)
            println("Termination status: ", termination_status)
            println("Has primal values: ", JuMP.has_values(EP))
            ############ Needs to be changed ############
            run_results_folder = string(GenX.get_default_output_folder(case), "_", unique_results_folder)
            if !isdir(run_results_folder)
                mkpath(run_results_folder)
            end
            outputs_path = joinpath(run_results_folder, string( string("_number_rep_period=", number_repr_period), string("_time_rep_period=",time_rep_period)))
            println("Outputs path", outputs_path)
            #############################################
            
            if !(isdir(outputs_path))
                mkpath(outputs_path)
            end

            accuracy_run_time = NamedTuple[]
            push!(accuracy_run_time, (number_repr_period, time_rep_period, total_cost= value(EP[:eObj]), Run_time=0, Time_Indexes=myinputs["T"], termination_Status=termination_status, has_primal=JuMP.has_values(EP)))

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
            CSV.write(file_name, data)

        end
    end
end

# # Write the saved accuracy vs. run_time results in a .csv file
# data = [
#     (time_rep_period, number_repr_period, dict["total_cost"], dict["Run_time"])
#     for (time_rep_period, inner_dict) in accuracy_run_time
#     for (number_repr_period, dict) in inner_dict
# ]

# # Create a DataFrame from the collected data
# df = DataFrame(
#     time_rep_period = [row[1] for row in data],
#     number_repr_period = [row[2] for row in data],
#     total_cost = [row[3] for row in data],
#     Run_time = [row[4] for row in data]
# )
# # Get the directory of the current script
# directory = @__DIR__

# # Base file name and extension
# base_name = "accuracy_run_time_TDR"
# extension = ".csv"
# version = 0

# # Generate a unique file name
# file_name = ""
# while true
#     version_suffix = version > 0 ? "_V$version" : ""
#     file_name = joinpath(directory, "$base_name$version_suffix$extension")
#     if !isfile(file_name)
#         break
#     end
#     version += 1
# end

# # Save the DataFrame to the unique file
# CSV.write(file_name, df)