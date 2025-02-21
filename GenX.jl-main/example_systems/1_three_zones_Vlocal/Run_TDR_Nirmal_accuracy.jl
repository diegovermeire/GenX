##
"""
Once the different scenarios using TDR (number_rep_period and UseExtremePeriods) have been computed, this file:
- rewrites one policies and resources folder for each TDR run by enforcing the installed capacities
of the TDR run into the rewritten policies and resources folder_min_cap
- runs GenX model for the rewritten policies and resources folder enforcing the installed capacities
  using MinCapReq, over the entire raw data (not using compressed raw data)
- saves the objective function of the above GenX run within the epsilon_values.csv folder of the
  corresponding TDR run
"""
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
include(raw"C:\Users\Diego\GenX\GenX.jl-main\Work directrory\create_resources_policies_enforce_cap.jl")

UseExtremePeriods =[1]
unique_results_folder_writing_enforce = "TDR_04_02_8.03pm_enforce_TDR"
unique_results_folder_reading_cap =  "TDR_28_01_11.03am_MinCAP=1_LDS=1"
base_folder = "resource_enforce_TDRV3"
Timesteps_per_period = [48, 168, 672]
number_rep_period = Dict(
    "48" => [8, 40, 80, 182],
    "168" => [8, 12, 24, 36, 48],
    "672" => [4, 8, 12]
    )


for time_rep_period in Timesteps_per_period
    # time_rep_period
    # ## Define variables
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
    mysetup["TimeDomainReduction"] = 0
    mysetup["MinCapReq"] = 1

    for extreme_period in UseExtremePeriods
        for number_repr_period in number_rep_period[string(time_rep_period)]
            run_results_folder = string(GenX.get_default_output_folder(case), "_", unique_results_folder_reading_cap)
            outputs_path_read = joinpath(run_results_folder, string( string("_number_rep_period=", number_repr_period), string("_time_rep_period=",time_rep_period), string("_ExtremePeriod=", extreme_period)))
            println("outputs_path_read", outputs_path_read)

            res_folder = create_resource_policies_enforce_cap(dirname(@__FILE__), outputs_path_read, base_folder)
            dir_res = extract_subpath(res_folder, base_folder)
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
                run_results_folder = string(GenX.get_default_output_folder(case), "_", unique_results_folder_writing_enforce)
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
                outputs_path = joinpath(run_results_folder, string( string("_number_rep_period=", number_repr_period), string("_time_rep_period=",time_rep_period), string("_ExtremePeriod=", extreme_period)))
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

                directory = outputs_path_read
                base_name = "epsilon_values.csv"
                filename = joinpath(directory, base_name)

                println("filename", filename)

                accuracy_run_time = CSV.read(filename, DataFrame)
                println("accuracy_run_time", accuracy_run_time)
                new_row = DataFrame(number_repr_period=[number_repr_period], time_rep_period= [time_rep_period], extreme_period=[extreme_period], total_cost=[value(EP[:eObj])], Run_time=[solve_time], 
                                    Time_Indexes=[myinputs["T"]], Build_model_time=[time_elapsed], write_outputs_time=[elapsed_time], 
                                    termination_Status=[string(termination_status)], has_primal=[JuMP.has_values(EP)])
                append!(accuracy_run_time, new_row)
                # Save the DataFrame to the unique file
                data = DataFrame(accuracy_run_time)
                CSV.write(filename, data)

            catch e 
                println("Writing Outputs")
                termination_status = JuMP.termination_status(EP)
                println("Termination status: ", termination_status)
                println("Has primal values: ", JuMP.has_values(EP))
                ############ Needs to be changed ############
                run_results_folder = string(GenX.get_default_output_folder(case), "_", unique_results_folder_writing_enforce)
                if !isdir(run_results_folder)
                    mkpath(run_results_folder)
                end
                outputs_path = joinpath(run_results_folder, string( string("_number_rep_period=", number_repr_period), string("_time_rep_period=",time_rep_period), string("_ExtremePeriod=", extreme_period)))
                println("Outputs path", outputs_path)
                #############################################
                
                if !(isdir(outputs_path))
                    mkpath(outputs_path)
                end

                directory = outputs_path_read
                base_name = "epsilon_values.csv"
                filename = joinpath(directory, base_name)

                println("filename", filename)

                accuracy_run_time = CSV.read(filename, DataFrame)
                println("accuracy_run_time", accuracy_run_time)
                new_row = DataFrame(epsilon=[eps], total_cost=[0], Run_time=[solve_time], 
                                    Time_Indexes=[myinputs["T"]], Build_model_time=[time_elapsed], write_outputs_time=[elapsed_time], 
                                    termination_Status=[string(termination_status)], has_primal=[JuMP.has_values(EP)])
                append!(accuracy_run_time, new_row)
                # Save the DataFrame to the unique file
                data = DataFrame(accuracy_run_time)
                CSV.write(filename, data)
            end
        end
    end
end