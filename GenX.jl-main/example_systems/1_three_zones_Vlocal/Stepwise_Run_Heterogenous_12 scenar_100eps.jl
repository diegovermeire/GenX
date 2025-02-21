##
ENV["GENX_PRECOMPILE"] = "false"
using GenX
using Gurobi
using JuMP
using MathOptInterface
include(raw"C:\Users\Diego\GenX\GenX.jl-main\Work directrory\TDR_HO_HE.jl")
include(raw"C:\Users\Diego\GenX\GenX.jl-main\Work directrory\TDR_HO_HO.jl")
include(raw"C:\Users\Diego\GenX\GenX.jl-main\Work directrory\Analyse_HO_HE.jl")

############ Needs to be changed ############
max_repetition = 1
unique_results_folder = "HE_CO2=0_UC=0_MC=1_eps=0.5"
# resource_folder= "resources_NGCC"
# policies_folder = "policies_NGCC"
resource_folder= "resources"
policies_folder = "policies"
# DEBUGINGGGGGGGGGGGGGGGGGGGGGGGGGG on line 66 ############

CO2Cap = [0]
UCommit = [0]
MinCapReq = [1]

epsilon_values = [0.50970053]

# epsilon_values = [0.001     , 0.00714915, 0.0132983 , 0.01944744, 0.02559659,
#        0.03174574, 0.03789489, 0.04404403, 0.05019318, 0.05634233,
#        0.06249148, 0.06864062, 0.07478977, 0.08093892, 0.08708807,
#        0.09323721, 0.09938636, 0.10553551, 0.11168466, 0.11783381,
#        0.12398295, 0.1301321 , 0.13628125, 0.1424304 , 0.14857954,
#        0.15472869, 0.16087784, 0.16702699, 0.17317613, 0.17932528,
#        0.18547443, 0.19162358, 0.19777273, 0.20392187, 0.21007102,
#        0.21622017, 0.22236932, 0.22851846, 0.23466761, 0.24081676,
#        0.24696591, 0.25311505, 0.2592642 , 0.26541335, 0.2715625 ,
#        0.27771164, 0.28386079, 0.29000994, 0.29615909, 0.30230824,
#        0.30845738, 0.31460653, 0.32075568, 0.32690483, 0.33305397,
#        0.33920312, 0.34535227, 0.35150142, 0.35765056, 0.36379971,
#        0.36994886, 0.37609801, 0.38224715, 0.3883963 , 0.39454545,
#        0.40463636, 0.42214706, 0.43965775, 0.45716845, 0.47467914,
#        0.49218984, 0.50970053, 0.52721123, 0.54472192, 0.56223262,
#        0.57974331, 0.59725401, 0.6147647 , 0.6322754 , 0.64978609,
#        0.66729679, 0.68480748, 0.70231818, 0.71982888, 0.73733957,
#        0.75485027, 0.77236096, 0.78987166, 0.80738235, 0.82489305,
#        0.84240374, 0.85991444, 0.87742513, 0.89493583, 0.91244652,
#        0.92995722, 0.94746791, 0.96497861, 0.9824893 , 1.        ]

#############################################

# accuracy_run_time = Dict(epsilon => Dict("total_cost" => 0.0, "Run_time" => 0.0) for epsilon in epsilon_values)


# for i in (epsilon_values)
#     eps = i
#     for mincap in MinCapReq
#         for commit in UCommit 
#             for CO2 in CO2Cap
#                 for repetition in 1:max_repetition


eps = epsilon_values[1]
mincap = MinCapReq[1]
commit = UCommit[1] 
CO2 = CO2Cap[1]
repetition = max_repetition[1]
accuracy_run_time = NamedTuple[]
## Define variables
path = dirname(@__FILE__)
genx_settings = GenX.get_settings_path(path, "genx_settings.yml") # Settings YAML file path
writeoutput_settings = GenX.get_settings_path(path, "output_settings.yml") # Write-output settings YAML file path
mysetup = GenX.configure_settings(genx_settings, writeoutput_settings) # mysetup dictionary stores settings and GenX-specific parameters

if eps <= 0.002
    mysetup["HeterogenousTimesteps"] = 0
else
    mysetup["HeterogenousTimesteps"] = 1
end


############ DEBUGINGGGGGGGGGGGGGGGGGGGGGGGGGG ############
mysetup["ParameterScale"] = 0
mysetup["SystemFolder"] = string("systems_epsilon/" , mysetup["SystemFolder"], "_epsilon_", eps)
mysetup["ResourcesFolder"] = resource_folder
mysetup["PoliciesFolder"] = policies_folder
mysetup["UCommit"] = commit
mysetup["CO2Cap"] = CO2
mysetup["MinCapReq"] = mincap
#############################################

case= path
optimizer = Gurobi.Optimizer
settings_path = GenX.get_settings_path(case)

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

termination_status = JuMP.termination_status(EP)
println("Termination status: ", termination_status)
println("Has primal values: ", JuMP.has_values(EP))


try
    println("Generating outputs for epsilon: $eps, CO2 cap: $CO2, and UCommit: $commit")

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
    outputs_path = joinpath(run_results_folder, string( string("_CO2=", CO2), string("_MinCapReq=",mincap), string("_COMMIT=", commit), string("_epsilon=", eps), string("_V",repetition)))
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
    push!(accuracy_run_time, (epsilon=eps, total_cost= value(EP[:eObj]), Run_time=solve_time, Time_Indexes=myinputs["T"], Build_model_time=time_elapsed, write_outputs_time=elapsed_time, termination_Status=termination_status, has_primal=JuMP.has_values(EP)))

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
    CSV.write(joinpath(directory,"mysetup.csv"), mysetup)
    CSV.write(joinpath(directory,"myinputs.csv"), myinputs)

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
    outputs_path = joinpath(run_results_folder, string( string("_CO2=", CO2), string("_MinCapReq=",mincap), string("_COMMIT=", commit), string("_epsilon=", eps), string("_V",repetition, "INFEASIBLE")))
    #############################################
    
    if !(isdir(outputs_path))
        mkpath(outputs_path)
    end

    accuracy_run_time = NamedTuple[]
    push!(accuracy_run_time, (epsilon=eps, total_cost= 0, Run_time=solve_time, Time_Indexes=myinputs["T"], Build_model_time=time_elapsed, termination_Status=termination_status, has_primal=JuMP.has_values(EP)))

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
    CSV.write(joinpath(directory,"mysetup.csv"), mysetup)
    CSV.write(joinpath(directory,"myinputs.csv"), myinputs)
end

#                 end
#             end
#         end
#     end
# end





