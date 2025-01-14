##
ENV["GENX_PRECOMPILE"] = "false"
using GenX
using Gurobi
using JuMP
include(raw"C:\Users\Diego\GenX\GenX.jl-main\Work directrory\TDR_HO_HE.jl")
include(raw"C:\Users\Diego\GenX\GenX.jl-main\Work directrory\TDR_HO_HO.jl")
include(raw"C:\Users\Diego\GenX\GenX.jl-main\Work directrory\Analyse_HO_HE.jl")

############ Needs to be changed ############
epsilon_values =[1.        , 0.98990909, 0.97981818, 0.96972727, 0.95963636, 
0.94954545, 0.93945455, 0.92936364, 0.91927273, 0.90918182, 
0.89909091, 0.889     , 0.87890909, 0.86881818, 0.85872727, 
0.84863636, 0.83854545, 0.82845455, 0.81836364, 0.80827273]
# epsilon_values = [0.667     , 0.65690909, 
# 0.64681818, 0.63672727, 0.62663636, 0.61654545, 0.60645455, 
# 0.59636364, 0.58627273, 0.57618182, 0.56609091, 0.556     , 
# 0.54590909, 0.53581818, 0.52572727, 0.51563636, 0.50554545, 
# 0.49545455, 0.48536364, 0.47527273, 0.46518182, 0.45509091, 
# 0.445     , 0.43490909, 0.42481818, 0.41472727, 0.40463636, 
# 0.39454545, 0.38445455, 0.37436364, 0.36427273, 0.35418182, 
# 0.34409091, 0.334     , 0.32390909, 0.31381818, 0.30372727, 
# 0.29363636, 0.28354545, 0.27345455, 0.26336364, 0.25327273, 
# 0.24318182, 0.23309091, 0.223     , 0.21290909, 0.20281818, 
# 0.19272727, 0.18263636, 0.17254545, 0.16245455, 0.15236364, 
# 0.14227273, 0.13218182, 0.12209091, 0.112     , 0.10190909, 
# 0.09181818, 0.08172727, 0.07163636, 0.06154545, 0.05145455, 
# 0.04136364, 0.03127273, 0.02118182, 0.01109091, 0.001     ]

#############################################

accuracy_run_time = Dict(epsilon => Dict("total_cost" => 0.0, "Run_time" => 0.0) for epsilon in epsilon_values)

for eps in epsilon_values
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


    ############ Changed ############
    mysetup["SystemFolder"] = string(mysetup["SystemFolder"], "_epsilon_", eps)
    #############################################

    case= path
    optimizer = Gurobi.Optimizer
    settings_path = GenX.get_settings_path(case)
    outputs_folder = ["", ""]
    TDRpath_homogenous = ""

    i=1

    if (mysetup["TimeDomainReduction"] == 1) && (mysetup["HeterogenousTimesteps"] == 0)
        mysetup["TimeDomainReductionFolder"] = "TDR_results_HO"
        TDRpath = joinpath(case, mysetup["TimeDomainReductionFolder"])
        TDRpath_homogenous = TDRpath
        system_path = joinpath(case, mysetup["SystemFolder"])
        GenX.prevent_doubled_timedomainreduction(system_path)

        #Verify whether TRD files exists at TRDpath
        if !GenX.time_domain_reduced_files_exist(TDRpath)
            println("Clustering Time Series Data (Grouped)...")
            GenX.cluster_inputs(case, settings_path, mysetup)
            TDR_HO_to_HO(TDRpath)
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
            TDR_HO_to_HE(TDRpath)
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

    println("Writing Outputs")
    ############ Needs to be changed ############
    outputs_path = string(GenX.get_default_output_folder(case), "_epsilon_", eps)
    if !isdir(outputs_path)
        mkdir(outputs_path)
    end
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

    outputs_folder[i] = outputs_path
    elapsed_time = @elapsed outputs_path = GenX.write_outputs(EP,
        outputs_path,
        mysetup,
        myinputs)
    println("Time elapsed for writing is")
    println(elapsed_time)

    # Save results 
    accuracy_run_time[eps]["total_cost"] = value(EP[:eObj])
    accuracy_run_time[eps]["Run_time"] = solve_time
end

using DataFrames       # For creating and working with DataFrame
using CSV              # For writing the DataFrame to a CSV file
using FilePathsBase    # Provides `joinpath` if you use `os.joinpath`
data = DataFrame(
    epsilon = [ key for key in  keys(accuracy_run_time)],
    total_cost = [get(dict, "total_cost", nothing) for dict in values(accuracy_run_time)],
    Run_time = [get(dict, "Run_time", nothing) for dict in values(accuracy_run_time)]
)

# Get the directory of the current script
directory = @__DIR__

# Base file name and extension
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

