##
ENV["GENX_PRECOMPILE"] = "false"
using GenX
using Gurobi
using JuMP
include(raw"C:\Users\Diego\GenX\GenX.jl-main\Work directrory\TDR_HO_HE.jl")
include(raw"C:\Users\Diego\GenX\GenX.jl-main\Work directrory\TDR_HO_HO.jl")

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

i=2
mysetup["HeterogenousTimesteps"] = i-1

### Cluster time series inputs if necessary and if specified by the user
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
    else
        println("Time Series Data Already Clustered.!")
    end
    TDR_HO_to_HO(TDRpath)
end

### Cluster time series inputs if necessary and if specified by the user
if (mysetup["TimeDomainReduction"] == 1) && (mysetup["HeterogenousTimesteps"] == 1)
    TDRpath = "C:\\Users\\Diego\\GenX\\GenX.jl-main\\test_systems\\example_systems_test\\1_three_zones\\TDR_results_HO"
    destination_path = replace(TDRpath_homogenous, "HO" => "HE")
    # Check if the source folder exists
    if isdir(TDRpath)
        try
            # Copy the folder
            cp(TDRpath, destination_path; force=true)
            println("Folder successfully copied and renamed from '$TDRpath' to '$destination_path'.")
        catch e
            println("An error occurred: ", e)
        end
    else
        println("The folder at '$TDRpath' does not exist.")
    end
    
    mysetup["TimeDomainReductionFolder"] = "TDR_results_HE"
    TDRpath = joinpath(case, mysetup["TimeDomainReductionFolder"])
    system_path = joinpath(case, mysetup["SystemFolder"])
    GenX.prevent_doubled_timedomainreduction(system_path)

    mysetup["HeterogenousTimesteps"] = 0
    #Verify whether TRD files exists at TRDpath
    if !GenX.time_domain_reduced_files_exist(TDRpath)
        println("Clustering Time Series Data (Grouped)...")
        GenX.cluster_inputs(case, settings_path, mysetup)
        TDR_HO_to_HE(TDRpath)
    else
        println("Time Series Data Already Clustered.!")
    end
    mysetup["HeterogenousTimesteps"] = 1
end

## Configure solver
println("Configuring Solver")
OPTIMIZER = configure_solver(settings_path, optimizer)

## Load inputs
println("Loading Inputs")
myinputs = GenX.load_inputs(mysetup, case)
myinputs

##Generate the model
println("Generating the Optimization Model")
time_elapsed = @elapsed EP = GenX.generate_model(mysetup, myinputs, OPTIMIZER)
println("Time elapsed for model building is")
println(time_elapsed)

println("Solving Model")
EP, solve_time = GenX.solve_model(EP, mysetup)
myinputs["solve_time"] = GenX.solve_time # Store the model solve time in myinputs

# GenX.write_energy_revenue(raw"C:\Users\Diego\GenX\GenX.jl-main\test_systems\3_simple_1week_easiest\results_1",myinputs, mysetup, EP)

println("Writing Outputss")
outputs_path = GenX.get_default_output_folder(case)
outputs_path = replace(outputs_path, '\\' => '/')
outputs_path

outputs_folder[i] = outputs_path
elapsed_time = @elapsed outputs_path = GenX.write_outputs(EP,
    outputs_path,
    mysetup,
    myinputs)
println("Time elapsed for writing is")
println(elapsed_time)

haskey(myinputs, "dfCO2Cap_slack")