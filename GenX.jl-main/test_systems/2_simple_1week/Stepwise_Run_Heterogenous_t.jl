##
ENV["GENX_PRECOMPILE"] = "false"
using GenX

## Define variables
path = dirname(@__FILE__)
genx_settings = GenX.get_settings_path(path, "genx_settings.yml") # Settings YAML file path
writeoutput_settings = GenX.get_settings_path(path, "output_settings.yml") # Write-output settings YAML file path
mysetup = GenX.configure_settings(genx_settings, writeoutput_settings) # mysetup dictionary stores settings and GenX-specific parameters

########### Adapt Manually when playing with HE/HO timsteps ###########
mysetup["HeterogenousTimesteps"] = 1
#######################################################################

# if mysetup["MultiStage"] == 0
#     run_genx_case_simple!(case, mysetup, optimizer)
# else
#     run_genx_case_multistage!(case, mysetup, optimizer)
# end

#Run: run_genx_case_simple

case= path
optimizer = GenX.HiGHS.Optimizer

settings_path = GenX.get_settings_path(case)

### Cluster time series inputs if necessary and if specified by the user
if mysetup["TimeDomainReduction"] == 1
    TDRpath = joinpath(case, mysetup["TimeDomainReductionFolder"])
    system_path = joinpath(case, mysetup["SystemFolder"])
    GenX.prevent_doubled_timedomainreduction(system_path)
    if !GenX.time_domain_reduced_files_exist(TDRpath)
        println("Clustering Time Series Data (Grouped)...")
        GenX.cluster_inputs(case, settings_path, mysetup)
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
myinputs

##Generate the model
println("Generating the Optimization Model")
time_elapsed = @elapsed EP = GenX.generate_model(mysetup, myinputs, OPTIMIZER)
println("Time elapsed for model building is")
println(time_elapsed)

println("Solving Model")
EP, solve_time = GenX.solve_model(EP, mysetup)
myinputs["solve_time"] = GenX.solve_time # Store the model solve time in myinputs

println("Writing Outputss")
outputs_path = GenX.get_default_output_folder(case)
elapsed_time = @elapsed outputs_path = GenX.write_outputs(EP,
    outputs_path,
    mysetup,
    myinputs)
println("Time elapsed for writing is")
println(elapsed_time)