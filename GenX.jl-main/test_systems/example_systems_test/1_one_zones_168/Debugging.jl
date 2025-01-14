##
using GenX

## Cell 0
case = dirname(@__FILE__)

##Define case and mysetup
genx_settings = GenX.get_settings_path(case, "genx_settings.yml") # Settings YAML file path
writeoutput_settings = GenX.get_settings_path(case, "output_settings.yml") # Write-output settings YAML file path
mysetup = GenX.configure_settings(genx_settings, writeoutput_settings) # mysetup dictionary stores settings and GenX-specific parameters
optimizer = GenX.HiGHS.Optimizer

##
time_elapsed = @elapsed EP = generate_model(mysetup, myinputs, OPTIMIZER)
println("Time elapsed for model building is")
println(time_elapsed)

println("Solving Model")
EP, solve_time = solve_model(EP, mysetup)

println("Writing Output Debug")
outputs_path = get_default_output_folder(case)
elapsed_time = @elapsed outputs_path = write_outputs(EP,
    outputs_path,
    mysetup,
    myinputs)
println("Time elapsed for writing is")
println(elapsed_time)
if mysetup["ModelingToGenerateAlternatives"] == 1
    println("Starting Model to Generate Alternatives (MGA) Iterations")
    mga(EP, case, mysetup, myinputs)
end

if mysetup["MethodofMorris"] == 1
    println("Starting Global sensitivity analysis with Method of Morris")
    morris(EP, case, mysetup, myinputs, outputs_path, OPTIMIZER)
end

###

