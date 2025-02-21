##################### Sample Runfile for GenX runs #####################
using GenX
using JuMP
using DataFrames
using CSV
using Gurobi

# The directory containing your settings folder and files
# settings_path = joinpath(@__DIR__, "Settings")

function get_settings_path(case::AbstractString)
    return joinpath(case, "Settings")
end

function get_settings_path(case::AbstractString, filename::AbstractString)
    return joinpath(get_settings_path(case), filename)
end

# The directory containing your input data
inputs_path = @__DIR__

## Initialize the settings (same for all of one type of run)
settings_path = GenX.get_settings_path(inputs_path)
print(settings_path)
genx_settings = GenX.get_settings_path(inputs_path, "genx_settings.yml")
mysetup = GenX.configure_settings(genx_settings)

# Load settings
TDRpath = joinpath(inputs_path, mysetup["TimeDomainReductionFolder"])
if mysetup["TimeDomainReduction"] == 1
    if !GenX.time_domain_reduced_files_exist(TDRpath)
        println("Clustering Time Series Data (Grouped)...")
        GenX.cluster_inputs(inputs_path, settings_path, mysetup)
    else
        println("Time Series Data Already Clustered.")
    end
end
# Setup logging 
# global_logger = setup_logging(mysetup)

### Load DOLPHYN
println("Loading packages")
# push!(LOAD_PATH, src_path)

# Setup time domain reduction and cluster inputs if necessary

# ### Configure solver
println("Configuring Solver")
OPTIMIZER = GenX.configure_solver(mysetup["Solver"], joinpath(inputs_path, "Settings"))

# Improve numerical stability by setting the BarHomogeneous attribute to 1
# This will slightly increase the solve time
# set_optimizer_attribute(OPTIMIZER, "Method", 1)
# set_optimizer_attribute(OPTIMIZER, "Method", 2)
GenX.set_optimizer_attribute(OPTIMIZER, "BarHomogeneous", 1)
# set_optimizer_attribute(OPTIMIZER, "Threads", 24)

println("Loading Inputs")
myinputs = GenX.load_inputs(mysetup, inputs_path)

## Set the load-based emissions cap
mysetup["CO2Cap"] = 2

## Set the number of CO2CapPeriods to the length of the scenario
mysetup["CO2CapPeriods"] = 23

## Set the scale factor
scale_factor = mysetup["ParameterScale"] == 1 ? ModelScalingFactor : 1

## The list of emissions constraints (in gCO2/kWh)
emiss_limit = 12.0


function meminfo()
    return Dict{String, Integer}(
        # "GC total" => Base.gc_total_bytes(Base.gc_num()),
        "GC live" => Base.gc_live_bytes(),
        "JIT" => Base.jit_total_bytes(),
        "Max. RSS" => Sys.maxrss()
    )
end

# function print_meminfo(meminfo_data::Dict{String, Integer})
#     for (key, value) in meminfo_data
#         @printf "%-10s: %9.3f MiB\n" key value/2^20
#     end
# end

output_path = joinpath(inputs_path, "Results") 

mkpath(output_path)

# Hard-coded to put all emissions in New Hampshire, but CO2 Cap is set to be system-wide
myinputs["dfMaxCO2Rate"][1] = emiss_limit / scale_factor ./ 1e3

# if isfile(joinpath(output_path, "costs.csv"))
#     println("Skipping Case for emiss limit = " * string(emiss_limit) * " because it already exists.")
#     continue
# end

## Generate model
println("Generating the Optimization Model")
memory = @allocated EP = GenX.generate_model(mysetup, myinputs, OPTIMIZER)

myinputs["memory"] = memory
########################
#### Add any additional constraints

## Solve model
println("Solving Model")
EP, solve_time = GenX.solve_model(EP, mysetup)
myinputs["solve_time"] = solve_time # Store the model solve time in myinputs

## Get the memory info dictionary


## Run MGA if the MGA flag is set to 1 else only save the least cost solution
println("Writing Output")

## Write only the outputs of interest: costs, capacity, and capacity_factor
GenX.write_capacity(output_path, myinputs, mysetup, EP)
println("Capacity written for greenfield")
GenX.write_capacityfactor(output_path, myinputs, mysetup, EP)
println("Capacity factor written for greenfield")
GenX.write_status(output_path, myinputs, mysetup, EP)
println("Status written for greenfield")
GenX.write_costs(output_path, myinputs, mysetup, EP)
println("Costs written for greenfield")
GenX.write_emissions(output_path, myinputs, mysetup, EP)
println("Emissions written for greenfield")
GenX.write_price(output_path, myinputs, mysetup, EP)
println("Price written for greenfield")


meminfo_dict = meminfo()

meminfo_df = DataFrame([:Key => collect(keys(meminfo_dict)), :Value => collect(values(meminfo_dict))])

## Write the memory info dataframe
meminfo_path = joinpath(output_path, "meminfo.csv")
CSV.write(meminfo_path, meminfo_df)

################################
## Add folder for brownfield data
## Add folder for brownfield results
## Add Nuclear to Generators_data.csv
## Finish making Generators_variability.csv
        