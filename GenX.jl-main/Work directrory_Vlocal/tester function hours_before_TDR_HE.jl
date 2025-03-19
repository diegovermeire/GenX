include(raw"C:\Users\Diego\GenX\GenX.jl-main\Work directrory\function hours_before_TDR_HE(inputs, hou.jl")

"""
file: HE_TDR_DBSCAN_InputData\\TDR_results_HO_number_repr_period=2time_rep_period=8736_epsilon\\system_epsilon_0.98990909

t   Demand      Rep_period, START_SUBPERIODS, INTERIOR_SUBPERIODS
1	12066.07	1	1	0
2	14773.64	1	0	1
3	14466.06	1	0	1
4	14728.61	1	0	1
5	11646.2	2	1	0
6	14716.11	2	0	1
7	14131.82	2	0	1
8	14649.99	2	0	1

"""
path_TDR_HE_DBSCAN = raw"C:\Users\Diego\GenX\GenX.jl-main\example_systems\11_New_England_200928h\HE_TDR_DBSCAN_InputData\TDR_results_HO_number_repr_period=2time_rep_period=8736_epsilon\system_epsilon_0.98990909"
demand_path = joinpath(path_TDR_HE_DBSCAN, "Demand_data.csv")
fuel_path = joinpath(path_TDR_HE_DBSCAN, "Fuels_data.csv")
capacity_path = joinpath(path_TDR_HE_DBSCAN, "Generators_variability.csv")

df_demand = CSV.read(demand_path, DataFrame)

input = Dict()
input["Start_Subperiods"] = df_demand[!, "START_SUBPERIODS"]
input["Rep_Period"] = df_demand[!, "Rep_Period"]
input["Rel_TimeStep"] = df_demand[!, "Rel_TimeStep"]

last_index = findlast(x -> x == input["Rep_Period"][5], input["Rep_Period"])
"""
Testing hour_before_TDR_HE
"""

function find_first_last_index_rep_period(inputs, current_rep_period)
    indexes = []
    last_index = findlast(x -> x == current_rep_period, inputs["Rep_Period"])
    first_index = findfirst(x -> x == current_rep_period, inputs["Rep_Period"])
    push!(indexes, first_index)
    push!(indexes, last_index)

    if ((inputs["Start_Subperiods"][last_index] == 0) && (inputs["Start_Subperiods"][first_index] == 1))
        return indexes
    end
end


# Compute hours_before for hours_backward = 1
hour_before_TDR_HE(input, 1, 1)
hour_before_TDR_HE(input, 1, 2)
hour_before_TDR_HE(input, 1, 4)
hour_before_TDR_HE(input, 1, 5)
hour_before_TDR_HE(input, 1, 6)

# Compute hours_before for hours_backward = 4
hour_before_TDR_HE(input, 4, 1)
hour_before_TDR_HE(input, 4, 5)

# Compute hours_before for hours_backward = 720
hour_before_TDR_HE(input, 720, 1)
hour_before_TDR_HE(input, 720, 5)

# Compute hours_before for hours_backward = 721
hour_before_TDR_HE(input, 721, 1)
hour_before_TDR_HE(input, 721, 5)

# Compute hours_before for hours_backward = 1416
hour_before_TDR_HE(input, 1416, 3)
hour_before_TDR_HE(input, 8015, 4)
hour_before_TDR_HE(input, 8016, 4)

"""
Testing set_hourS_before_TDR_HE
"""
# Compute hours_before for hours_backward = 1
set_hourS_before_TDR_HE(input, 1, 1)
set_hourS_before_TDR_HE(input, 1, 2)
set_hourS_before_TDR_HE(input, 1, 4)
set_hourS_before_TDR_HE(input, 1, 5)
set_hourS_before_TDR_HE(input, 1, 6)

# Compute hours_before for hours_backward = 4
set_hourS_before_TDR_HE(input, 4, 1)
set_hourS_before_TDR_HE(input, 4, 5)

# Compute hours_before for hours_backward = 720
set_hourS_before_TDR_HE(input, 720, 1)
set_hourS_before_TDR_HE(input, 720, 5)

# Compute hours_before for hours_backward = 721
set_hourS_before_TDR_HE(input, 721, 1)
set_hourS_before_TDR_HE(input, 721, 5)

# Compute hours_before for hours_backward = 1416
set_hourS_before_TDR_HE(input, 1416, 3)
set_hourS_before_TDR_HE(input, 8015, 4)
set_hourS_before_TDR_HE(input, 8016, 4)

find_first_last_index_rep_period(input, 2)[2]
