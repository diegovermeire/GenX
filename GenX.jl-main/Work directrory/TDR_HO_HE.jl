using CSV
using DataFrames
# using YAML

function TDR_HO_to_HE(path::String)
    # Define filenames
    p_demand = "Demand_data.csv"
    p_fuel = "Fuels_data.csv"
    p_CF = "Generators_variability.csv"

    # Construct full file paths
    file_demand = joinpath(path, p_demand)
    file_fuel = joinpath(path, p_fuel)
    file_CF = joinpath(path, p_CF)

    # Load data
    df_demand = CSV.read(file_demand, DataFrame)
    df_fuel = CSV.read(file_fuel, DataFrame)
    df_CF = CSV.read(file_CF, DataFrame)

    # Function to calculate group number
    function get_group_number(row_number)
        
        rows_per_group = df_demand[1, :Timesteps_per_Rep_Period]

        if row_number == 1
            return 1
        else
            return div(row_number - 1, rows_per_group) + 1
        end
    end

    # Process demand data
    timesteps_per_rep_period = df_demand[1, :Timesteps_per_Rep_Period]
    df_demand[!, :omega] = [df_demand[get_group_number(i), :Sub_Weights] / timesteps_per_rep_period for i in 1:nrow(df_demand)]

    midpoint = div(nrow(df_demand), 2)
    if !(midpoint % 2 == 0)
        midpoint += 1
    end
    modified_demand = deepcopy(df_demand)
    modified_demand[!, :Rel_TimeStep] .= 1

    rows_to_keep = collect(1:midpoint)
    for i in midpoint+1:2:nrow(df_demand)
        if i + 1 <= nrow(df_demand)
            push!(rows_to_keep, i)
            modified_demand[i, :omega] += df_demand[i + 1, :omega]
            modified_demand[i, :Rel_TimeStep] = 2
        end
    end

    final_demand = modified_demand[rows_to_keep, :]
    final_demand[!, :Time_Index] = 1:nrow(final_demand)

    println("Demand Data: Original rows = ", nrow(df_demand), ", Final rows = ", nrow(final_demand))

    # Process fuel data
    midpoint = div(nrow(df_fuel)-1, 2)
    if !(midpoint % 2 == 0)
        midpoint += 2
    else   
        midpoint += 1
    end
    rows_to_keep = collect(1:midpoint)
    for i in midpoint+1:2:nrow(df_fuel)
        if i + 1 <= nrow(df_fuel)
            push!(rows_to_keep, i)
        end
    end

    final_fuel = df_fuel[rows_to_keep, :]
    final_fuel[!, :Time_Index] = 1:nrow(final_fuel)

    println("Fuel Data: Original rows = ", nrow(df_fuel), ", Final rows = ", nrow(final_fuel))

    # Process capacity factor data
    midpoint = div(nrow(df_CF), 2)
    if !(midpoint % 2 == 0)
        midpoint += 1
    end
    rows_to_keep = collect(1:midpoint)
    for i in midpoint+1:2:nrow(df_CF)
        if i + 1 <= nrow(df_CF)
            push!(rows_to_keep, i)
        end
    end

    final_CF = df_CF[rows_to_keep, :]
    final_CF[!, :Time_Index] = 1:nrow(final_CF)

    println("CF Data: Original rows = ", nrow(df_CF), ", Final rows = ", nrow(final_CF))

    # Save the final DataFrames
    CSV.write(file_demand, final_demand)
    CSV.write(file_fuel, final_fuel)
    CSV.write(file_CF, final_CF)

    println("Data processing complete. Results saved.")
end
