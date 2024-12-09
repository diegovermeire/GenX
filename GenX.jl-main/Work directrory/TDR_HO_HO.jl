using CSV
using DataFrames
using YAML

function TDR_HO_to_HO(path::String)
    """
    Process demand, fuel, and generator variability data.
    
    Args:
    - path (String): Directory path containing the CSV files.

    The function modifies the data and saves the updated CSV files.
    """

    # Define file names
    p_demand = "Demand_data.csv"
    p_fuel = "Fuels_data.csv"
    p_CF = "Generators_variability.csv"

    # Construct file paths
    path_demand = joinpath(path, p_demand)
    path_fuel = joinpath(path, p_fuel)
    path_CF = joinpath(path, p_CF)

    # Load data
    df_demand = CSV.read(path_demand, DataFrame)
    df_fuel = CSV.read(path_fuel, DataFrame)
    df_CF = CSV.read(path_CF, DataFrame)

    
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
    df_demand[!, :omega] = [
        df_demand[get_group_number(i), :Sub_Weights] / timesteps_per_rep_period
        for i in 1:nrow(df_demand)
    ]

    midpoint = div(nrow(df_demand), 2)
    midpoint = div(nrow(df_demand), 2)
    if !(midpoint % 2 == 0)
        midpoint+=1
    end
    modified_demand = deepcopy(df_demand)
    modified_demand[!, :Rel_TimeStep] .= 1

    for i in midpoint+1:2:nrow(df_demand)
        if i + 1 <= nrow(df_demand)  # Ensure there is a next row
            modified_demand[i + 1, :] .= collect(modified_demand[i, :])
        end
    end

    modified_demand[!, :Time_Index] = 1:nrow(modified_demand)

    # Process fuel data
    # Account for splitting uneven timeseries
    midpoint = div(nrow(df_fuel)-1, 2)
    if !(midpoint % 2 == 0)
        midpoint += 2
    else   
        midpoint += 1
    end

    modified_fuel = deepcopy(df_fuel)

    for i in midpoint+1:2:nrow(df_fuel)
        if i + 1 <= nrow(df_fuel)  # Ensure there is a next row
            modified_fuel[i + 1, :] .= collect(modified_fuel[i, :])
        end
    end

    modified_fuel[!, :Time_Index] = 1:nrow(modified_fuel)

    # Process capacity factor data
    midpoint = div(nrow(df_CF), 2)
    if !(midpoint % 2 == 0)
        midpoint += 2
    end
    modified_CF = deepcopy(df_CF)

    for i in midpoint+1:2:nrow(df_CF)
        if i + 1 <= nrow(df_CF)  # Ensure there is a next row
            modified_CF[i + 1, :] .= collect(modified_CF[i, :])
        end
    end

    modified_CF[!, :Time_Index] = 1:nrow(modified_CF)

    # Save the updated dataframes
    CSV.write(path_demand, modified_demand)
    CSV.write(path_fuel, modified_fuel)
    CSV.write(path_CF, modified_CF)

    println("Data processing complete. Updated files saved.")
end
