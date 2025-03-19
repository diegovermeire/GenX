using CSV
using DataFrames
using Statistics
using Clustering
using DelimitedFiles
using Printf
using Dates
using FilePathsBase
using FilePaths
using Distances

function HO_TDR_to_HE_DBSCAN_seperate_max(epsilon_values, TDR_HO_path, path_to_save)
    # Define file names
    p_demand = "Demand_data.csv"
    p_fuel   = "Fuels_data.csv"
    p_CF     = "Generators_variability.csv"

    # Load data
    path_demand = joinpath(TDR_HO_path, p_demand)
    path_fuel   = joinpath(TDR_HO_path, p_fuel)
    path_CF     = joinpath(TDR_HO_path, p_CF)
    path_TDR_period_map = joinpath(TDR_HO_path, "Period_map.csv")

    df_demand = CSV.read(path_demand, DataFrame)
    df_fuel   = CSV.read(path_fuel, DataFrame)
    df_CF     = CSV.read(path_CF, DataFrame)

    # Get representative period and length
    rep_period = first(skipmissing(df_demand[!, "Rep_Periods"]))
    length_rep_period = first(skipmissing(df_demand[!, "Timesteps_per_Rep_Period"]))

    df_demand[!, "Rep_Period"] .= 0  # Use broadcasting to fill the entire column
    df_demand[!, "START_SUBPERIODS"] .= 0
    df_demand[!, "INTERIOR_SUBPERIODS"] .= 0

    # Assign representative period numbers to rows
    for i in 0:(Int(rep_period) - 1)
        start_index = i * Int(length_rep_period) + 1  # Adjust for Julia 1-indexing
        last_index = start_index + Int(length_rep_period) - 1
        df_demand[start_index:last_index, "Rep_Period"] .= i + 1
    end 

    # Identify starting subperiods: first occurrence of each Rep_Period
    start_subperiods_indexes = Int[]
    seen = Set{Any}()
    for (i, rp) in enumerate(df_demand[!, "Rep_Period"])
        if !(rp in seen)
            push!(start_subperiods_indexes, i)
            push!(seen, rp)
        end
    end

    # Identify interior subperiod indexes (all other row indices)
    all_indexes = collect(1:nrow(df_demand))
    interior_subperiod_indexes = setdiff(all_indexes, start_subperiods_indexes)

    df_demand[start_subperiods_indexes, "START_SUBPERIODS"] .= 1
    df_demand[interior_subperiod_indexes, "INTERIOR_SUBPERIODS"] .= 1

    # Select demand columns based on criteria
    selected_demand_columns = [col for col in names(df_demand) if ((startswith(col, "Demand") || (col in ["Time_Index", "Rep_Period", "START_SUBPERIODS", "INTERIOR_SUBPERIODS"])) && (col != "Demand_Segment"))]
    selected_fuel_columns = [col for col in names(df_fuel) if !(col in ["Time_Index", "None"])]
    selected_CF_columns   = [col for col in names(df_CF) if !(col in ["Time_Index", "None"])]


    df_demand_values = deepcopy(df_demand[:, selected_demand_columns])
    # Drop first row (index 1) for fuel and reset index by taking rows 2:end
    df_fuel_values = deepcopy(df_fuel[2:end, selected_fuel_columns])
    df_CF_values   = deepcopy(df_CF[:, selected_CF_columns])

    # Concatenate dataframes horizontally
    df_tot = hcat(df_demand_values, df_fuel_values, df_CF_values)

    # Assuming df_tot is your original DataFrame
    grouped_dfs = groupby(df_tot, :Rep_Period)
    # Dictionaries for storing results for each epsilon value

    # df_tot_eps = Dict{Any, Any, DataFrame}()
    df_tot_eps = Dict{Tuple{Int, Float64}, DataFrame}()
    df_DBSCAN_stat_eps = Dict{Tuple{Int, Float64},DataFrame}()

    for df_tot_rep_period in grouped_dfs

        rep_period_nr =  first(df_tot_rep_period[!, "Rep_Period"])

        # Standardize (normalize) df_tot values excluding specific columns
        norm_exclude = Set(["Time_Index", "Rep_Period", "START_SUBPERIODS", "INTERIOR_SUBPERIODS"])

        # norm_cols = [col for col in names(df_tot) if !(col in norm_exclude)]
        norm_cols = [col for col in names(df_tot_rep_period) if !(col in norm_exclude)]

        M = Matrix(df_tot_rep_period[:, norm_cols])
        M = float.(M)  # Convert values to Float64 if needed
        means = mean(M, dims=1)
        stds = std(M, dims=1; corrected=true)
        stds[stds .== 0] .= 1
        df_tot_normalized = (M .- means) ./ stds

        # Process each epsilon value for DBSCAN clustering
        for eps in epsilon_values
            row_indexes_to_drop = Int[]  # Will accumulate row indices to drop

            # Apply DBSCAN clustering.
            # Clustering.jl expects observations as columns, so we transpose the normalized matrix.
            dbscan_result = dbscan(df_tot_normalized', eps; metric=Euclidean(), min_neighbors=2)
            # In Clustering.jl, noise is labeled as 0. Convert to -1 for consistency.
            labels = [ (l == 0 ? -1 : l) for l in dbscan_result.assignments ]

            # Add cluster labels and initialize relative timestep
            df_tot_rep_period[!, "Cluster"] = labels
            df_tot_rep_period[!, "Rel_TimeStep"] .= 1

            # Store a copy of the dataset for this eps value
            df_tot_eps[rep_period_nr, eps] = deepcopy(df_tot_rep_period)

            # Remove added columns from the original dataframe
            select!(df_tot_rep_period, Not(["Cluster", "Rel_TimeStep"]))

            # Process each identified cluster (ignore noise with label -1)
            for label in unique(labels)
                if label != -1
                    # Obtain the indices of rows in df_tot_eps[eps] with this cluster label.
                    cluster_indices = findall(x -> x == label, df_tot_eps[rep_period_nr, eps][!, "Cluster"])
                    # Sort the indices by "Time_Index" to ensure chronological order.
                    sorted_indices = sort(cluster_indices, by = i -> df_tot_eps[rep_period_nr, eps][i, "Time_Index"])
                    
                    index_consecutive_indices_to_merge = Int[]
                    n_sorted = length(sorted_indices)
                    for i in 1:(n_sorted-1)
                        # Check if rows are consecutive and belong to the same representative period
                        current_idx = sorted_indices[i]
                        next_idx = sorted_indices[i+1]
                        if (df_tot_eps[rep_period_nr, eps][current_idx, "START_SUBPERIODS"] != 1) &&
                            (df_tot_eps[rep_period_nr, eps][current_idx, "Rep_Period"] == df_tot_eps[rep_period_nr, eps][next_idx, "Rep_Period"]) &&
                            (df_tot_eps[rep_period_nr, eps][current_idx, "Time_Index"] + 1 == df_tot_eps[rep_period_nr, eps][next_idx, "Time_Index"])
                            push!(index_consecutive_indices_to_merge, current_idx)
                            # Check if the next row is not consecutive, finalize the merge set
                            if ((i <= n_sorted - 2 && (df_tot_eps[rep_period_nr, eps][sorted_indices[i+1], "Time_Index"] + 1 != df_tot_eps[rep_period_nr, eps][sorted_indices[i+2], "Time_Index"])) ||
                                (i <= n_sorted - 2 && (df_tot_eps[rep_period_nr, eps][sorted_indices[i+1], "Rep_Period"] != df_tot_eps[rep_period_nr, eps][sorted_indices[i+2], "Rep_Period"])))
                                push!(index_consecutive_indices_to_merge, next_idx)
                                # Compute mean values for the merged rows (for numeric columns)
                                mean_row = Dict{String,Any}()
                                for col in names(df_tot_eps[rep_period_nr, eps])
                                    # Attempt to compute mean if the column element is Number
                                    vals = df_tot_eps[rep_period_nr, eps][index_consecutive_indices_to_merge, col]
                                    if all(x -> x isa Number, vals)
                                        mean_row[col] = maximum(skipmissing(vals))
                                    else
                                        # For non-numeric columns, take the first value
                                        mean_row[col] = vals[1]
                                    end
                                end
                                # Preserve representative period explicitly
                                mean_row["Rep_Period"] = df_tot_eps[rep_period_nr, eps][current_idx, "Rep_Period"]
                                # Set relative timestep to count merged rows
                                mean_row["Rel_TimeStep"] = length(index_consecutive_indices_to_merge)
                                
                                # Replace the first row with the merged values
                                index_row_to_keep = index_consecutive_indices_to_merge[1]
                                for col in keys(mean_row)
                                    # Check if the column expects an integer and convert if necessary
                                    if eltype(df_tot_eps[rep_period_nr, eps][!, col]) <: Integer
                                        df_tot_eps[rep_period_nr, eps][index_row_to_keep, col] = round(Int, mean_row[col])  # ✅ Ensure integer storage
                                    else
                                        df_tot_eps[rep_period_nr, eps][index_row_to_keep, col] = mean_row[col]
                                    end
                                end

                                # Schedule other rows for removal
                                append!(row_indexes_to_drop, index_consecutive_indices_to_merge[2:end])
                                empty!(index_consecutive_indices_to_merge)
                            end
                            # Handle the last row in the cluster
                            if i == n_sorted - 1
                                push!(index_consecutive_indices_to_merge, next_idx)
                                mean_row = Dict{String,Any}()
                                for col in names(df_tot_eps[rep_period_nr, eps])
                                    vals = df_tot_eps[rep_period_nr, eps][index_consecutive_indices_to_merge, col]
                                    if all(x -> x isa Number, vals)
                                        mean_row[col] = maximum(skipmissing(vals))
                                    else
                                        mean_row[col] = vals[1]
                                    end
                                end
                                mean_row["Rel_TimeStep"] = length(index_consecutive_indices_to_merge)
                                index_row_to_keep = index_consecutive_indices_to_merge[1]
                                for col in keys(mean_row)
                                    # Check if the column expects an integer and convert if necessary
                                    if eltype(df_tot_eps[rep_period_nr, eps][!, col]) <: Integer
                                        df_tot_eps[rep_period_nr, eps][index_row_to_keep, col] = round(Int, mean_row[col])  # ✅ Ensure integer storage
                                    else
                                        df_tot_eps[rep_period_nr, eps][index_row_to_keep, col] = mean_row[col]
                                    end
                                end
                                append!(row_indexes_to_drop, index_consecutive_indices_to_merge[2:end])
                                empty!(index_consecutive_indices_to_merge)
                            end
                        end
                    end
                end
            end
            # Remove merged rows that were replaced by averaged values
            # Ensure unique row indices are dropped and sort in descending order to avoid reindexing issues.
            row_indexes_to_drop = unique(row_indexes_to_drop)
            row_indexes_to_drop = sort(row_indexes_to_drop, rev=true)
            for idx in row_indexes_to_drop
                delete!(df_tot_eps[rep_period_nr, eps], idx)
            end

            # Store statistics about clustering results
            qt_data_points = nrow(df_tot_eps[rep_period_nr, eps])
            max_rel_timestep = maximum(df_tot_eps[rep_period_nr, eps][!, "Rel_TimeStep"])
            df_DBSCAN_stat_eps[rep_period_nr, eps] = DataFrame(Qt_Data_Points = [qt_data_points], Max_Rel_TimeStep = [max_rel_timestep])
        end
    end    

    # Step 1: Group keys by epsilon
    epsilon_groups = Dict{Float64, Vector{Tuple{Int, Float64}}}()

    for key in keys(df_tot_eps)
        rep_period[1]
        eps = key[2]
        if !haskey(epsilon_groups, eps)
            epsilon_groups[eps] = []
        end
        push!(epsilon_groups[eps], key)
    end

    # Step 2: Sort keys within each epsilon group by rep_period
    for eps in keys(epsilon_groups)
        epsilon_groups[eps] = sort(epsilon_groups[eps], by=x -> x[1])
    end


    # Step 3: Merge DataFrames for each epsilon
    df_tot_eps_merged = Dict{Float64, DataFrame}()

    for (eps, keys_list) in epsilon_groups
        # Concatenate DataFrames in sorted order
        df_tot_eps_merged[eps] = vcat([df_tot_eps[key] for key in keys_list]...)
    end

    for eps in epsilon_values
        df_tot_eps_merged[eps][!, "Time_Index"] = collect(1:nrow(df_tot_eps_merged[eps]))
    end


    # Step 3: Merge DataFrames for each epsilon
    df_tot_DBSCAN_stat_merged = Dict{Float64, DataFrame}()

    for (eps, keys_list) in epsilon_groups
        # Concatenate DataFrames in sorted order
        df_tot_DBSCAN_stat_merged[eps] = vcat([df_DBSCAN_stat_eps[key] for key in keys_list]...)
    end

    ########## New trial #########

    # println("Available keys: ", keys(df_tot_eps))
    # df_tot_eps[0.98990909][!, "START_SUBPERIODS"]
    # Append additional columns to selected_demand_columns
    push!(selected_demand_columns, "Rel_TimeStep")
    # push!(selected_demand_columns, "Rep_Period")
    # push!(selected_demand_columns, "START_SUBPERIODS")
    # push!(selected_demand_columns, "INTERIOR_SUBPERIODS")
    # Create core demand dataframe: difference between all columns and selected_demand_columns
    core_demand_cols = setdiff(names(df_demand), selected_demand_columns)
    df_demand_core = deepcopy(df_demand[:, core_demand_cols])

    push!(selected_fuel_columns, "Time_Index")
    df_fuel_core = deepcopy(df_fuel[1:1, :])

    push!(selected_CF_columns, "Time_Index")

    for eps in epsilon_values
        ## Step 1: Save modified data under the same format as the original df
        # Demand.csv
        df_demand_modified = deepcopy(dropmissing(df_tot_eps_merged[eps][:, selected_demand_columns]))
        df_demand_modified[!, :row] = 1:nrow(df_demand_modified)
        df_demand_core[!, :row] = 1:nrow(df_demand_core) 
        df_demand_eps = innerjoin(df_demand_core, df_demand_modified, makeunique=true, on="row")
        select!(df_demand_eps, Not(:row))

        # df_demand_eps = dropmissing(df_demand_eps; disallowmissing=false)

        # Fuel.csv
        df_fuel_modified = deepcopy(df_tot_eps_merged[eps][:, selected_fuel_columns])
        if "None" in names(df_fuel_modified)
            select!(df_fuel_modified, Not("None"))
        end
        if "None" in names(df_fuel_core)
            select!(df_fuel_core, Not("None"))
        end
        # Reset index by just taking the rows as they are and then adjust row numbering implicitly
        df_fuel_eps = vcat(df_fuel_core, df_fuel_modified)
        
        df_fuel_eps = dropmissing(df_fuel_eps; disallowmissing=false)

        # Generators_variability.csv
        df_CF_eps = deepcopy(dropmissing(df_tot_eps_merged[eps][:, selected_CF_columns]; disallowmissing=false))

        ## Step 2: Create a directory and save the data
        path_eps = joinpath(string(path_to_save), "system_epsilon_" * string(eps))
        mkpath(path_eps)

        # Save the DataFrame as a CSV file in the directory
        path_demand_eps = joinpath(path_eps, "Demand_data.csv")
        path_fuel_eps = joinpath(path_eps, "Fuels_data.csv")
        path_CF_eps   = joinpath(path_eps, "Generators_variability.csv")
        path_DBSCAN_stat = joinpath(path_eps, "DBSCAN_stats.csv")
        path_period_map = joinpath(path_eps, "Period_map.csv")

        println("Unique Rep_Period values in merged DataFrame for eps = ", eps, ":")
        println(unique(df_tot_eps_merged[eps][!, "Rep_Period"]))
        CSV.write(path_demand_eps, df_demand_eps)
        CSV.write(path_fuel_eps, df_fuel_eps)
        CSV.write(path_CF_eps, df_CF_eps)
        CSV.write(path_DBSCAN_stat, df_tot_DBSCAN_stat_merged[eps])
        cp(path_TDR_period_map, path_period_map, force=true)

    end

    println("DBSCAN processing complete!")
end