using CSV, DataFrames, FilePaths, FileIO
using Base.Filesystem: cp, mkpath

function create_resource_policies_enforce_cap(current_working_direct::String, path_results_to_enforce::String, source_name::String)
    # Load capacity.csv and drop the last row
    path_capacity_enforce = joinpath(path_results_to_enforce, "capacity.csv")
    capacity_to_enforce = CSV.read(path_capacity_enforce, DataFrame)
    delete!(capacity_to_enforce, nrow(capacity_to_enforce))  # Drop last row

    # Define resource paths
    folder_resource_to_copy = "resources"
    source_folder = joinpath(current_working_direct, folder_resource_to_copy)
    new_source_folder = joinpath(current_working_direct, source_name, folder_resource_to_copy * "_" * basename(path_results_to_enforce))

    # Create new directory and copy files
    if !isdir(new_source_folder)
        mkpath(new_source_folder)
        cp(source_folder, new_source_folder, force=true)
    end

    # Define data storage
    dict_tech_to_df = DataFrame(technology=String[], df_name=String[])
    dataframes_dict = Dict{String, DataFrame}()

    # Load resource data
    resource_files = Dict(
        "df_thermal" => "Thermal.csv",
        "df_hydro" => "Hydro.csv",
        "df_vre" => "Vre.csv",
        "df_electrolyzer" => "Electrolyzer.csv",
        "df_storage" => "Storage.csv",
        "df_vre_storage" => "Vre_stor.csv"
    )

    for (df_name, file_name) in resource_files
        file_path = joinpath(new_source_folder, file_name)
        if isfile(file_path)
            df = CSV.read(file_path, DataFrame)
            dataframes_dict[df_name] = df
            append!(dict_tech_to_df, DataFrame(technology=df.Resource, df_name=df_name))
        end
    end

    # Update capacity constraints
    for (df_name, df_resource) in dataframes_dict
        if "Min_Cap_MW" in names(df_resource) && "Max_Cap_MW" in names(df_resource)
            df_resource.Min_Cap_MW .= Float64.(df_resource.Min_Cap_MW)
            df_resource.Max_Cap_MW .= Float64.(df_resource.Max_Cap_MW)
        end
        if "Min_Cap_MWh" in names(df_resource) && "Max_Cap_MWh" in names(df_resource)
            df_resource.Min_Cap_MWh .= Float64.(df_resource.Min_Cap_MWh)
            df_resource.Max_Cap_MWh .= Float64.(df_resource.Max_Cap_MWh)
        end

        for row in eachrow(df_resource)
            cap_row = capacity_to_enforce[capacity_to_enforce.Resource .== row.Resource, :]
            if nrow(cap_row) > 0
                cap_MW = Float64(cap_row.EndCap[1])
                df_resource[df_resource.Resource .== row.Resource, [:Min_Cap_MW, :Max_Cap_MW]] .= cap_MW
    
                if "Max_Cap_MWh" in names(df_resource)
                    cap_MWh = Float64(cap_row.EndEnergyCap[1])
                    df_resource[df_resource.Resource .== row.Resource, [:Min_Cap_MWh, :Max_Cap_MWh]] .= cap_MWh
                end
            end
        end
    end

    # Save updated files
    for (df_name, file_name) in resource_files
        file_path = joinpath(new_source_folder, file_name)
        if haskey(dataframes_dict, df_name)
            CSV.write(file_path, dataframes_dict[df_name])
        end
    end

    return new_source_folder
end


function extract_subpath(path::String, base_folder::String)
    idx = findfirst(base_folder, path)
    if idx !== nothing
        return path[first(idx):end]  # Use `first(idx)` to get the start index
    else
        return ""  # Return empty string if base_folder is not found
    end
end
