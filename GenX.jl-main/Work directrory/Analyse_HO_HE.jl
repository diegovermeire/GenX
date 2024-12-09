using DataFrames
using CSV
using Statistics
using Printf

"""
Load all files from a given folder into a dictionary-like structure.
For CSV files, the content is loaded as a DataFrame.
"""
function load_files(folder_path::String)
    files_data = Dict{String, Any}()
    for entry in readdir(folder_path)
        file_path = joinpath(folder_path, entry)
        if isfile(file_path)
            file_key = splitext(basename(entry))[1]
            if endswith(entry, ".csv")
                files_data[file_key] = CSV.read(file_path, DataFrame)
            else
                files_data[file_key] = read(file_path, String)
            end
        end
    end
    return files_data
end

"""
Add a row to a DataFrame.
"""
function add_row(df::DataFrame, row_type::String, total_homogeneous::Float64, total_heterogeneous::Float64)
    absolute_difference = total_homogeneous - total_heterogeneous
    percentage_difference = if total_homogeneous != 0
        (absolute_difference / total_homogeneous) * 100
    else
        0.0
    end
    
    new_row = DataFrame(
        type = [row_type],
        total_homogeneous = [total_homogeneous],
        total_heterogeneous = [total_heterogeneous],
        absolute_difference = [absolute_difference],
        percentage_decrease_HO_to_HE_percent = [percentage_difference]
    )
    
    return vcat(df, new_row)
end

"""
Main function to process and compare results from two folders and save output.
"""
function compare_results(folder_ho::String, folder_he::String, output_folder::String)
    # Load results
    homogenous_results = load_files(folder_ho)
    heterogenous_results = load_files(folder_he)
    
    # Initialize results DataFrame
    total_results = DataFrame(
        type = String[],
        total_homogeneous = Float64[],
        total_heterogeneous = Float64[],
        absolute_difference = Float64[],
        percentage_decrease_HO_to_HE_percent = Float64[]
    )
    
    # Capacity comparison
    for (row_ho, row_he) in zip(
        eachrow(homogenous_results["capacity"][:, ["Resource", "NewCap"]]),
        eachrow(heterogenous_results["capacity"][:, ["Resource", "NewCap"]])
    )
        if row_ho.NewCap > 0 || row_he.NewCap > 0
            total_results = add_row(
                total_results,
                "Installed Capacity [MW]: $(row_ho.Resource)",
                row_ho.NewCap,
                row_he.NewCap
            )
        end
    end
    
    # Costs comparison
    for (row_ho, row_he) in zip(
        eachrow(homogenous_results["costs"][:, ["Costs", "Total"]]),
        eachrow(heterogenous_results["costs"][:, ["Costs", "Total"]])
    )

        if row_ho.Total > 0 || row_he.Total > 0
            total_results = add_row(
                total_results,
                String(row_ho.Costs),
                row_ho.Total,
                row_he.Total
            )
        end
    end
    
    # Power comparison
    homogenous_annual_sum = homogenous_results["power"][homogenous_results["power"].Resource .== "AnnualSum", :]
    heterogenous_annual_sum = heterogenous_results["power"][heterogenous_results["power"].Resource .== "AnnualSum", :]
    
    for (i, j) in zip(2:size(homogenous_annual_sum, 2), 2:size(heterogenous_annual_sum, 2))
        homogenous_value = homogenous_annual_sum[1, i]
        heterogenous_value = heterogenous_annual_sum[1, j]
        
        if homogenous_value > 0 || heterogenous_value > 0
            total_results = add_row(
                total_results,
                "Power [MWh]: $(names(homogenous_annual_sum)[i])",
                homogenous_value,
                heterogenous_value
            )
        end
    end
    
    # Save results
    if !isdir(output_folder)
        mkpath(output_folder)
    end
    output_path = joinpath(output_folder, "HO_HE_comparison.csv")
    CSV.write(output_path, total_results)
    println("Results saved to $output_path")
end
