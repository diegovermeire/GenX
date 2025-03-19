@doc raw"""
    hoursbefore(t::Int, b::Int, timesteps_per_hour::Vector{Int})

    Determines the time index b timesteps before index t in
    a landscape starting from t=1, where each hour has a
    variable number of timesteps defined by timesteps_per_hour,
    and periods wrap every p timesteps.

    For example, if p = 10 and timesteps_per_hour = [2, 3, 1, 4],
    5 timesteps before t=1 wraps around within the period.
"""
function hoursbefore_index( t::Int, b::Int, timesteps_per_hour::Vector{Int})::Int
    remaining_timesteps = b
    current_timestep = t  # Start at the current timestep

    while remaining_timesteps > 0
        # Determine the current hour index (1-based, wrap around)
        hour_index = mod1(current_timestep - 1, length(timesteps_per_hour))
        
        # Get the number of timesteps in the current hour
        timesteps_in_hour = timesteps_per_hour[hour_index]

        if remaining_timesteps <= timesteps_in_hour
            # If the remaining timesteps fit within this hour, return the hour index
            return hour_index
        else
            # Subtract the timesteps for this hour
            remaining_timesteps -= timesteps_in_hour
            current_timestep -= 1
            
            # Wrap around within the period
            if current_timestep <= 0
                current_timestep = length(timesteps_per_hour)
            end
        end
    end

    return -1  # Fallback (should never reach here in proper inputs)
end


@doc raw"""
    cumulative_sum(inputs::Dict)
    Returns the cumulative sum of Rel_Timesteps 
"""
function cumulative_sum(inputs::Dict)
    total = 0
    cumulative_rel_timestep = Vector{Int}(undef, length(inputs["Rel_TimeStep"]))

    for i in eachindex(inputs["Rel_TimeStep"])
        total += inputs["Rel_TimeStep"][i]
        cumulative_rel_timestep[i] = total
    end

    return cumulative_rel_timestep
end

@doc raw"""
    indexes_in_same_period(cumul_sum::Vector{Int}, H::Int, k::Int)::Vector{Int}
    
    Determines which indexes in the cumulative sum vector (`cumul_sum`) belong to the same period
    as the given index `k`, based on a period size `H`.

    # Arguments:
    - `cumul_sum::Vector{Int}`: A vector containing the cumulative sums of the dataset.
    - `H::Int`: The period size, defining the range of cumulative sums in each period.
    - `k::Int`: The target index for which we want to find other indexes in the same period.

    # Returns:
    - `Vector{Int}`: A vector of indexes in the cumulative sum that fall within the same period as index `k`. 
"""
function indexes_in_same_period(cumul_sum::Vector{Int}, H::Int, k::Int)::Vector{Int}
    # Find the period number of the given index k
    period_number = ceil(Int, cumul_sum[k] / H)
    
    # Determine the range of cumulative sums for this period
    period_start = (period_number - 1) * H + 1
    period_end = period_number * H

    println("period_start", period_start)
    println("period_end", period_end)

    # Find all indexes in the same period
    result = [i for i in eachindex(cumul_sum) if cumul_sum[i] >= period_start && cumul_sum[i] <= period_end]
    
    return result
end

function hours_before_HE(t, HO_timesteps_backward, inputs::Dict, hours_per_subperiod)
    cumulative_sum_RelTimesteps = cumulative_sum(inputs)
    indices_cumul_sum = indexes_in_same_period(cumulative_sum_RelTimesteps, hours_per_subperiod, t)
    t_in_same_period = (t-indices_cumul_sum[1]) + 1
    t_start_same_period = indices_cumul_sum[1] - 1
    return t_start_same_period + hoursbefore_index(t_in_same_period, HO_timesteps_backward, inputs["Rel_TimeStep"][indices_cumul_sum])
end

##################################################################################################
my_dict = Dict{String, Vector{Int}}() 
my_dict["Rel_TimeStep"] = [1, 2, 2, 1, 1, 1, 2]
hours_before_HE(4, 3, my_dict, 5)

hoursbefore_index(2, 6, [2, 2, 2, 2])


my_dict = Dict{String, Vector{Int}}() 
my_dict["Rel_TimeStep"] = [1, 2, 1, 1, 1]
cumulative_sum(my_dict)

indexes_in_same_period(cumulative_sum(my_dict), 3, 2)
indexes_in_same_period(cumulative_sum(my_dict), 2, 4)[1]


hoursbefore_index( 1, 2, my_dict["Rel_TimeStep"][indexes_in_same_period(cumulative_sum(my_dict), 2, 4)])

# Example usage of the defined functions:

# 1. Calculate the index `b` timesteps before a given time `t`
println("Example 1: Calculating the hour index `b` timesteps before `t` using `hoursbefore_index`")

# Define timesteps per hour
timesteps_per_hour = [2, 3, 1, 4]

# Find the hour index 4 timesteps before index 1
result1 = hoursbefore_index(1, 4, timesteps_per_hour)
println("Timesteps per hour: $timesteps_per_hour")
println("4 timesteps before index 1 corresponds to hour index: $result1")
println("------------------------------------------------------------")

# 2. Compute the cumulative sum of a dataset
println("Example 2: Computing the cumulative sum of a dataset using `cumulative_sum`")

# Create a dictionary with relative timesteps
my_dict = Dict{String, Vector{Int}}()
my_dict["Rel_TimeStep"] = [1, 1, 2, 1, 1]  # Original dataset

# Compute the cumulative sum
cumulative_result = cumulative_sum(my_dict)
println("Original dataset (Rel_TimeStep): $(my_dict["Rel_TimeStep"])")
println("Cumulative sum of dataset: $cumulative_result")
println("------------------------------------------------------------")

# 3. Find indexes in the same period for a given index `k`
println("Example 3: Finding all indexes in the same period using `indexes_in_same_period`")

# Define period size (H)
H = 2

# Use cumulative sum to find indexes in the same period as index 1
indexes_in_period = indexes_in_same_period(cumulative_result, H, 1)
println("Period size (H): $H")
println("Indexes in the same period as index 1: $indexes_in_period")
println("------------------------------------------------------------")

# 4. Combine the functions for a complete analysis
println("Example 4: Combining functions to analyze periods and find timesteps before")

# Using the indexes from the same period, find the hour index 2 timesteps before
hoursbefore_index( 1, 2, my_dict["Rel_TimeStep"][indexes_in_same_period(cumulative_sum(my_dict), 2, 4)])
println("Indexes in the same period: $indexes_in_period")
# println("2 timesteps before index 1 (based on the same period): $hours_before_result")
println("------------------------------------------------------------")

#5.
my_dict = Dict{String, Vector{Int}}() 
my_dict["Rel_TimeStep"] = [1, 1, 1, 1, 2, 2]
hours_per_subperiod = 4
HO_timesteps_backward = 2
current_t = 3

hours_before_HE(current_t, HO_timesteps_backward, my_dict, hours_per_subperiod)

#6.
my_dict = Dict{String, Vector{Int}}() 
my_dict["Rel_TimeStep"] = [1, 1, 1, 1, 2, 2, 1, 3]
hours_per_subperiod = 4
HO_timesteps_backward = 3
current_t = 2

hours_before_HE(current_t, HO_timesteps_backward, my_dict, hours_per_subperiod)