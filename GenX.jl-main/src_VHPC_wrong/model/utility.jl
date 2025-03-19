using CSV
using DataFrames

@doc raw"""
    hoursbefore(p::Int, t::Int, b::Int)

Determines the time index b hours before index t in
a landscape starting from t=1 which is separated
into distinct periods of length p.

For example, if p = 10,
1 hour before t=1 is t=10,
1 hour before t=10 is t=9
1 hour before t=11 is t=20
"""
function hoursbefore(p::Int, t::Int, b::Int)::Int
    period = div(t - 1, p)
    return period * p + mod1(t - b, p)
end

@doc raw"""
    hoursbefore(p::Int, t::Int, b::UnitRange)

This is a generalization of hoursbefore(... b::Int)
to allow for example b=1:3 to fetch a Vector{Int} of the three hours before
time index t.
"""
function hoursbefore(p::Int, t::Int, b::UnitRange{Int})::Vector{Int}
    period = div(t - 1, p)
    return period * p .+ mod1.(t .- b, p)
end

@doc raw"""
    hoursafter(p::Int, t::Int, a::Int)

Determines the time index a hours after index t in
a landscape starting from t=1 which is separated
into distinct periods of length p.

For example, if p = 10,
1 hour after t=9 is t=10,
1 hour after t=10 is t=1,
1 hour after t=11 is t=2
"""
function hoursafter(p::Int, t::Int, a::Int)::Int
    period = div(t - 1, p)
    return period * p + mod1(t + a, p)
end

@doc raw"""
    hoursafter(p::Int, t::Int, b::UnitRange)

This is a generalization of hoursafter(... b::Int)
to allow for example a=1:3 to fetch a Vector{Int} of the three hours after
time index t.
"""
function hoursafter(p::Int, t::Int, a::UnitRange{Int})::Vector{Int}
    period = div(t - 1, p)
    return period * p .+ mod1.(t .+ a, p)
end

@doc raw"""
    is_nonzero(df::DataFrame, col::Symbol)::BitVector

This function checks if a column in a dataframe is all zeros.
"""
function is_nonzero(df::DataFrame, col::Symbol)::BitVector
    convert(BitVector, df[!, col] .> 0)::BitVector
end

function is_nonzero(rs::Vector{<:AbstractResource}, col::Symbol)
    !isnothing(findfirst(r -> get(r, col, 0) ≠ 0, rs))
end

@doc raw""" 
    by_rid_res(rid::Integer, sym::Symbol, rs::Vector{<:AbstractResource})
    
    This function returns the value of the attribute `sym` for the resource given by the ID `rid`.
"""
function by_rid_res(rid::Integer, sym::Symbol, rs::Vector{<:AbstractResource})
    r = rs[findfirst(resource_id.(rs) .== rid)]
    # use getter function for attribute `sym` if exists in GenX, otherwise get the attribute directly
    f = isdefined(GenX, sym) ? getfield(GenX, sym) : x -> getproperty(x, sym)
    return f(r)
end

@doc raw""" 
    iterate_backward(inputs::Dict, current_t::Int, tao::Int)

    This function iterates backward from the current time index `current_t`, summing relative timesteps from `inputs["Rel_Timesteps"]` until the cumulative sum exceeds the specified threshold `tao`. It returns the index where the sum exceeds `tao` or 1 if the sum does not exceed the threshold at any point.
"""
function iterate_backward(inputs::Dict, current_t::Int, up_down_tao::Int)
    # Array containing all time_indexes up_down_tao from 
    time_set = []
    p = inputs["hours_per_subperiod"] #total number of hours per subperiod
    #Add current_t, in HE_timesteps start-up assumed to occur during first HO time index of th HE time index 
    push!(time_set, current_t)

    #Decrease up_down_tao to account for having already pushed current_t
    for t_back in 1:(up_down_tao-1)
        push!(time_set, hours_before_HE(current_t, t_back, inputs, p))
    end
    return unique(time_set)
end

@doc raw""" 
    index_toa_upwards(inputs::Dict, current_t::Int, tao_up::Int)

    This function defines how many time indexes correspond to tao_up hohmogenous timesteps upwards
"""
function index_toa_upwards(inputs::Dict, current_t::Int, tao_up::Int)
    tao_tot = 0.0
    upward_index = 0

    Rel_TimeStep = inputs["Rel_TimeStep"]

    # Iterate forwards from current_t
    for t in current_t:length(Rel_TimeStep)
        # Get the value from Rel_TimeStep[t]
        rel_timestep = Rel_TimeStep[t]
        tao_tot += rel_timestep

        # Adjust rel_timestep if t is the current_t (heterogeneous timestep)
        # if t == current_t
        #     rel_timestep -= 1.0
        # end

        # Check if tao_tot + rel_timestep is less than tao_up
        if tao_tot < tao_up
            # Add to tao_tot if condition is satisfied
            upward_index += 1
        else
            # Return the index where the condition was not satisfied
            return upward_index
        end
    end

    # Return the last index if tao_up is never exceeded
    return upward_index
end

@doc raw"""
    hoursbefore(p::Int, t::Int, b::Int, timesteps_per_hour::Vector{Int})

    Determines the time index b timesteps before index t in
    a landscape starting from t=1, where each hour has a
    variable number of timesteps defined by timesteps_per_hour,
    and periods wrap every p timesteps.

    For example, if p = 10 and timesteps_per_hour = [2, 3, 1, 4],
    5 timesteps before t=1 wraps around within the period.
"""
function hoursbefore_index( t::Int, b::Int, timesteps_per_hour::Vector{Float64})::Int
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

    # Find all indexes in the same period
    result = [i for i in eachindex(cumul_sum) if cumul_sum[i] >= period_start && cumul_sum[i] <= period_end]
    
    return result
end

@doc raw"""
    hours_before_HE(t, HO_timesteps_backward, inputs::Dict, hours_per_subperiod)
    
    Determines which index from the TDR temporal corresponds to HO_timesteps_backward homogenous timesteps backward  
    from current index t

    # Returns:
    - index: The index corresponding to HO_timesteps_backward homogenous timesteps backward  
    from current index t
"""
function hours_before_HE(current_index, HO_timesteps_backward, inputs::Dict, hours_per_subperiod)
    cumulative_sum_RelTimesteps = cumulative_sum(inputs)
    indices_cumul_sum = indexes_in_same_period(cumulative_sum_RelTimesteps, hours_per_subperiod, current_index)
    t_in_same_period = (current_index - indices_cumul_sum[1]) + 1
    t_start_same_period = indices_cumul_sum[1] - 1
    return t_start_same_period + hoursbefore_index(t_in_same_period, HO_timesteps_backward, float.(inputs["Rel_TimeStep"][indices_cumul_sum]))
end


@doc raw"""
This function calculates the multiplicative increase between the average value of the heterogeneous `timesteps_t_min_1` and the maximum average of the heterogeneous `timesteps_t` that meets the ramp constraints.

In this calculation:
- We assume minimal variance, meaning that each timestep within the homogeneous indices of `timesteps_t_min_1` is equal to the overall average value of hetereogenous `timesteps_t_min_1`.
- `k_homog_ramp_up` is the homogeneous ramp-up/down rate.
- `timesteps_t` is the number of timesteps behind the consecutive timeindex t over which the ramp-up occurs.

The function computes the cumulative effect of `k_homog_ramp_up` applied over `timesteps_t` and returns the average increase factor relative to the starting value.
"""
function heterogenous_ramp_small_variance(k_homog_ramp_up, timesteps_t, time_step_t_min_1)
    k_he_ramp_up_t = (sum((1 + k_homog_ramp_up)^i for i in 1:timesteps_t) / timesteps_t) - 1
    k_he_ramp_up_t_min_1 = (sum((1 + k_homog_ramp_up)^i for i in 1:time_step_t_min_1) / time_step_t_min_1) - 1

    k_he_final = max(k_he_ramp_up_t, k_he_ramp_up_t_min_1)

    # When timesteps_t or timesteps_t_min_1 are too large, k_he_final = Inf which the model cannot use
    if isinf(k_he_final)
        k_he_final = k_homog_ramp_up^max(timesteps_t, time_step_t_min_1)
    end

    return k_he_final
end



@doc raw"""
    hour_before_TDR_HE(inputs, hours_backward, t)

Determines the time index 'hours_backward' hours before index t in
within the same representative period, though this uses homogenous
TDR data that was transformed to hetereogenous data using DBSCAN

This function uses the inputs read from the Demand_data.csv

Please refer to working directory to better understand how the function works
C:\Users\Diego\GenX\GenX.jl-main\Work directrory\tester function hours_before_TDR_HE.jl

For example, if p = 10,
1 hour before t=1 is t=10,
1 hour before t=10 is t=9
1 hour before t=11 is t=20
"""
function hour_before_TDR_HE(inputs, hours_backward, t)
    last_index = findlast(x -> x == inputs["Rep_Period"][t], inputs["Rep_Period"])
    T = last_index

    if hours_backward == 1
        if inputs["Start_Subperiods"][t] == 1
            return T
        else
            return t - 1
        end
    else
        backward_iteration = 0
        current_index = t

        while backward_iteration < hours_backward
            if inputs["Start_Subperiods"][current_index] == 1
                preceding_index = T
            else
                preceding_index = current_index - 1
            end

            current_index = preceding_index
            backward_iteration += inputs["Rel_TimeStep"][current_index]
        end

        return current_index
    end
end

@doc raw"""
    set_hourS_before_TDR_HE(inputs, hours_backward, t)

Determines the SET OF time index 'hours_backward' hours before index t in
within the same representative period, though this uses homogenous
TDR data that was transformed to hetereogenous data using DBSCAN

This function uses the inputs read from the Demand_data.csv

Please refer to working directory to better understand how the function works
C:\Users\Diego\GenX\GenX.jl-main\Work directrory\tester function hours_before_TDR_HE.jl

For example, if p = 10,
3 hours before t=1 is t=[10, 9, 8]
3 hours before t=10 is t=[9, 8, 7]
4 hours before t=11 is t=[20,19,18,17]
"""
function set_hourS_before_TDR_HE(inputs, hours_backward, t)
    last_index = findlast(x -> x == inputs["Rep_Period"][t], inputs["Rep_Period"])
    T = last_index
    set_hours_before = []
    push!(set_hours_before, t)

    if hours_backward == 1
        if inputs["Start_Subperiods"][t] == 1
            push!(set_hours_before, T)
            return set_hours_before
        else
            push!(set_hours_before, t-1)
            return set_hours_before
        end
    else
        backward_iteration = 0
        current_index = t

        while backward_iteration < hours_backward
            if inputs["Start_Subperiods"][current_index] == 1
                preceding_index = T
            else
                preceding_index = current_index - 1
            end

            current_index = preceding_index
            push!(set_hours_before, current_index)
            backward_iteration += inputs["Rel_TimeStep"][current_index]
        end

        return set_hours_before
    end
end

@doc raw"""
    function find_first_last_index_rep_period(inputs, current_rep_period)

Return the first and last index of a rep_period
Used in long_duration_storage.jl
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