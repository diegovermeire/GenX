using CSV
using DataFrames

function hour_before_TDR_HE(inputs, hours_backward, t)
    last_index = findlast(x -> x == inputs["Rep_Period"][t], inputs["Rep_Period"])
    T = last_index

    if hours_backward == 1
        if inputs["START_SUBPERIODS"][t] == 1
            return T
        else
            return t - 1
        end
    else
        backward_iteration = 0
        current_index = t

        while backward_iteration < hours_backward
            if inputs["START_SUBPERIODS"][current_index] == 1
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

function set_hourS_before_TDR_HE(inputs, hours_backward, t)
    last_index = findlast(x -> x == inputs["Rep_Period"][t], inputs["Rep_Period"])
    T = last_index
    set_hours_before = []

    if hours_backward == 1
        if inputs["START_SUBPERIODS"][t] == 1
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
            if inputs["START_SUBPERIODS"][current_index] == 1
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

