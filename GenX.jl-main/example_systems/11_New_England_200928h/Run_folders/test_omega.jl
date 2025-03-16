inputs = Dict()
inputs["Weights"] = [10, 5]
inputs["Rel_TimeStep"] =[1,2,2,2,3]
inputs["Rep_Period"] =[1,1,1,2,2]
inputs["hours_per_subperiod"] = 5

inputs["omega"] = [inputs["Weights"][rp] * rel_t / inputs["hours_per_subperiod"] 
                   for (rp, rel_t) in zip(inputs["Rep_Period"], inputs["Rel_TimeStep"])]

