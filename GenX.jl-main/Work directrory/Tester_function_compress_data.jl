using DataFrames, OrderedCollections

# Example dictionary with keys (rep_period, eps) and DataFrame values
df_tot_eps = Dict(
    (1, 0.1) => DataFrame(A=1:3, B=rand(3)),
    (2, 0.1) => DataFrame(A=4:6, B=rand(3)),
    (3, 0.1) => DataFrame(A=7:9, B=rand(3)),
    (1, 0.2) => DataFrame(A=10:12, B=rand(3)),
    (2, 0.2) => DataFrame(A=13:15, B=rand(3))
)

df_tot_eps[1, 0.1]

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

epsilon_groups[0.1]

# Step 2: Sort keys within each epsilon group by rep_period
for eps in keys(epsilon_groups)
    epsilon_groups[eps] = sort(epsilon_groups[eps], by=x -> x[1])
end

# Step 3: Merge DataFrames for each epsilon
merged_dfs = Dict{Float64, DataFrame}()

for (eps, keys_list) in epsilon_groups
    # Concatenate DataFrames in sorted order
    merged_dfs[eps] = vcat([df_tot_eps[key] for key in keys_list]...)
end

df_tot_eps[1, 0.1]
df_tot_eps[2, 0.1]
df_tot_eps[3, 0.1]

merged_dfs[0.1]

# Step 4: Access merged DataFrames
for (eps, df) in merged_dfs
    println("Epsilon: ", eps)
    println(df)
    println("-------------")
end
