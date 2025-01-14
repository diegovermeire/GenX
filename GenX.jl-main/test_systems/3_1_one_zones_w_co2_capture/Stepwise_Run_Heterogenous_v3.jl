##
ENV["GENX_PRECOMPILE"] = "false"
using GenX
using Gurobi
using JuMP
include(raw"C:\Users\Diego\GenX\GenX.jl-main\Work directrory\TDR_HO_HE.jl")
include(raw"C:\Users\Diego\GenX\GenX.jl-main\Work directrory\TDR_HO_HO.jl")
include(raw"C:\Users\Diego\GenX\GenX.jl-main\Work directrory\Analyse_HO_HE.jl")

df_HO = DataFrame(
        Metric = String[],
        HO_Value = Float64[],
        )

df_HE = DataFrame(
    Metric = String[],
    HE_Value = Float64[],
    )

## Define variables
path = dirname(@__FILE__)
genx_settings = GenX.get_settings_path(path, "genx_settings.yml") # Settings YAML file path
writeoutput_settings = GenX.get_settings_path(path, "output_settings.yml") # Write-output settings YAML file path
mysetup = GenX.configure_settings(genx_settings, writeoutput_settings) # mysetup dictionary stores settings and GenX-specific parameters

case= path
optimizer = Gurobi.Optimizer
settings_path = GenX.get_settings_path(case)
outputs_folder = ["", ""]
TDRpath_homogenous = ""

for i in 1:2
    mysetup["HeterogenousTimesteps"] = i-1
    println("HeterogenousTimesteps is EQUAL TO ",i)
    ### Cluster time series inputs if necessary and if specified by the user
    if (mysetup["TimeDomainReduction"] == 1) && (mysetup["HeterogenousTimesteps"] == 0)
        mysetup["TimeDomainReductionFolder"] = "TDR_results_HO"
        TDRpath = joinpath(case, mysetup["TimeDomainReductionFolder"])
        TDRpath_homogenous = TDRpath
        system_path = joinpath(case, mysetup["SystemFolder"])
        GenX.prevent_doubled_timedomainreduction(system_path)

        #Verify whether TRD files exists at TRDpath
        if !GenX.time_domain_reduced_files_exist(TDRpath)
            println("Clustering Time Series Data (Grouped)...")
            GenX.cluster_inputs(case, settings_path, mysetup)
            TDR_HO_to_HO(TDRpath)
        else
            println("Time Series Data Already Clustered.!")
        end
    end

    ### Cluster time series inputs if necessary and if specified by the user
    if (mysetup["TimeDomainReduction"] == 1) && (mysetup["HeterogenousTimesteps"] == 1)        
        mysetup["TimeDomainReductionFolder"] = "TDR_results_HE"
        TDRpath = joinpath(case, mysetup["TimeDomainReductionFolder"])
        system_path = joinpath(case, mysetup["SystemFolder"])
        GenX.prevent_doubled_timedomainreduction(system_path)

        #Verify whether TRD files exists at TRDpath
        if !GenX.time_domain_reduced_files_exist(TDRpath)
            println("Clustering Time Series Data (Grouped)...")
            destination_path = replace(TDRpath_homogenous, "HO" => "HE")
            # Check if the source folder exists
            if isdir(TDRpath_homogenous)
                try
                    # Copy the folder
                    cp(TDRpath_homogenous, destination_path; force=true)
                    println("Folder successfully copied and renamed from '$TDRpath' to '$destination_path'.")
                catch e
                    println("An error occurred: ", e)
                end
            else
                println("The folder at '$TDRpath' does not exist.")
            end
            TDR_HO_to_HE(TDRpath)
        else
            println("Time Series Data Already Clustered.!")
        end
    end

    ## Configure solver
    println("Configuring Solver")
    OPTIMIZER = configure_solver(settings_path, optimizer)

    ## Load inputs
    println("Loading Inputs")
    myinputs = GenX.load_inputs(mysetup, case)

    ##Generate the model
    println("Generating the Optimization Model")
    time_elapsed = @elapsed EP = GenX.generate_model(mysetup, myinputs, OPTIMIZER)
    println("Time elapsed for model building is")
    println(time_elapsed)

    println("Solving Model")
    EP, solve_time = GenX.solve_model(EP, mysetup)
    myinputs["solve_time"] = GenX.solve_time # Store the model solve time in myinputs

    myinputs["T"]

    println("Writing Outputss")
    outputs_path = GenX.get_default_output_folder(case)
    outputs_path = replace(outputs_path, '\\' => '/')

    if mysetup["OverwriteResults"] == 1
        # Overwrite existing results if dir exists
        # This is the default behaviour when there is no flag, to avoid breaking existing code
        if !(isdir(path))
            mkpath(path)
        end
    else
        # Find closest unused ouput directory name and create it
        outputs_path = choose_output_dir(outputs_path)
    end

    outputs_folder[i] = outputs_path
    elapsed_time = @elapsed outputs_path = GenX.write_outputs(EP,
        outputs_path,
        mysetup,
        myinputs)
    println("Time elapsed for writing is")
    println(elapsed_time)

    CSV.write(joinpath(outputs_path, "mysetup.csv"), DataFrame(Key = collect(keys(mysetup)), Value = collect(values(mysetup))))
    CSV.write(joinpath(outputs_path, "myinputs.csv"), DataFrame(Key = collect(keys(myinputs)), Value = collect(values(myinputs))))

    # ################# DEBUGGING on NSE costs #################
    println("#############################################")
    println("Scenario $i ")
    println("#############################################")

    # T = myinputs["T"]     # Number of time steps
    # Z = myinputs["Z"]     # Number of zones
    # SEG = myinputs["SEG"] # Number of demand curtailment segments
    # omega = myinputs["omega"]
    # STOR_ALL = myinputs["STOR_ALL"]
    # gen = myinputs["RESOURCES"]
    # THERM_ALL = myinputs["THERM_ALL"]

    # if i == 1
    #     t_test = [905, 906]
    # end 

    # if i ==2
    #     t_test = [873]
    # end
    # s_test = [4]
    # c_NSE = 0
    # for s in s_test, t in t_test, z in 1:Z
    #     println("Segment:", s ," Unit_cost_nSE at time ", t, "is ", myinputs["omega"][t] * myinputs["pC_D_Curtail"][s])
    #     println("The Qt. of NSE in zone $s is ", (value(EP[:vNSE][s, t, z])))
    #     println("The max allowed curtailment in zone 4 is $(myinputs["pMax_D_Curtail"][s] * myinputs["pD"][t, z]))")
    #     println("                  ")
    #     c_NSE += value(omega[t]*EP[:eCNSE][s, t, z]) #entails omega
    # end
    # #############################################
    
    # ################# DEBUGGING on cost comparison btw thermal and battery #################
    # println("#############################################")
    # # NGCC_var_cost = value(sum(value(EP[:eCVar_out][y,t]) for y in THERM_ALL, t in t_test)) #entails omega
    # # NGCC_fuel_cost = sum(omega[t]*EP[:eCFuelOut][y, t] for y in THERM_ALL, t in t_test)

    # # println("Costs of $i scenario: ", value(cap_stor_cost) + value(var_stor_cost) + c_NSE + value(NGCC_var_cost) + value(NGCC_fuel_cost))
    # println("Total NSE cost ", c_NSE)

    # for s in s_test, t in t_test, z in 1:Z
    #     println(" NGCC variable cost at time $t ", value(sum(value(EP[:eCVar_out][y,t]) for y in THERM_ALL)))
    #     println(" NGCC fuel cost at time $t ", value(sum(omega[t]*EP[:eCFuelOut][y, t] for y in THERM_ALL)))
    #     println("                  ")
    # end
    # println("                  ")
    # println("                  ")
    # #############################################

    ################# DEBUGGING on fixed cost composotion and installed capacities #################
    
    # if i ==1
    #     tot_storage = 0
        G = myinputs["G"] # Number of resources (generators, storage, DR, and DERs)
        STOR_ALL = myinputs["STOR_ALL"] # Set of all storage resources

    #     using JuMP
    #     for y in 1:G
    #         for y in STOR_ALL
    #             println(value(EP[:eCFix][y]))
    #             push!(df_HO, ("Investment Discharge $y", value(EP[:eCFix][y])))
    #             tot_storage += value(EP[:eCFix][y])
    #         end
    #     end 
    #     ##############################################

    #     ################# DEBUGGING #################
    #     using JuMP
    #     STOR_ASYMMETRIC = myinputs["STOR_ASYMMETRIC"]
    #     for y in STOR_ASYMMETRIC
    #         println(value(EP[:eCFixCharge][y]))
    #         push!(df_HO, ("Investment Charge $y", value(EP[:eCFixCharge][y])))
    #         tot_storage += value(EP[:eCFixCharge][y])
    #     end 
    #     ##############################################

        ################# DEBUGGING #################
        STOR_ALL = myinputs["STOR_ALL"] # Set of all storage resources
        using JuMP
        for y in STOR_ALL
            println(value(EP[:eCFixEnergy][y]))
            push!(df_HO, ("Investment Energy $y", value(EP[:eCFixEnergy][y])))
            println("Inputs STORAGE_ALL", STOR_ALL)
            println("For storage $y the new energy_capacity [kWh] vCAPENERGY ", value(EP[:vCAPENERGY][y]))
            # println("vCAPCHARGE", value(EP[:vCAPCHARGE][y]))
            println("vCAPCHARGE", value(EP[:vCAP][y]))
        end 

    #     push!(df_HO, ("Total fixed cost", tot_storage))
    # end

    # if i ==2
    #     tot_storage = 0
    #     G = myinputs["G"] # Number of resources (generators, storage, DR, and DERs)
    #     STOR_ALL = myinputs["STOR_ALL"] # Set of all storage resources
    #     for y in 1:G
    #         for y in STOR_ALL
    #             println(value(EP[:eCFix][y]))
    #             push!(df_HE, ("Investment Discharge $y", value(EP[:eCFix][y])))
    #             tot_storage += value(EP[:eCFix][y])
    #         end
    #     end 
    #     ##############################################

    #     ################# DEBUGGING #################
    #     using JuMP
    #     STOR_ASYMMETRIC = myinputs["STOR_ASYMMETRIC"]
    #     for y in STOR_ASYMMETRIC
    #         println(value(EP[:eCFixCharge][y]))
    #         push!(df_HE, ("Investment Charge $y", value(EP[:eCFixCharge][y])))
    #         tot_storage += value(EP[:eCFixCharge][y])
    #     end 
    #     ##############################################

    #     ################# DEBUGGING #################
    #     STOR_ALL = myinputs["STOR_ALL"] # Set of all storage resources
    #     using JuMP
    #     for y in STOR_ALL
    #         println(value(EP[:eCFixEnergy][y]))
    #         push!(df_HE, ("Investment Energy $y", value(EP[:eCFixEnergy][y])))
    #         tot_storage += value(EP[:eCFixEnergy][y])
    #         println("Inputs STORAGE_ALL", STOR_ALL)
    #         println("For storage $y the new energy_capacity [kWh] vCAPENERGY ", value(EP[:vCAPENERGY][y]))
    #         # println("vCAPCHARGE", value(EP[:vCAPCHARGE][y]))
    #         println("vCAPCHARGE ", value(EP[:vCAP][y]))

    #     end 

    #     push!(df_HE, ("Total fixed cost", tot_storage))
    # end

    # println("omega ",sum(myinputs["omega"][t] for t in myinputs["T"]))
   
    # ##############################################
    # println("For scenario $i total omega over ", sum( myinputs["T"]), " ", sum(myinputs["omega"][t] for t in myinputs["T"]))
    # df_HE


    ################# DEBUGGING on index tao upwards within storage  #################
    START_SUBPERIODS = myinputs["START_SUBPERIODS"]
    INTERIOR_SUBPERIODS = myinputs["INTERIOR_SUBPERIODS"]
    hours_per_subperiod = myinputs["hours_per_subperiod"] #total number of hours per subperiod

    for t in START_SUBPERIODS
        if t == 1177 && i==2
            println("Index tao upwards from t 1176 should be 1260 but is ", GenX.index_toa_upwards(myinputs, 1176, hours_per_subperiod))
        elseif t == 1513 && i==1
            println("Index tao upwards from t 1512 should be 1680 but is ", 1513 + GenX.index_toa_upwards(myinputs, 1512, hours_per_subperiod))
        end
    end
    ################################################################################################

end

compare_results(outputs_folder[1], outputs_folder[2], outputs_folder[2])

################# DEBUGGING on fixed cost composotion and installed capacities #################
df_cFix = outerjoin(df_HO, df_HE, on = :Metric, makeunique = true)
df_cFix.HO_min_HE = df_cFix.HO_Value .- df_cFix.HE_Value
df_cFix
################################################################################################

df_HO[1,"HO_Value"] - df_HO[2,"HO_Value"]

