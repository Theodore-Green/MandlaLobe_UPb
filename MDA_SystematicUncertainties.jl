## --- Load required packages
    using Statistics, StatsBase, DelimitedFiles, SpecialFunctions, Distributions
    using Chron
    using CSV, DataFrames

## --- Add systematic uncertainties for the U-Pb maximum depositional ages (RBEC and RBFQ are the only ones interpretted as MDAs but do the calculation for all of them in case others want to interpret differently)
    #Load MDA sample files
        cd(@__DIR__)
        filelocation = joinpath(@__DIR__, "MandlaData/UPbData")

        SampleData = DataFrame(Name = [], Age = [], Uncertainty = [])
        RBCT = readdlm(joinpath(filelocation,"RBCT.csv"),',',header=false); push!(SampleData, ["RBCT", RBCT[:,1], RBCT[:,2]]);
        RBCU = readdlm(joinpath(filelocation,"RBCU.csv"),',',header=false); push!(SampleData, ["RBCU", RBCU[:,1], RBCU[:,2]]);
        RBCV = readdlm(joinpath(filelocation,"RBCV.csv"),',',header=false); push!(SampleData, ["RBCV", RBCV[:,1], RBCV[:,2]]);
        RBCZ = readdlm(joinpath(filelocation,"RBCZ.csv"),',',header=false); push!(SampleData, ["RBCZ", RBCZ[:,1], RBCZ[:,2]]);
        RBCY = readdlm(joinpath(filelocation,"RBCY.csv"),',',header=false); push!(SampleData, ["RBCY", RBCY[:,1], RBCY[:,2]]);
        RBCW = readdlm(joinpath(filelocation,"RBCW.csv"),',',header=false); push!(SampleData, ["RBCW", RBCW[:,1], RBCW[:,2]]);
        RBFS = readdlm(joinpath(filelocation,"RBFS.csv"),',',header=false); push!(SampleData, ["RBFS", RBFS[:,1], RBFS[:,2]]);
        RBFT = readdlm(joinpath(filelocation,"RBFT.csv"),',',header=false); push!(SampleData, ["RBFT", RBFT[:,1], RBFT[:,2]]);
        RBEC = readdlm(joinpath(filelocation,"RBEC.csv"),',',header=false); push!(SampleData, ["RBEC", RBEC[:,1], RBEC[:,2]]);
        RBFR = readdlm(joinpath(filelocation,"RBFR.csv"),',',header=false); push!(SampleData, ["RBFR", RBFR[:,1], RBFR[:,2]]);
        RBEB = readdlm(joinpath(filelocation,"RBEB.csv"),',',header=false); push!(SampleData, ["RBEB", RBEB[:,1], RBEB[:,2]]);
        RBFQ = readdlm(joinpath(filelocation,"RBFQ.csv"),',',header=false); push!(SampleData, ["RBFQ", RBFQ[:,1], RBFQ[:,2]]);
    # Tracer (ET2535) uncertainty converted from per cent to relative
        unc_tracer = 0.03/2/100
        
    # U-238 Decay constant and uncertainty, Myr^-1, Jaffey et al. 1971
        lambda238 = 1.55125e-10 * 1e6
        unc_lambda238 = 0.107/2/100 # converted from per cent to relative
        
    # Tracer and systematic uncertainites for RBEC and RBFQ youngest zircon MDA ages
        n = 1000000 #Same length as the Bayesian age distributions used for the other dates  

        age_dist_X = Array{Float64}(undef, length(SampleData.Name), n); 
    
        for i = 1:length(SampleData.Name)
            sample_name = SampleData.Name[i]
            ages = SampleData.Age[i]
            uncert = SampleData.Uncertainty[i] ./ 2 #1sigma uncertainties for calculations
                
            #Sort the data by age (youngest zircon ages)
                sorted = sortperm(ages, rev=true)
                ages = ages[sorted];
                uncert = uncert[sorted];
    
                age_dist_X[i,:] = rand(Normal(ages[end], uncert[end]), n)
        end

    # Convert ages to 206Pb/238U ratios of the distribution
        ratio_dist = exp.(age_dist_X.*lambda238) .- 1;

    # Add tracer uncertainty
        ratio_dist_tracerunc = Array{Float64}(undef, length(SampleData.Name), size(ratio_dist,2));
        for i=1:size(ratio_dist,2)
            ratio_dist_tracerunc[:,i] = ratio_dist[:,i].*(1 + unc_tracer*randn());
        end

    # Convert 206/238 ratios back to ages, in Ma
        age_dist_XY = log.(ratio_dist_tracerunc .+ 1) ./ lambda238
        
        # Add decay constant uncertainty
        age_dist_XYZ = Array{Float64}(undef,size(ratio_dist))
        for i=1:size(ratio_dist,2)
            age_dist_XYZ[:,i] = log.(ratio_dist_tracerunc[:,i] .+ 1)./(lambda238.*( 1 .+ unc_lambda238.*randn()))
        end

    # Calculate the means and 95% confidence intervals for different levels of systematic uncertainties
        age_dist_X_mean = nanmean(age_dist_X,dim=2) # Mean age
        age_dist_X_std =  nanstd(age_dist_X,dim=2) # Standard deviation

        age_dist_XY_mean = nanmean(age_dist_XY,dim=2) # Mean age
        age_dist_XY_std =  nanstd(age_dist_XY,dim=2) # Standard deviation
        
        age_dist_XYZ_mean = nanmean(age_dist_XYZ,dim=2) # Mean age
        age_dist_XYZ_std =  nanstd(age_dist_XYZ,dim=2) # Standard deviation
        
        age_X_2σ = [age_dist_X_mean age_dist_X_std.*2]
        age_XY_2σ = [age_dist_XY_mean age_dist_XY_std.*2]
        age_XYZ_2σ = [age_dist_XYZ_mean age_dist_XYZ_std.*2]

    #Save tracer and systematic uncertainties
        Uncertainties = vcat(["Sample" "Age_X" "Age_X_plus" "Age_X_minus" "Age_Y" "Age_Y_plus" "Age_Y_minus" "Age_Z" "Age_Z_plus" "Age_Z_minus"], hcat(collect(SampleData.Name), collect(age_dist_X_mean), collect(age_dist_X_std .* 2), collect(age_dist_X_std .* 2), collect(age_dist_XY_mean), collect(age_dist_XY_std .* 2), collect(age_dist_XY_std .* 2), collect(age_dist_XYZ_mean), collect(age_dist_XYZ_std .* 2), collect(age_dist_XYZ_std .* 2)))
        writedlm(joinpath(@__DIR__, "MandlaData/UPbData/BayesianModeling/SystematicUncertainties_MDA.csv"), Uncertainties, ',')

## --- End of file