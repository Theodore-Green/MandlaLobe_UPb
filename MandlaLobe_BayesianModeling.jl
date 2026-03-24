## --- Load required packages
using Statistics, StatsBase, DelimitedFiles, SpecialFunctions
using Chron
using Plots; gr();
using ColorSchemes

## --- Sample Properties. Keep this sample order or change all the other codes to match another order, since outputs will be saved in this order

#NOTE: RBEC and RBFQ are max depositional ages, so don't actually use this interpretation for them
    nSamples = 12 # The number of samples you have data for
    smpl = ChronAgeData(nSamples)
    smpl.Name =  ("RBCT", "RBCU", "RBCV", "RBCZ", "RBCY", "RBCW", "RBFS", "RBFT", "RBEC", "RBFR", "RBEB", "RBFQ")
    smpl.Height .=  [430, 497, 547, 596, 600, 605, 654, 655, 665, 682, 789, 818]
    smpl.Height_sigma .= [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
    smpl.Age_Sidedness .= zeros(nSamples) # Sidedness (zeros by default: geochron constraints are two-sided). Use -1 for a maximum age and +1 for a minimum age, 0 for two-sided
    smpl.Path = joinpath(@__DIR__, "MandlaData/UPbData") # Where are the data files?
    smpl.inputSigmaLevel = 2 # i.e., are the data files 1-sigma or 2-sigma. Integer.
    smpl.Age_Unit = "Ma" # Unit of measurement for ages and errors in the data files
    smpl.Height_Unit = "m" # Unit of measurement for Height and Height_sigma

# IMPORTANT: smpl.Height must increase with increasing stratigraphic height
# -- i.e., stratigraphically younger samples must be more positive. For this reason, it is convenient to represent depths below surface as negative numbers.

# For each sample in smpl.Name, we must have a csv file at smpl.Path which contains two columns of data, namely: Age, Age sigma where uncertainty (sigma) is absolute uncertainty.
# If you are using U-Pb data and want Pb-loss-aware eruption estimation, simply provide five columns of data instea of two, corresponding to ²⁰⁷Pb/²³⁵U, ²⁰⁷Pb/²³⁵U sigma, ²⁰⁶Pb/²³⁸U, ²⁰⁶Pb/²³⁸U sigma, correlation


## --- Bootstrap pre-eruptive distribution - - - - - - - - - - - - - - - - - - -

    # Bootstrap a KDE of the pre-eruptive (or pre-depositional) mineral age distribution using a KDE of stacked sample data from each data file
        BootstrappedDistribution = BootstrapCrystDistributionKDE(smpl)
        x = range(0,1,length=length(BootstrappedDistribution))
        h = plot(x, BootstrappedDistribution,
            label="Bootstrapped distribution",
            xlabel="Time (arbitrary units)",
            ylabel="Probability Density",
            framestyle=:box
        )
        savefig(h, joinpath(smpl.Path,"BayesianModeling/BootstrappedDistribution.pdf"))
        display(h)

## --- Estimate the eruption age distributions for each sample  - - - - - - - -

    # Configure distribution model
        distSteps = 1*10^6 # Number of steps to run in distribution MCMC
        distBurnin = distSteps÷10 # Number to discard

    # Choose the form of the prior closure/crystallization distribution to use (uncomment the chosen distribution line)
        dist = BootstrappedDistribution; dist_name = "bootstrapped"         # Requires running the bootstrapping lines above
        #dist = UniformDistribution; dist_name = "uniform"              # A reasonable default
        #dist = MeltsVolcanicZirconDistribution; dist_name = "meltsvolcanic"  # A single magmatic pulse, truncated by eruption
        #dist = ExponentialDistribution; dist_name = "exponential"          # Applicable for survivorship processes, potentially including inheritance/dispersion in Ar-Ar dates
        #dist = HalfNormalDistribution; dist_name = "halfnormal"
        #dist = TriangularDistribution; dist_name = "triangular"

    # Run MCMC to estimate saturation and eruption/deposition age distributions
        @time tMinDistMetropolis(smpl,distSteps,distBurnin,dist)
    
    #Rename results file to contain the type of prior used and put into BayesianModeling results folder
        mv(joinpath(smpl.Path, "distresults.csv"), joinpath(smpl.Path, "BayesianModeling/distresults_$dist_name.csv"), force=true)

    # This will save rank-order and distribution plots, and print results to a csv file -- you can find them in smpl.Path

    # Save distribution parameters for plotting and hiatus likelihood calculations
        distribution_Parameters = vcat(["Sample" "normalizationconstant" "location" "scl" "shp" "skew"], hcat(collect(smpl.Name), smpl.Params[1,:], smpl.Params[2,:],smpl.Params[3,:],smpl.Params[4,:],smpl.Params[5,:]))
        writedlm(joinpath(smpl.Path, "BayesianModeling/distparameters_$dist_name.csv"), distribution_Parameters, ',')
        writedlm(joinpath(@__DIR__, "MandlaData", "distparameters_$dist_name.csv"), distribution_Parameters, ',') #Save in MandlaData folder for later use and comparison to other types of priors
        
    #Save results in main MandlaData folder for easy later use and comparison between different types of priors
        cp(joinpath(@__DIR__, "MandlaData/UPbData/BayesianModeling", "distresults_$dist_name.csv"), joinpath(@__DIR__, "MandlaData", "distresults_$dist_name.csv"), force=true)
        cp(joinpath(@__DIR__, "MandlaData/UPbData/BayesianModeling", "distparameters_$dist_name.csv"), joinpath(@__DIR__, "MandlaData", "distparameters_$dist_name.csv"), force=true)


## --- Add systematic uncertainties for U-Pb data

    # Tracer (ET2535) uncertainty converted from per cent to relative
        unc_tracer = 0.03/2/100
        
    # U-238 Decay constant and uncertainty, Myr^-1, Jaffey et al. 1971
        lambda238 = 1.55125e-10 * 1e6
        unc_lambda238 = 0.107/2/100 # converted from per cent to relative
        
    # Compile Bayesian eruption age distributions for each sample into a single matrix for easier calculations
        age_dist_X = Array{Float64}(undef, length(smpl.Name), size(smpl.Age_Distribution[1],1));
        for i = 1:length(smpl.Name)
            age_dist_X[i,:] = smpl.Age_Distribution[i];
        end

    # Convert ages to 206Pb/238U ratios of the distribution
        ratio_dist = exp.(age_dist_X.*lambda238) .- 1;

    # Add tracer uncertainty
        ratio_dist_tracerunc = Array{Float64}(undef, length(smpl.Name), size(ratio_dist,2));
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
        age_dist_X_median = nanmedian(age_dist_X,dim=2) # Median age
        age_dist_X_025p = nanpctile(age_dist_X,2.5,dim=2) # 2.5th percentile
        age_dist_X_975p = nanpctile(age_dist_X,97.5,dim=2) # 97.5th percentile
        
        age_dist_XY_mean = nanmean(age_dist_XY,dim=2) # Mean age
        age_dist_XY_std =  nanstd(age_dist_XY,dim=2) # Standard deviation
        age_dist_XY_median = nanmedian(age_dist_XY,dim=2) # Median age
        age_dist_XY_025p = nanpctile(age_dist_XY,2.5,dim=2) # 2.5th percentile
        age_dist_XY_975p = nanpctile(age_dist_XY,97.5,dim=2) # 97.5th percentile
        
        age_dist_XYZ_mean = nanmean(age_dist_XYZ,dim=2) # Mean age
        age_dist_XYZ_std =  nanstd(age_dist_XYZ,dim=2) # Standard deviation
        age_dist_XYZ_median = nanmedian(age_dist_XYZ,dim=2) # Median age
        age_dist_XYZ_025p = nanpctile(age_dist_XYZ,2.5,dim=2) # 2.5th percentile
        age_dist_XYZ_975p = nanpctile(age_dist_XYZ,97.5,dim=2) # 97.5th percentile
        
        age_X_95p = [age_dist_X_mean age_dist_X_975p.-age_dist_X_mean age_dist_X_mean.-age_dist_X_025p]
        age_XY_95p = [age_dist_XY_mean age_dist_XY_975p.-age_dist_XY_mean age_dist_XY_mean.-age_dist_XY_025p]
        age_XYZ_95p = [age_dist_XYZ_mean age_dist_XYZ_975p.-age_dist_XYZ_mean age_dist_XYZ_mean.-age_dist_XYZ_025p]
    
    #Save tracer and uncertainty distributions
        Uncertainties = vcat(["Sample" "Age_X" "Age_X_95plus" "Age_X_95minus" "Age_Y" "Age_Y_95plus" "Age_Y_95minus" "Age_Z" "Age_Z_95plus" "Age_Z_95minus"], hcat(collect(smpl.Name), collect(age_dist_X_mean), collect(age_dist_X_975p.-age_dist_X_mean), collect(age_dist_X_mean.-age_dist_X_025p), collect(age_dist_XY_mean), collect(age_dist_XY_975p.-age_dist_XY_mean), collect(age_dist_XY_mean.-age_dist_XY_025p), collect(age_dist_XYZ_mean), collect(age_dist_XYZ_975p.-age_dist_XYZ_mean), collect(age_dist_XYZ_mean.-age_dist_XYZ_025p)))
        writedlm(joinpath(smpl.Path, "BayesianModeling/SystematicUncertainties_$dist_name.csv"), Uncertainties, ',')

## --- Not running the stratigraphic model for the Mandla Lobe samples because they do not follow stratigraphic superposition. Can uncomment following lines to run it just to see, but does not actually work for Mandla samples

    ## --- Run stratigraphic model  - - - - - - - - - - - - - - - - - - - - - - - - -

        # Configure the stratigraphic Monte Carlo model
        #config = StratAgeModelConfiguration()
        # If you in doubt, you can probably leave these parameters as-is
        #config.resolution = 0.5 # Same units as sample height. Smaller is slower!
        #config.bounding = 0.5 # how far away do we place runaway bounds, as a fraction of total section height
        #(bottom, top) = extrema(smpl.Height)
        #npoints_approx = round(Int,length(bottom:config.resolution:top) * (1 + 2*config.bounding))
        #config.nsteps = 15000 # Number of steps to run in distribution MCMC
        #config.burnin = 10000*npoints_approx # Number to discard
        #config.sieve = round(Int,npoints_approx) # Record one out of every nsieve steps
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        # Run the stratigraphic MCMC model
        #@time (mdl, agedist, lldist) = StratMetropolisDist(smpl, config)


    ## --- Plot stratigraphic model - - - - - - - - - - - - - - - - - - - - - - - -

        # Plot results (mean and 95% confidence interval for both model and data)
        #hdl = plot(framestyle=:box,
        #    xlabel="Age ($(smpl.Age_Unit))",
        #    ylabel="Height ($(smpl.Height_Unit))",
        #)
        #plot!(hdl, [mdl.Age_025CI; reverse(mdl.Age_975CI)],[mdl.Height; reverse(mdl.Height)], fill=(round(Int,minimum(mdl.Height)),0.5,:blue), label="model") # Age-depth model CI
        #plot!(hdl, mdl.Age, mdl.Height, linecolor=:blue, label="") # Center line
        #t = smpl.Age_Sidedness .== 0 # Two-sided constraints (plot in black)
        #any(t) && plot!(hdl, smpl.Age[t], smpl.Height[t], xerror=(smpl.Age[t]-smpl.Age_025CI[t],smpl.Age_975CI[t]-smpl.Age[t]),label="data",seriestype=:scatter,color=:black)
        #t = smpl.Age_Sidedness .== 1 # Minimum ages (plot in cyan)
        #any(t) && plot!(hdl, smpl.Age[t], smpl.Height[t], xerror=(smpl.Age[t]-smpl.Age_025CI[t],zeros(count(t))),label="",seriestype=:scatter,color=:cyan,msc=:cyan)
        #any(t) && zip(smpl.Age[t], smpl.Age[t].+nanmean(smpl.Age_sigma[t])*4, smpl.Height[t]) .|> x-> plot!([x[1],x[2]],[x[3],x[3]], arrow=true, label="", color=:cyan)
        #t = smpl.Age_Sidedness .== -1 # Maximum ages (plot in orange)
        #any(t) && plot!(hdl, smpl.Age[t], smpl.Height[t], xerror=(zeros(count(t)),smpl.Age_975CI[t]-smpl.Age[t]),label="",seriestype=:scatter,color=:orange,msc=:orange)
        #any(t) && zip(smpl.Age[t], smpl.Age[t].-nanmean(smpl.Age_sigma[t])*4, smpl.Height[t]) .|> x-> plot!([x[1],x[2]],[x[3],x[3]], arrow=true, label="", color=:orange)
        #savefig(hdl,"AgeDepthModel_Mandla.pdf")
        #display(hdl)


    ## --- Systematic uncertainties for Strat Model ages
        # Tracer (ET2535) uncertainty converted from per cent to relative
        #unc_tracer = 0.03/2/100;

        # U-238 Decay constant and uncertainty, Myr^-1
        #lambda238 = 1.55125e-10 * 1e6;
        #unc_lambda238 = 0.107/2/100; # converted from per cent to relative

        # Consider only the distribution of ages at model nodes where we have an ash bed
        #age_dist_X = Array{Float64}(length(smpl.Height),size(agedist,2));
        #for i = 1:length(smpl.Height)
        #closest_model_node = indmin(abs.(mdl.Height-smpl.Height[i]))
        #age_dist_X[i,:] = agedist[closest_model_node,:];
        #end

        # Convert ages to 206Pb/238U ratios of the distribution
        #ratio_dist = exp.(age_dist_X.*lambda238)-1;

        # Add tracer uncertainty
        #ratio_dist_tracerunc = Array{Float64}(size(ratio_dist));
        #for i=1:size(ratio_dist,2)
        #    ratio_dist_tracerunc[:,i] = ratio_dist[:,i].*(1 + unc_tracer*randn());
        #end

        # Convert 206/238 ratios back to ages, in Ma
        #age_dist_XY = log.(ratio_dist_tracerunc+1)./lambda238;

        # Add decay constant uncertainty
        #age_dist_XYZ = Array{Float64}(size(ratio_dist));
        #for i=1:size(ratio_dist,2)
        #    age_dist_XYZ[:,i] = log.(ratio_dist_tracerunc[:,i]+1)./(lambda238.*(1 + unc_lambda238.*randn()));
        #end

        # Calculate the means and 95% confidence intervals for different levels of systematic uncertainties

        #age_dist_X_mean = mean(age_dist_X,2); # Mean age
        #age_dist_X_std =  std(age_dist_X,2); # Standard deviation
        #age_dist_X_median = median(age_dist_X,2); # Median age
        #age_dist_X_025p = pctile(age_dist_X,2.5,dim=2); # 2.5th percentile
        #age_dist_X_975p = pctile(age_dist_X,97.5,dim=2); # 97.5th percentile

        #age_dist_XY_mean = mean(age_dist_XY,2); # Mean age
        #age_dist_XY_std =  std(age_dist_XY,2); # Standard deviation
        #age_dist_XY_median = median(age_dist_XY,2); # Median age
        #age_dist_XY_025p = pctile(age_dist_XY,2.5,dim=2); # 2.5th percentile
        #age_dist_XY_975p = pctile(age_dist_XY,97.5,dim=2); # 97.5th percentile

        #age_dist_XYZ_mean = mean(age_dist_XYZ,2); # Mean age
        #age_dist_XYZ_std =  std(age_dist_XYZ,2); # Standard deviation
        #age_dist_XYZ_median = median(age_dist_XYZ,2); # Median age
        #age_dist_XYZ_025p = pctile(age_dist_XYZ,2.5,dim=2); # 2.5th percentile
        #age_dist_XYZ_975p = pctile(age_dist_XYZ,97.5,dim=2); # 97.5th percentile

        #age_X_95p = [age_dist_X_mean age_dist_X_975p-age_dist_X_mean age_dist_X_mean-age_dist_X_025p];
        #age_XY_95p = [age_dist_XY_mean age_dist_XY_975p-age_dist_XY_mean age_dist_XY_mean-age_dist_XY_025p];
        #age_XYZ_95p = [age_dist_XYZ_mean age_dist_XYZ_975p-age_dist_XYZ_mean age_dist_XYZ_mean-age_dist_XYZ_025p];


## --- End of file