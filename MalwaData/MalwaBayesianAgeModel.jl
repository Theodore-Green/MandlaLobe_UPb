## --- Load required packages
using Statistics, StatsBase, DelimitedFiles, SpecialFunctions
using Chron
using Plots; gr();
using CSV, DataFrames
using QuadGK, Dierckx, DSP, LsqFit, Distributions, ColorSchemes

## --- Bilinear Exponential Distribution function from Chron.jl
    function bilinear_exponential(x, p)
        i₀ = firstindex(p)
        d = BilinearExponential(p[i₀], p[i₀+1], abs(p[i₀+2]), abs(p[i₀+3]), abs(p[i₀+4]))
        return pdf.(d, x)
    end

## --- Sample Properties. Keep this sample order or change all the other files to match another order 

    nSamples = 7 
    smpl = ChronAgeData(nSamples)
    smpl.Name = ("RBCF", "RBBZ", "RBCI", "RBCJ2", "RBCL", "RBCN", "RBCQ")
    smpl.Height .= [237, 263, 383, 417, 447, 554, 606]
    smpl.Height_sigma .= zeros(nSamples)
    smpl.Age_Sidedness .= zeros(nSamples) # Sidedness (zeros by default: geochron constraints are two-sided). Use -1 for a maximum age and +1 for a minimum age, 0 for two-sided
    smpl.Path = "/Users/tg5176/Documents/GitHub/MandlaLobe_UPb/MalwaData/MalwaUPbData" # Where are the data files?
    smpl.inputSigmaLevel = 2; # i.e., are the data files 1-sigma or 2-sigma. Integer.
    smpl.Age_Unit = "Ma"
    smpl.Height_Unit = "m"

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
    savefig(h, joinpath(smpl.Path,"BootstrappedDistribution_MalwaPlateau.pdf"))
    display(h)
    
    ## --- Estimate the eruption age distributions for each sample  - - - - - - - -

    # Configure distribution model here
    distSteps = 1*10^6 # Number of steps to run in distribution MCMC
    distBurnin = distSteps÷10 # Number to discard? Should it be this (2019 version) or distSteps÷100 (current Chron.jl example)

    # Choose the form of the prior closure/crystallization distribution to use

    ## You might alternatively consider:
     dist = BootstrappedDistribution         # Requires running the bootstrapping lines above
    # dist = UniformDistribution              # A reasonable default
    # dist = MeltsVolcanicZirconDistribution  # A single magmatic pulse, truncated by eruption
    # dist = ExponentialDistribution          # Applicable for survivorship processes, potentially including inheritance/dispersion in Ar-Ar dates
    # dist = HalfNormalDistribution
    # dist = TriangularDistribution
    # dist = [2.0,1.0,0.0] # Triangular distribution used in Schoene et al. 2019
    # dist = [1.0, 1.0] # Uniform distribution used in Schoene et al. 2019

    # Run MCMC to estimate saturation and eruption/deposition age distributions
    @time tMinDistMetropolis(smpl,distSteps,distBurnin,dist)

    # This will save rank-order and distribution plots, and print results to a csv file -- you can find them in smpl.Path
    #distribution_Parameters = vcat(["Sample" "normalizationconstant" "location" "scl" "shp" "skew"], hcat(collect(smpl.Name), smpl.Params[1,:], smpl.Params[2,:],smpl.Params[3,:],smpl.Params[4,:],smpl.Params[5,:]))
    #writedlm(joinpath("MalwaData/MalwaUPbData/MalwaPlateau_distparameters_notstratigraphic.csv"), distribution_Parameters, ',')
    
## --- Run stratigraphic model  - - - - - - - - - - - - - - - - - - - - - - - - -

    # Configure the stratigraphic Monte Carlo model
    config = StratAgeModelConfiguration()
    # If you in doubt, you can probably leave these parameters as-is
    config.resolution = 1.0 # Same units as sample height. Smaller is slower!
    config.bounding = 0.1 # how far away do we place runaway bounds, as a fraction of total section height
    (bottom, top) = extrema(smpl.Height)
    npoints_approx = round(Int,length(bottom:config.resolution:top) * (1 + 2*config.bounding))
    config.nsteps = 15000 # Number of steps to run in distribution MCMC
    config.burnin = 10000*npoints_approx # Number to discard
    config.sieve = round(Int,npoints_approx) # Record one out of every nsieve steps
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # Run the stratigraphic MCMC model
    @time (mdl, agedist, lldist) = StratMetropolisDist(smpl, config)


## --- Plot stratigraphic model - - - - - - - - - - - - - - - - - - - - - - - -

    # Plot results (mean and 95% confidence interval for both model and data)
    hdl = plot(framestyle=:box,
        xlabel="Age ($(smpl.Age_Unit))",
        ylabel="Height ($(smpl.Height_Unit))",
        )
    plot!(hdl, [mdl.Age_025CI; reverse(mdl.Age_975CI)],[mdl.Height; reverse(mdl.Height)], fill=(round(Int,minimum(mdl.Height)),0.5,:blue), label="model") # Age-depth model CI
    plot!(hdl, mdl.Age, mdl.Height, linecolor=:blue, label="") # Center line
    t = smpl.Age_Sidedness .== 0 # Two-sided constraints (plot in black)
    any(t) && plot!(hdl, smpl.Age[t], smpl.Height[t], xerror=(smpl.Age[t]-smpl.Age_025CI[t],smpl.Age_975CI[t]-smpl.Age[t]),label="data",seriestype=:scatter,color=:black)
    t = smpl.Age_Sidedness .== 1 # Minimum ages (plot in cyan)
    any(t) && plot!(hdl, smpl.Age[t], smpl.Height[t], xerror=(smpl.Age[t]-smpl.Age_025CI[t],zeros(count(t))),label="",seriestype=:scatter,color=:cyan,msc=:cyan)
    any(t) && zip(smpl.Age[t], smpl.Age[t].+nanmean(smpl.Age_sigma[t])*4, smpl.Height[t]) .|> x-> plot!([x[1],x[2]],[x[3],x[3]], arrow=true, label="", color=:cyan)
    t = smpl.Age_Sidedness .== -1 # Maximum ages (plot in orange)
    any(t) && plot!(hdl, smpl.Age[t], smpl.Height[t], xerror=(zeros(count(t)),smpl.Age_975CI[t]-smpl.Age[t]),label="",seriestype=:scatter,color=:orange,msc=:orange)
    any(t) && zip(smpl.Age[t], smpl.Age[t].-nanmean(smpl.Age_sigma[t])*4, smpl.Height[t]) .|> x-> plot!([x[1],x[2]],[x[3],x[3]], arrow=true, label="", color=:orange)
    savefig(hdl,joinpath(smpl.Path, "AgeDepthModel_Malwa_bootstrapped.pdf"))
    display(hdl)

    ## --- output ages and uncertainties at heights wanted

    #height_wanted = 0;
    Ash_Ages = linterp1s(mdl.Height,mdl.Age, smpl.Height);
    Ash_Ages_min = linterp1s(mdl.Height,mdl.Age_025CI, smpl.Height);
    Ash_Ages_max = linterp1s(mdl.Height,mdl.Age_975CI, smpl.Height);
    Ashes = [Ash_Ages Ash_Ages_max-Ash_Ages Ash_Ages-Ash_Ages_min]
    #@printf("Ashbed age estimates: %0.3f +%0.3f/-%0.3f Ma", Ash_Ages, Ash_Ages_max-Ash_Ages, Ash_Ages-Ash_Ages_min)
    Bayezage = [smpl.Age smpl.Age_975CI-smpl.Age smpl.Age-smpl.Age_025CI]

    #plot!(hdl, Ash_Ages, smpl.Height, xerror=(Ash_Ages - Ash_Ages_min, Ash_Ages_max - Ash_Ages),label="strat data",seriestype=:scatter,color=:red)
    #savefig(hdl, "AgeDepthModel_MalwaPlateau_withstratconstraints.pdf")

## --- Parameters for age distribution of each sample from stratigraphic model outputs (age distribution at a given height) 
    StratigraphicDistributionParameters = DataFrame(Sample = [], normalizationconstant = [], location = [], scl = [], shp = [], skew = [])

    offset = round((top-bottom)*config.bounding/config.resolution)*config.resolution
    model_heights = (bottom-offset):config.resolution:(top+offset)
    active_height_t = bottom .<= model_heights .<= top
    model_heights = model_heights[active_height_t]  
    
    for i = 1:length(model_heights)
        for j = 1:length(smpl.Height)
            if model_heights[i] == smpl.Height[j]
                name = smpl.Name[j]
                l, u = nanminimum(agedist[i,:]), nanmaximum(agedist[i,:])
                    if isfinite(u-l)
                        # Fit custom many-parametric distribution function to histogram
                        binedges = range(l, u, length=101)
                        bincounts = histcounts(agedist[i,:], binedges)

                        t = bincounts.>0 # Only look at bins with one or more results
                        N = bincounts[t] ./ config.nsteps .* count(t) # Normalized number of MCMC steps per bin
                        bincenters = cntr(binedges)[t] # Vector of bin centers

                        # Initial guess for parameters
                        p = ones(5)
                        p[2] = nanmean(agedist[i,:])
                        p[3] = nanstd(agedist[i,:])

                        # Fit nonlinear model
                        fobj = curve_fit(bilinear_exponential,bincenters,N,p)
                        fobj.param[3:end] .= abs.(fobj.param[3:end]) # Ensure positive
                        μ, σ, shp = fobj.param[2:4]
                        area, e = quadgk(x->bilinear_exponential(x,fobj.param), μ-200σ/shp, μ+200σ/shp)
                        fobj.param[1] -= log(area) # Normalize
                        push!(StratigraphicDistributionParameters, [name, fobj.param[1], fobj.param[2], fobj.param[3], fobj.param[4], fobj.param[5]])
                    end
            end
        end
    end

    CSV.write(joinpath(smpl.Path, "distparameters_bootstrapped.csv"), StratigraphicDistributionParameters)

## --- End of file