## --- Load required packages
using Statistics, StatsBase, DelimitedFiles, SpecialFunctions
using Chron
using Plots; gr();
using ColorSchemes, Measures
using CSV, DataFrames

## --- Load required data
    cd(@__DIR__)
    filelocation = joinpath(@__DIR__, "MandlaData/UPbData")

    BayesianAges = DataFrame(CSV.File(joinpath(filelocation, "BayesianModeling/distresults_bootstrapped.csv")));

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

    samples = ("RBCT", "RBCU", "RBCV", "RBCZ", "RBCY", "RBCW", "RBFS", "RBFT", "RBEC", "RBFR", "RBEB", "RBFQ")
    sample_heights = [430, 497, 547, 596, 600, 605, 654, 655, 665, 682, 789, 818] #elevation in meters

    #Color palette for plotting
    colors = ColorSchemes.twelvebitrainbow;

## --- Height vs elevation plot (Bayesian ages)
    h1 = plot(xlabel="Age (Ma)", ylabel="Elevation (m)", xlims=(65.3,66.7), ylims=(400,850), framestyle=:box, margins=5mm,  ytickfontsize=9, xtickfontsize=9)
    for i = 1:length(SampleData.Name)
        sample_name = SampleData.Name[i]
        
        if sample_name != "RBFQ" && sample_name != "RBEC"
            plot!(h1, [BayesianAges.Age[i]], [sample_heights[i]], xerror=([BayesianAges.Age[i] .- BayesianAges."2.5% CI"[i]], [BayesianAges."97.5% CI"[i] .- BayesianAges.Age[i]]), seriestype=:scatter, color=colors[i], markerstrokecolor=colors[i], label="")
            annotate!(h1, [BayesianAges.Age[i]], [sample_heights[i]], text.(["$sample_name"], 9, :top))
        elseif sample_name == "RBEC"
            RBEC_age = SampleData.Age[i]
            RBEC_uncert = SampleData.Uncertainty[i] 
            sorted = sortperm(RBEC_age, rev=true)
            RBEC_age = RBEC_age[sorted];
            RBEC_uncert = RBEC_uncert[sorted];
            plot!(h1, [RBEC_age[end]], [sample_heights[i]], xerror=([RBEC_uncert[end]], [RBEC_uncert[end]]), seriestype=:scatter, color=colors[i], markerstrokecolor=colors[i], label="")
            annotate!(h1, [RBEC_age[end]], [sample_heights[i]], text.(["$sample_name"], 9, :top))
        elseif sample_name == "RBFQ"
            RBFQ_age = SampleData.Age[i]
            RBFQ_uncert = SampleData.Uncertainty[i] 
            sorted = sortperm(RBFQ_age, rev=true)
            RBFQ_age = RBFQ_age[sorted];
            RBFQ_uncert = RBFQ_uncert[sorted];
            plot!(h1, [RBFQ_age[end]], [sample_heights[i]], xerror=([RBFQ_uncert[end]], [RBFQ_uncert[end]]), seriestype=:scatter, color=colors[i], markerstrokecolor=colors[i], label="")
            annotate!(h1, [RBFQ_age[end]], [sample_heights[i]], text.(["$sample_name"], 9, :top))
        end
    end

    display(h1)
    savefig(h1, "Height_Age_plot_eastern_redboles.pdf")
    savefig(h1, "Height_Age_plot_eastern_redboles.svg")
   
## --- Comparison to Ar-Ar data from Shrivastava et al. 2015 (must include the systematic uncertainties)
   #Set up plot
        h2 = plot(xlabel="Age (Ma)", ylabel="Elevation (m)", framestyle=:box, margins=5mm,  ytickfontsize=9, xtickfontsize=9)

    #U-Pb red bole data with systematic uncertainties (this study)
        #Load sample information with systematic uncertainties
            BayesianAges_Systematic = DataFrame(CSV.File(joinpath(filelocation, "BayesianModeling/SystematicUncertainties_bootstrapped.csv")));
            MDA_Systematic = DataFrame(CSV.File(joinpath(filelocation, "BayesianModeling/SystematicUncertainties_MDA.csv")));

        #Plotting
            for i = 1:length(BayesianAges_Systematic.Sample)
                 sample_name = BayesianAges_Systematic.Sample[i]
        
                if sample_name != "RBFQ" && sample_name != "RBEC"
                    plot!(h2, [BayesianAges_Systematic.Age_Z[i]], [sample_heights[i]], xerror=([BayesianAges_Systematic.Age_Z_95minus[i]], [BayesianAges_Systematic.Age_Z_95plus[i]]), seriestype=:scatter, color=colors[i], markerstrokecolor=colors[i], label="")
                    annotate!(h2, [BayesianAges_Systematic.Age_Z[i]], [sample_heights[i]], text.(["$sample_name"], 9, :top))
                else
                    plot!(h2, [MDA_Systematic.Age_Z[i]], [sample_heights[i]], xerror=([MDA_Systematic.Age_Z_plus[i]], [MDA_Systematic.Age_Z_minus[i]]), seriestype=:scatter, color=colors[i], markerstrokecolor=colors[i], label="")
                    annotate!(h2, [MDA_Systematic.Age_Z[i]], [sample_heights[i]], text.(["$sample_name"], 9, :top))
                end
            end

    #Shrivastava et al. 2015 Ar-Ar data
        #Ar-Ar sample information
            Shrivastava_Ar = DataFrame(CSV.File(joinpath(@__DIR__, "Shrivastava_ArAr_Data/Shrivastava_ArAr_SystematicUncertainties.csv")));

        # Plateau ages (Ma) and 2 sigma external uncertainties for plagioclase separates in Shrivastava et al. 2015
            samples_Shrivastava = Shrivastava_Ar.Sample
            ages_Shrivastava = Shrivastava_Ar.Plateau #Plateau ages (preferred ages from Shrivastava et al. 2015, recalculated)
            uncert_Shrivastava = Shrivastava_Ar.Plateau1sigma .* 2 #2 sigma for plateau ages (with external uncertainties) recalculated 
            height_Shrivastava = [1103, 1022, 680, 675, 610] #Elevation (m) Double check that these are all correct

        #Sort by age
            t_Shrivastava = sortperm(ages_Shrivastava, rev=true)
            ages_Shrivastava = ages_Shrivastava[t_Shrivastava];
            uncert_Shrivastava = uncert_Shrivastava[t_Shrivastava];
        
        #Plotting
            for i = 1:length(samples_Shrivastava)
                sample_name = samples_Shrivastava[i]
                plot!(h2, [ages_Shrivastava[i]], [height_Shrivastava[i]], xerror=[uncert_Shrivastava[i]], seriestype=:scatter, color=:tan4, markerstrokecolor=:tan4, markershape=:rect, label="")
                annotate!(h2, [ages_Shrivastava[i]], [height_Shrivastava[i]], text.(["$sample_name"], 9, :top))
            end

            display(h2)
            savefig(h2, "Height_Age_plot_eastern_redboles_with_ArArdata.pdf")
            savefig(h2, "Height_Age_plot_eastern_redboles_with_ArArdata.svg")

## --- End of file