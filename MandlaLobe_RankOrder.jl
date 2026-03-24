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

    #Color palette for plotting
    colors = ColorSchemes.twelvebitrainbow;

## --- Set up rank order plot
    h1 = plot(ylabel="Age (Ma)", xlabel="", legend=false, xaxis=false, xticks=false, yflip=true, ylims=(64.90, 67.60), ygridalpha=0.3, size=(1100,500), margin=5mm, ytickfontsize=9)

## --- Magnetic reversals (Ogg 2021, Geomagnetic polarity timescale)
    rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h]) #Define plotting a rectangle (coordinates of the four vertices. Negative height here because y axis decreases upwards)
    plot!(h1, rectangle(3, -(68.178-66.380), 0, 68.178), opacity=1, color=:black, linecolor=:black, label=""); annotate!(h1, [0], [66.380], text.(["C30N"], 9, :top, color=:gray)) #C30n
    plot!(h1, rectangle(3, -(66.380-65.700), 0, 66.380), opacity=1, color=:white, linecolor=:black, label=""); annotate!(h1, [0], [66.380], text.(["C29R"], 9, :bottom, color=:black)) #C29r
    plot!(h1, rectangle(3, -(65.700-64.862), 0, 65.700), opacity=1, color=:black, linecolor=:black, label=""); annotate!(h1, [0], [65.700], text.(["C29N"], 9, :bottom, color=:gray)) #C29n
    plot!(h1, rectangle(3, -(64.862-64.645), 0, 64.862), opacity=1, color=:white, linecolor=:black, label=""); annotate!(h1, [0], [64.645], text.(["C28R"], 9, :top, color=:black)) #C28r
    plot!(h1, rectangle(3, -(64.645-63.537), 0, 64.645), opacity=1, color=:black, linecolor=:black, label=""); annotate!(h1, [0], [64.645], text.(["C28N"], 9, :bottom, color=:gray)) #C28n
    plot!(h1, rectangle(3, -(63.537-62.530), 0, 63.537), opacity=1, color=:white, linecolor=:black, label=""); annotate!(h1, [0], [63.537], text.(["C27R"], 9, :bottom, color=:black)) #C27r
    plot!(h1, rectangle(3, -(62.530-62.278), 0, 62.530), opacity=1, color=:black, linecolor=:black, label=""); annotate!(h1, [0], [62.530], text.(["C27N"], 9, :bottom, color=:gray)) #C27n

    #Starting x point for sample data
    x2 = 4;

## --- Cretaceous-Paleogene boundary (Schoene et al. 2019 recalculation of Clyde et al. 2016 U-Pb date: 66.016±0.050 Ma (0.099Ma for full systematic uncertainty))
    plot!(h1, [0,3], [66.016, 66.016], ribbon=0.050, seriestype=:line, color=:firebrick2, markerstrokecolor=:firebrick2, label="", linealpha=1, fillalpha=0.3, series_annotations=["",text("KPB",9,:bottom, :firebrick2)])


## --- Rank order plot in elevation order (left to right is lower to higher)
    #Desired samples (Match to those put into Bayesian modeling)
        samples = ("RBCT", "RBCU", "RBCV", "RBCZ", "RBCY", "RBCW", "RBFS", "RBFT", "RBEC", "RBFR", "RBEB", "RBFQ");
        sample_heights = [430, 497, 547, 596, 600, 605, 654, 655, 665, 682, 789, 818]; #elevation in meters
        weightedmean_numberofzircons = [2, 4, 4, 3, 3, 4, 2, 2, 1, 2, 4, 1] #Change these based on which zircons should be included in weighted mean age of the youngest grains

    #Data frames for age outputs
        Youngest_Zircon = DataFrame(Sample = [], AgeSummary = [], Age = [], Uncertainty_2σ = [])
        Weighted_Mean = DataFrame(Sample =[], AgeSummary = [], WeightedMeanAge = [],  Uncertainty_2σ = [], MSWD = [], N = [])

    #Rank order plot for the redboles 
        number_of_samples = 0
        startpoint = 0
        endpoint = 0

        for i = 1:length(SampleData.Name)

            sample_name = SampleData.Name[i]
            ages = SampleData.Age[i]
            uncert = SampleData.Uncertainty[i]
                
            #Sort the data by age for the rank order plot
                sorted = sortperm(ages, rev=true)
                ages = ages[sorted];
                uncert = uncert[sorted];

            #Youngest zircon age
                y1 = round(ages[end], digits=3)
                y2 = round(uncert[end], digits = 4)
                push!(Youngest_Zircon, [sample_name, "$y1 ± $y2", ages[end], uncert[end]])

            #Weighted mean calculation and age
                j = weightedmean_numberofzircons[i]
                if sample_name != "RBFQ" && sample_name != "RBEC"
                    weightedmean_age = wmean(ages[(end-j+1):end], (uncert[(end-j+1):end] ./ 2), corrected=false)
                    w1 = round(weightedmean_age[1], digits=3)
                    w2 = round(weightedmean_age[2] .* 2, digits=4)
                    push!(Weighted_Mean, [sample_name, "$w1 ± $w2", weightedmean_age[1], weightedmean_age[2] .* 2, weightedmean_age[3], j])
                else
                    push!(Weighted_Mean, [sample_name, "MDA", 0, 0, 0, j])
                end

            #Plotting 
                number_of_samples += length(ages)
                
                if i .== 1
                    startpoint = x2 + 1
                else
                    startpoint = 1 + endpoint
                end

                endpoint = x2 + number_of_samples
                
                if sample_name != "RBFQ" && sample_name != "RBEC"
                    plot!(h1, startpoint:endpoint, ages, yerror=uncert, seriestype=:scatter, color=colors[i], markerstrokecolor=colors[i], label="$sample_name")
                    plot!(h1, [startpoint, endpoint], BayesianAges.Age[i].*[1,1], ribbon=(BayesianAges.Age[i] .- BayesianAges."2.5% CI"[i], BayesianAges."97.5% CI"[i] .- BayesianAges.Age[i]), seriestype=:line, color=colors[i], markerstrokecolor=colors[i], label="", linealpha=1, fillalpha=0.3, series_annotations=[text("$sample_name",9,:bottom, color=colors[i]), ""])
                
                else #RBEC and RBFQ are maximum depositional age of the youngest zircon rather than Bayesian age like others, so don't plot the Bayesian age data. Both have only two Deccan-age zircons that are not overlapping, so treated as max depostional age instead of Bayesian eruption age
                    plot!(h1,startpoint:endpoint,ages,yerror=uncert,seriestype=:scatter, color=colors[i], markerstrokecolor=colors[i], label="$sample_name", xlims=(0,endpoint+1))
                    plot!(h1, [startpoint, endpoint], ages[end].*[1,1], ribbon=(uncert[end]), seriestype=:line, color=colors[i], markerstrokecolor=colors[i], label="", linealpha=1, fillalpha=0.3, series_annotations=[text("$sample_name",9,:bottom, color=colors[i]), ""])
                end
        end

    #Save figure
        display(h1)
        savefig(h1, "Rank_order_plot_eastern_redboles.pdf")
        savefig(h1, "Rank_order_plot_eastern_redboles.svg")

    #Save age interpretation files
        CSV.write("MandlaData/YoungestZirconAges.csv", Youngest_Zircon)
        CSV.write("MandlaData/WeightedMeanAges.csv", Weighted_Mean)


## --- Comparison to Ar-Ar data from Shrivastava et al. 2015 (must include the systematic uncertainties)
    #Set up plot
        h2 = plot(ylabel="Age (Ma)", xlabel="", legend=false, xaxis=false, xticks=false, yflip=true, ygridalpha=0.3, ylims=(62.3, 67.6), size=(1100,500), margin=5mm, ytickfontsize=9)

        #Magnetic reversals (Ogg 2021, Geomagnetic polarity timescale)
            rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h]) #Define plotting a rectangle (coordinates of the four vertices. Negative height here because y axis decreases upwards)
            plot!(h2, rectangle(3, -(68.178-66.380), 0, 68.178), opacity=1, color=:black, linecolor=:black, label=""); annotate!(h2, [0], [66.380], text.(["C30N"], 8, :top, color=:gray)) #C30n
            plot!(h2, rectangle(3, -(66.380-65.700), 0, 66.380), opacity=1, color=:white, linecolor=:black, label=""); annotate!(h2, [0], [66.380], text.(["C29R"], 8, :bottom, color=:black)) #C29r
            plot!(h2, rectangle(3, -(65.700-64.862), 0, 65.700), opacity=1, color=:black, linecolor=:black, label=""); annotate!(h2, [0], [65.700], text.(["C29N"], 8, :bottom, color=:gray)) #C29n
            plot!(h2, rectangle(3, -(64.862-64.645), 0, 64.862), opacity=1, color=:white, linecolor=:black, label=""); annotate!(h2, [0], [64.645], text.(["C28R"], 8, :top, color=:black)) #C28r
            plot!(h2, rectangle(3, -(64.645-63.537), 0, 64.645), opacity=1, color=:black, linecolor=:black, label=""); annotate!(h2, [0], [64.645], text.(["C28N"], 8, :bottom, color=:gray)) #C28n
            plot!(h2, rectangle(3, -(63.537-62.530), 0, 63.537), opacity=1, color=:white, linecolor=:black, label=""); annotate!(h2, [0], [63.537], text.(["C27R"], 8, :bottom, color=:black)) #C27r
            plot!(h2, rectangle(3, -(62.530-62.278), 0, 62.530), opacity=1, color=:black, linecolor=:black, label=""); annotate!(h2, [0], [62.530], text.(["C27N"], 8, :bottom, color=:gray)) #C27n

            x2 = 4; #Start rank order plots to the right hand side of the magnetic reversal bars
            
        #Cretaceous-Paleogene boundary (Schoene et al. 2019 recalculation of Clyde et al. 2016 U-Pb date: 66.016±0.050/0.69/0.099 Ma (full systematic uncertainty used here))
            plot!(h2, [0,3], [66.016, 66.016], ribbon=0.099, seriestype=:line, color=:firebrick2, markerstrokecolor=:firebrick2, label="", linealpha=1, fillalpha=0.3, series_annotations=["",text("KPB",9,:bottom, :firebrick2)])
        

    #U-Pb red bole dates and plotting (with systematic uncertainties)
        #Load systematic uncertainty datasets
            BayesianAges_Systematic = DataFrame(CSV.File(joinpath(filelocation, "BayesianModeling/SystematicUncertainties_bootstrapped.csv")));
            MDA_Systematic = DataFrame(CSV.File(joinpath(filelocation, "BayesianModeling/SystematicUncertainties_MDA.csv")));
        
        #Desired samples (Match to those put into Bayesian modeling)
            samples = ("RBCT", "RBCU", "RBCV", "RBCZ", "RBCY", "RBCW", "RBFS", "RBFT", "RBEC", "RBFR", "RBEB", "RBFQ");
            sample_heights = [430, 497, 547, 596, 600, 605, 654, 655, 665, 682, 789, 818]; #elevation in meters

        #Rank order plot for the redboles 
            number_of_samples = 0
            startpoint = 0
            endpoint = 0

            for i = 1:length(SampleData.Name)

                sample_name = SampleData.Name[i]
                ages = SampleData.Age[i]
                uncert = SampleData.Uncertainty[i]
                    
                #Sort the data by age for the rank order plot
                    sorted = sortperm(ages, rev=true)
                    ages = ages[sorted];
                    uncert = uncert[sorted];

                #Plotting 
                    number_of_samples += length(ages)
                    
                    if i .== 1
                        startpoint = x2 + 1
                    else
                        startpoint = 1 + endpoint
                    end

                    endpoint = x2 + number_of_samples
                    
                    if sample_name != "RBFQ" && sample_name != "RBEC"
                        plot!(h2, startpoint:endpoint, ages, yerror=uncert, seriestype=:scatter, color=colors[i], markerstrokecolor=colors[i], label="$sample_name")
                        plot!(h2, [startpoint, endpoint], BayesianAges_Systematic.Age_Z[i].*[1,1], ribbon=(BayesianAges_Systematic.Age_Z_95minus[i], BayesianAges_Systematic.Age_Z_95plus[i]), seriestype=:line, color=colors[i], markerstrokecolor=colors[i], label="", linealpha=1, fillalpha=0.3, series_annotations=[text("$sample_name",9,:bottom, color=colors[i]), ""])
                    
                    else #RBEC and RBFQ are maximum depositional age of the youngest zircon rather than Bayesian age like others, so don't plot the Bayesian age data. Both have only two Deccan-age zircons that are not overlapping, so treated as max depostional age instead of Bayesian eruption age
                        plot!(h2,startpoint:endpoint,ages,yerror=uncert,seriestype=:scatter, color=colors[i], markerstrokecolor=colors[i], label="$sample_name", xlims=(0,endpoint+1))
                        plot!(h2, [startpoint, endpoint], MDA_Systematic.Age_Z[i].*[1,1], ribbon=[MDA_Systematic.Age_Z_minus[i]], seriestype=:line, color=colors[i], markerstrokecolor=colors[i], label="", linealpha=1, fillalpha=0.3, series_annotations=[text("$sample_name",9,:bottom, color=colors[i]), ""])
                    end
            end
        
    #Ar-Ar dates and plotting
        #Load dataset
            Shrivastava_Ar = DataFrame(CSV.File(joinpath(@__DIR__, "Shrivastava_ArAr_Data/Shrivastava_ArAr_SystematicUncertainties.csv")));

        # Plateau ages (Ma) and 2 sigma systematic uncertainties for plagioclase separates in Shrivastava et al. 2015
            samples_Shrivastava = Shrivastava_Ar.Sample
            ages_Shrivastava = Shrivastava_Ar.Plateau #Plateau ages (preferred ages from Shrivastava et al. 2015)
            uncert_Shrivastava = Shrivastava_Ar.Plateau1sigma .* 2 #2 sigma for plateau ages (with systematic uncertainties) 
            number_of_samples +=length(ages_Shrivastava)

        #Sort by age
            t_Shrivastava = sortperm(ages_Shrivastava, rev=true)
            ages_Shrivastava = ages_Shrivastava[t_Shrivastava];
            uncert_Shrivastava = uncert_Shrivastava[t_Shrivastava];

        #Weighted mean with systematic uncertainties
            Shrivastava_Ar_WeightedMean_Systematic = DataFrame(CSV.File(joinpath(@__DIR__, "Shrivastava_ArAr_Data/Shrivastava_ArAr_WeightedMean_Systematic.csv")));

        #Plot
            startpoint = 2 + endpoint
            vline!(h2, [startpoint - 1], linestyle=:dash, color=:gray) #Line to separate the Ar-Ar and U-Pb data 
            endpoint = x2 + number_of_samples + 1 
            plot!(h2,startpoint:endpoint,ages_Shrivastava,yerror=uncert_Shrivastava,seriestype=:scatter,color=:tan4, markerstrokecolor=:tan4, markershape=:rect, label="", xlims=(0,endpoint+2))
            plot!(h2, [startpoint, endpoint], [Shrivastava_Ar_WeightedMean_Systematic.Plateau[1], Shrivastava_Ar_WeightedMean_Systematic.Plateau[1]], ribbon=(Shrivastava_Ar_WeightedMean_Systematic.Plateau1sigma[1] * 2), seriestype=:line, color=:tan4, markerstrokecolor=:tan4, label="", linealpha=0.7, fillalpha=0.3)
            
        #Label the two datasets
            plot!(h2, [startpoint, endpoint], [67, 67], linealpha=0, label="", series_annotations=[text("U-Pb data \n This study",9,:bottom), ""])
            plot!(h2, [startpoint, endpoint], [67.5, 67.5], linealpha=0, label="", series_annotations=[text("Shrivastava \n et al. 2015 \n Ar-Ar data",9,:bottom), ""])

    #Save and display figure
        display(h2)
        savefig(h2, "Rank_order_plot_eastern_redboles_with_ArArdata.pdf")
        savefig(h2, "Rank_order_plot_eastern_redboles_with_ArArdata.svg")

## --- End of file