## --- Load required packages
using Statistics, StatsBase
using Plots; gr();
using ColorSchemes, Measures
using CSV, DataFrames

## --- Load required data
    cd(@__DIR__)

    Bayesian_bootstrapped = DataFrame(CSV.File("MandlaData/distresults_bootstrapped.csv",header=true));
    Bayesian_uniform = DataFrame(CSV.File("MandlaData/distresults_uniform.csv",header=true));
    Bayesian_meltsvolcanic = DataFrame(CSV.File("MandlaData/distresults_meltsvolcanic.csv",header=true));
    Bayesian_triangular = DataFrame(CSV.File("MandlaData/distresults_triangular.csv",header=true));
    Weighted_mean = DataFrame(CSV.File("MandlaData/WeightedMeanAges.csv",header=true));
    Youngest_zircon = DataFrame(CSV.File("MandlaData/YoungestZirconAges.csv",header=true));

## --- Color palette for plotting
    colors = ColorSchemes.twelvebitrainbow;

## --- Plotting all Mandla samples
    #Note: RBEC and RBFQ are maximum depositional ages, so plot only youngest zircon interpretation for those

    h1 = plot(ylabel="Age (Ma)", xlabel="", xaxis=false, xticks=false, yflip=true, ygridalpha=0.3, margin=5mm, ytickfontsize=9, size=(1100,500), legend=:topright)

    #KPB
        hline!(h1, [66.016], ribbon=[0.050], label="", color=:firebrick2, markerstrokecolor=:firebrick2, fillalpha=0.3)
        annotate!(h1, [0], [66.016], text.(["KPB"], 9, :top, color=:firebrick2))

    #Plot age interpretations    
        startpoint = 1

        for i = 1:length(Bayesian_bootstrapped.Sample)
            sample_name = Bayesian_bootstrapped.Sample[i]
            if sample_name != "RBFQ" && sample_name != "RBEC" 
                plot!(h1, [startpoint], [Bayesian_bootstrapped.Age[i]], yerror = ([Bayesian_bootstrapped.Age[i] .- Bayesian_bootstrapped."2.5% CI"[i]], [Bayesian_bootstrapped."97.5% CI"[i] .- Bayesian_bootstrapped.Age[i]]), seriestype=:scatter, color=colors[i], markerstrokecolor=colors[i], label="", markershape=:circle, markersize=4.5, markerstrokewidth=2)
                plot!(h1, [startpoint+1], [Bayesian_uniform.Age[i]], yerror = ([Bayesian_uniform.Age[i] .- Bayesian_uniform."2.5% CI"[i]], [Bayesian_uniform."97.5% CI"[i] .- Bayesian_uniform.Age[i]]), seriestype=:scatter, color=colors[i], markerstrokecolor=colors[i], label="", markershape=:square, markersize=4.5, markerstrokewidth=2)
                plot!(h1, [startpoint+2], [Bayesian_meltsvolcanic.Age[i]], yerror = ([Bayesian_meltsvolcanic.Age[i] .- Bayesian_meltsvolcanic."2.5% CI"[i]], [Bayesian_meltsvolcanic."97.5% CI"[i] .- Bayesian_meltsvolcanic.Age[i]]), seriestype=:scatter, color=colors[i], markerstrokecolor=colors[i], label="", markershape=:utriangle, markersize=4.5, markerstrokewidth=2)
                plot!(h1, [startpoint+3], [Bayesian_triangular.Age[i]], yerror = ([Bayesian_triangular.Age[i] .- Bayesian_triangular."2.5% CI"[i]], [Bayesian_triangular."97.5% CI"[i] .- Bayesian_triangular.Age[i]]), seriestype=:scatter, color=colors[i], markerstrokecolor=colors[i], label="", markershape=:diamond, markersize=4.5, markerstrokewidth=2)
                plot!(h1, [startpoint+4], [Weighted_mean.WeightedMeanAge[i]], yerror = (Weighted_mean.Uncertainty_2σ[i]), seriestype=:scatter, color=colors[i], markerstrokecolor=colors[i], label="", markershape=:xcross, markersize=4.5, markerstrokewidth=2)
                plot!(h1, [startpoint+5], [Youngest_zircon.Age[i]], yerror = (Youngest_zircon.Uncertainty_2σ[i]), seriestype=:scatter, color=colors[i], markerstrokecolor=colors[i], label="", markershape=:hline, markersize=4.5, markerstrokewidth=2)
                annotate!(h1, [startpoint], [Bayesian_bootstrapped."97.5% CI"[i]], text.(["$sample_name"], 10, :top, color=colors[i]))
                startpoint += 7
            else
                plot!(h1, [startpoint], [Youngest_zircon.Age[i]], yerror = (Youngest_zircon.Uncertainty_2σ[i]), seriestype=:scatter, color=colors[i], markerstrokecolor=colors[i], label="", markershape=:hline, markersize=4.5, markerstrokewidth=2)
                annotate!(h1, [startpoint], [Youngest_zircon.Age[i] + Youngest_zircon.Uncertainty_2σ[i]], text.(["$sample_name"], 10, :top, color=colors[i]))
                startpoint += 2
            end
        end

    #Legend
        plot!(h1, [NaN], [NaN], label="Bayesian bootstrapped prior", seriestype=:scatter, color=:black, markershape=:circle);
        plot!(h1, [NaN], [NaN], label="Bayesian uniform prior", seriestype=:scatter, color=:black, markershape=:square);
        plot!(h1, [NaN], [NaN], label="Bayesian MELTS prior", seriestype=:scatter, color=:black, markershape=:utriangle);
        plot!(h1, [NaN], [NaN], label="Bayesian triangular prior", seriestype=:scatter, color=:black, markershape=:diamond);
        plot!(h1, [NaN], [NaN], label="Low-N weighted mean", seriestype=:scatter, color=:black, markershape=:xcross);
        plot!(h1, [NaN], [NaN], label="Youngest zircon", seriestype=:scatter, color=:black, markershape=:hline);

    #Save and display figure
        display(h1)
        savefig(h1, "AgeInterpretations_All.pdf")
        savefig(h1, "AgeInterpretations_All.svg")

## --- Plotting individual samples
    #Note: RBEC and RBFQ are maximum depositional ages, so plot only youngest zircon interpretation for those   

    for i = 1:length(Bayesian_bootstrapped.Sample)
        sample_name = Bayesian_bootstrapped.Sample[i]
        
        h2 = plot(title="", ylabel="Age (Ma)", xlabel="", xlims=(0.5,6.5), xaxis=false, xticks=false, yflip=true, margin=5mm, ytickfontsize=18, ylabelfontsize=18, background_color = :transparent)

        #KPB
            #hline!(h2, [66.016], ribbon=[0.050], label="", color=:firebrick2, markerstrokecolor=:firebrick2, fillalpha=0.3)
            #annotate!(h2, [0], [66.016], text.(["KPB"], 9, :top, color=:firebrick2))

        #Plot samples indiviudally 
            startpoint = 1
    
            if sample_name != "RBFQ" && sample_name != "RBEC" 
                plot!(h2, [startpoint], [Bayesian_bootstrapped.Age[i]], yerror = ([Bayesian_bootstrapped.Age[i] .- Bayesian_bootstrapped."2.5% CI"[i]], [Bayesian_bootstrapped."97.5% CI"[i] .- Bayesian_bootstrapped.Age[i]]), seriestype=:scatter, color=colors[i], markerstrokecolor=colors[i], label="", markershape=:circle, markersize=12)
                plot!(h2, [startpoint+1], [Bayesian_uniform.Age[i]], yerror = ([Bayesian_uniform.Age[i] .- Bayesian_uniform."2.5% CI"[i]], [Bayesian_uniform."97.5% CI"[i] .- Bayesian_uniform.Age[i]]), seriestype=:scatter, color=colors[i], markerstrokecolor=colors[i], label="", markershape=:square, markersize=12)
                plot!(h2, [startpoint+2], [Bayesian_meltsvolcanic.Age[i]], yerror = ([Bayesian_meltsvolcanic.Age[i] .- Bayesian_meltsvolcanic."2.5% CI"[i]], [Bayesian_meltsvolcanic."97.5% CI"[i] .- Bayesian_meltsvolcanic.Age[i]]), seriestype=:scatter, color=colors[i], markerstrokecolor=colors[i], label="", markershape=:utriangle, markersize=12)
                plot!(h2, [startpoint+3], [Bayesian_triangular.Age[i]], yerror = ([Bayesian_triangular.Age[i] .- Bayesian_triangular."2.5% CI"[i]], [Bayesian_triangular."97.5% CI"[i] .- Bayesian_triangular.Age[i]]), seriestype=:scatter, color=colors[i], markerstrokecolor=colors[i], label="", markershape=:diamond, markersize=12)
                plot!(h2, [startpoint+4], [Weighted_mean.WeightedMeanAge[i]], yerror = (Weighted_mean.Uncertainty_2σ[i]), seriestype=:scatter, color=colors[i], markerstrokecolor=colors[i], label="", markershape=:xcross, markersize=12)
                plot!(h2, [startpoint+5], [Youngest_zircon.Age[i]], yerror = (Youngest_zircon.Uncertainty_2σ[i]), seriestype=:scatter, color=colors[i], markerstrokecolor=colors[i], label="", markershape=:hline, markersize=12)
                
                #Legend
                    plot!(h2, [NaN], [NaN], label="Bayesian bootstrapped prior", seriestype=:scatter, color=:black, markershape=:circle);
                    plot!(h2, [NaN], [NaN], label="Bayesian uniform prior", seriestype=:scatter, color=:black, markershape=:square);
                    plot!(h2, [NaN], [NaN], label="Bayesian MELTS prior", seriestype=:scatter, color=:black, markershape=:utriangle);
                    plot!(h2, [NaN], [NaN], label="Bayesian triangular prior", seriestype=:scatter, color=:black, markershape=:diamond);
                    plot!(h2, [NaN], [NaN], label="Low-N weighted mean", seriestype=:scatter, color=:black, markershape=:xcross);
                    plot!(h2, [NaN], [NaN], label="Youngest zircon", seriestype=:scatter, color=:black, markershape=:hline);
                    
            else
                plot!(h2, [startpoint], [Youngest_zircon.Age[i]], yerror = (Youngest_zircon.Uncertainty_2σ[i]), seriestype=:scatter, color=colors[i], markerstrokecolor=colors[i], label="", markershape=:hline, markersize=12)
                #Legend
                    plot!(h2, [NaN], [NaN], label="MDA: Youngest zircon", seriestype=:scatter, color=:black, markershape=:hline);
            end
        
        #Save and display figure
            display(h2)
            savefig(h2, "MandlaData/AgeInterpretations_$sample_name.pdf")
    end


## --- End of file
