## --- Load required packages
using Statistics, StatsBase, DelimitedFiles, SpecialFunctions
using Chron
using Plots; gr();
using ColorSchemes, Measures
using CSV, DataFrames

## --- Malwa Plateau ages (Eddy et al. 2020). Bayesian MCMC stratigraphic age model outputs, bootstrapped prior, 95% CI with internal uncertainties
    Malwa_samples = ["RBCF", "RBBZ", "RBCI", "RBCJ2", "RBCL", "RBCN", "RBCQ"]
    Malwa_eruptive_ages = [66.358, 66.320, 66.275, 66.256, 66.250, 66.199, 66.150]
    Malwa_95CI_p = [0.065, 0.057, 0.066, 0.034, 0.039, 0.055, 0.065]
    Malwa_95CI_m = [0.070, 0.057, 0.086, 0.055, 0.059, 0.074, 0.084]

    t_Malwa = sortperm(Malwa_eruptive_ages, rev=true)
    Malwa_eruptive_ages = Malwa_eruptive_ages[t_Malwa]
    Malwa_95CI_p = Malwa_95CI_p[t_Malwa]
    Malwa_95CI_m = Malwa_95CI_m[t_Malwa]

## --- Western Ghats ages (Schoene et al. 2019). Bayesian MCMC stratigraphic age model outputs, triangular prior, 95% CI with internal uncertainties
    WesternGhats_samples = ["DEC13-30", "RBAB", "RBAG", "RBBQ", "RBBL", "RBBI", "RBBM", "RBBR", "RBBJ", "RBBS", "RBAW", "RBX", "RBBH", "RBBF", "RBB2", "RBAY", "RBAO", "RBAN", "RBP", "RBO", "RBE", "RBF", "DEC13-09", "RBG"]
    WesternGhats_eruptive_ages = [66.296, 66.225, 66.185, 66.161, 66.132, 66.088, 66.070, 66.055, 66.052, 66.047, 66.044, 66.039, 65.926, 65.920, 65.905, 65.895, 65.889, 65.885, 65.879, 65.875, 65.631, 65.620, 65.614, 65.590]
    WesternGhats_95CI_p = [0.037, 0.077, 0.061, 0.066, 0.069, 0.032, 0.031, 0.017, 0.015, 0.017, 0.019, 0.022, 0.035, 0.032, 0.033, 0.032, 0.029, 0.027, 0.024, 0.022, 0.053, 0.028, 0.015, 0.026]
    WesternGhats_95CI_m = [0.030, 0.071, 0.056, 0.069, 0.058, 0.026, 0.028, 0.018, 0.018, 0.020, 0.022, 0.030, 0.028, 0.031, 0.033, 0.030, 0.028, 0.029, 0.029, 0.042, 0.030, 0.021, 0.017, 0.027]

    t_WesternGhats = sortperm(WesternGhats_eruptive_ages, rev=true)
    WesternGhats_eruptive_ages = WesternGhats_eruptive_ages[t_WesternGhats]
    WesternGhats_95CI_p = WesternGhats_95CI_p[t_WesternGhats]
    WesternGhats_95CI_m = WesternGhats_95CI_m[t_WesternGhats]

## --- Mandla Lobe ages. Bayesian eruption ages, bootstrapped prior (not stratigraphic age model), 95% CI with internal uncertainties. RBEC and RBFQ are max depositonal ages, single youngest zircon with 2σ ucnertainties (they are treated separately so the plus and minus uncertainties are correctly symmetrical, but put into the same data frame as the others for convenience)
    cd(@__DIR__)
    filelocation = joinpath(@__DIR__, "MandlaData/UPbData")

    Mandla_sample_data = DataFrame(CSV.File(joinpath(filelocation, "BayesianModeling/distresults_bootstrapped.csv")))
    RBEC = DataFrame(CSV.File(joinpath(filelocation, "RBEC.csv"),header=false)); sort!(RBEC, order(:Column1,rev=true));
    RBFQ = DataFrame(CSV.File(joinpath(filelocation, "RBFQ.csv"),header=false)); sort!(RBFQ, order(:Column1,rev=true));

    Mandla_samples = ["RBCT", "RBCU", "RBCV", "RBCZ", "RBCY", "RBCW", "RBFS", "RBFT", "RBEC", "RBFR", "RBEB", "RBFQ"]
    Mandla_eruptive_ages = []
    Mandla_95CI_p = []
    Mandla_95CI_m = []

    for i = 1:length(Mandla_samples)
        if Mandla_samples[i] == "RBEC" #Max depositional age, not bayesian. Single youngest zircon, 2sigma uncertainties
            push!(Mandla_eruptive_ages, RBEC[end,1])
            push!(Mandla_95CI_p, RBEC[end,2])
            push!(Mandla_95CI_m, RBEC[end,2])

        elseif Mandla_samples[i] == "RBFQ" #Max depositional age, not bayesian. Single youngest zircon, 2sigma uncertainties
            push!(Mandla_eruptive_ages, RBFQ[end,1])
            push!(Mandla_95CI_p, RBFQ[end,2])
            push!(Mandla_95CI_m, RBFQ[end,2])

        else      #Bayesian eruption ages, 95% confidence interval 
            push!(Mandla_eruptive_ages, Mandla_sample_data.Age[i])
            push!(Mandla_95CI_p, Mandla_sample_data."97.5% CI"[i] .- Mandla_sample_data.Age[i])
            push!(Mandla_95CI_m, Mandla_sample_data.Age[i] .- Mandla_sample_data."2.5% CI"[i])

        end
    end

    #Organize in age order for rank order plot
        t_Mandla = sortperm(Mandla_eruptive_ages, rev=true)
        Mandla_eruptive_ages = Mandla_eruptive_ages[t_Mandla];
        Mandla_95CI_p = Mandla_95CI_p[t_Mandla];
        Mandla_95CI_m = Mandla_95CI_m[t_Mandla];
        Mandla_samples = Mandla_samples[t_Mandla];

## --- Plotting
    h1 = plot(ylabel="Age (Ma)", xlabel="", legend=:topleft, xaxis=false, xticks=false, yflip=true, xlims=(-3, 45), ylims=(65.3, 66.70), ygridalpha=0.3, margins=5mm, yticksfontsize=9)

    #Magnetic reversals (Ogg 2021, Geomagnetic polarity timescale)
        rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h]) #Define plotting a rectangle (coordinates of the four vertices. Negative height here because y axis decreases upwards)
        plot!(h1, rectangle(-3, -(68.178-66.380), 0, 68.178), opacity=1, color=:black, linecolor=:black, label=""); annotate!(h1, [-3], [66.380], text.(["C30N"], 7, :top, color=:gray)) #C30n
        plot!(h1, rectangle(-3, -(66.380-65.700), 0, 66.380), opacity=1, color=:white, linecolor=:black, label=""); annotate!(h1, [-3], [66.380], text.(["C29R"], 7, :bottom, color=:black)) #C29r
        plot!(h1, rectangle(-3, -(65.700-64.862), 0, 65.700), opacity=1, color=:black, linecolor=:black, label=""); annotate!(h1, [-3], [65.700], text.(["C29N"], 7, :bottom, color=:gray)) #C29n
        plot!(h1, rectangle(-3, -(64.862-64.645), 0, 64.862), opacity=1, color=:white, linecolor=:black, label=""); annotate!(h1, [-3], [64.645], text.(["C28R"], 7, :top, color=:black)) #C28r
        plot!(h1, rectangle(-3, -(64.645-63.537), 0, 64.645), opacity=1, color=:black, linecolor=:black, label=""); annotate!(h1, [-3], [64.645], text.(["C28N"], 7, :bottom, color=:gray)) #C28n
        plot!(h1, rectangle(-3, -(63.537-62.530), 0, 63.537), opacity=1, color=:white, linecolor=:black, label=""); annotate!(h1, [-3], [63.537], text.(["C27R"], 7, :bottom, color=:black)) #C27r
        plot!(h1, rectangle(-3, -(62.530-62.278), 0, 62.530), opacity=1, color=:black, linecolor=:black, label=""); annotate!(h1, [-3], [62.530], text.(["C27N"], 7, :bottom, color=:gray)) #C27n


    #Cretaceous-Paleogene boundary (Schoene et al. 2019 recalculation of Clyde et al. 2016 U-Pb date: 66.016±0.050 Ma (0.099Ma for full systematic uncertainty))
        plot!(h1, [-3,45], [66.016, 66.016], ribbon=0.050, seriestype=:line, color=:firebrick2, markerstrokecolor=:firebrick2, label="", linealpha=1, fillalpha=0.3, series_annotations=[text("KPB",9,:bottom, :firebrick2),""])
    
    #Eruption ages
        for i in 1:length(Malwa_samples)
            if i == 1
                plot!(h1, [i], [Malwa_eruptive_ages[i]], yerror=([Malwa_95CI_m[i]], [Malwa_95CI_p[i]]), color=ColorSchemes.tol_PRGn[2], markerstrokecolor=ColorSchemes.tol_PRGn[2], seriestype=:scatter, label="Malwa Plateau Bayesian MCMC Stratigraphic Model Ages")
            else 
                plot!(h1, [i], [Malwa_eruptive_ages[i]], yerror=([Malwa_95CI_m[i]], [Malwa_95CI_p[i]]), color=ColorSchemes.tol_PRGn[2], markerstrokecolor=ColorSchemes.tol_PRGn[2], seriestype=:scatter, label="")
            end
        end

        for j in 1:length(WesternGhats_samples)
            if j == 1
                plot!(h1, (j+length(Malwa_samples))*ones(1), [WesternGhats_eruptive_ages[j]], yerror=([WesternGhats_95CI_m[j]], [WesternGhats_95CI_p[j]]), color=:black, markerstrokecolor=:black, seriestype=:scatter, label="Western Ghats Bayesian MCMC Stratigraphic Model Ages")
            else 
                plot!(h1, (j+length(Malwa_samples))*ones(1), [WesternGhats_eruptive_ages[j]], yerror=([WesternGhats_95CI_m[j]], [WesternGhats_95CI_p[j]]), color=:black, markerstrokecolor=:black, seriestype=:scatter, label="")
            end
        end

        for k in 1:length(Mandla_samples) #In figure, add in a star to distinguish the MDAs (RBFQ and RBEC, first two in rank order) from the Bayesian eruption ages
            if k == 1
                plot!(h1, (k+length(Malwa_samples)+length(WesternGhats_samples))*ones(1), [Mandla_eruptive_ages[k]], yerror=([Mandla_95CI_m[k]], [Mandla_95CI_p[k]]), color=ColorSchemes.YlGnBu[6], markerstrokecolor=ColorSchemes.YlGnBu[6], seriestype=:scatter, label="Mandla Lobe Bayesian Eruption Ages and MDAs")
            else 
                plot!(h1, (k+length(Malwa_samples)+length(WesternGhats_samples))*ones(1), [Mandla_eruptive_ages[k]], yerror=([Mandla_95CI_m[k]], [Mandla_95CI_p[k]]), color=ColorSchemes.YlGnBu[6], markerstrokecolor=ColorSchemes.YlGnBu[6], seriestype=:scatter, label="")
            end
        end
        

    #Save and display figure
        display(h1)
        savefig(h1, "AllDeccanAges_RankOrder.pdf")
        savefig(h1, "AllDeccanAges_RankOrder.svg")

## --- End of file