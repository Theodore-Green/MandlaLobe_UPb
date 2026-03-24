## --- Load required packages
using Statistics, StatsBase, DelimitedFiles, SpecialFunctions
using Chron
using Plots; gr();
using ColorSchemes, Measures
using CSV, DataFrames

## --- Load data
    cd(@__DIR__)
    
    #TEA data
    #Uncertainties on TEA element concentration measurements are all percent uncertainties
    Redbole_TEA = DataFrame(CSV.File(joinpath(@__DIR__, "MandlaData/TEAData/Mandla_TEA_redbole_samples_all.csv"))); #TEA data for all redbole zircons except nbhs. Needs to be filtered for just those actually included in age interpretations
    GZ7_STD = DataFrame(CSV.File(joinpath(@__DIR__, "MandlaData/TEAData/Mandla_TEA_GZ7_standard.csv"))); #TEA data for the GZ7 standard
    ZrHf_STD = DataFrame(CSV.File(joinpath(@__DIR__, "MandlaData/TEAData/Mandla_TEA_ZrHf_standard.csv"))); #TEA data for the Zr/Hf standard
    TPBs = DataFrame(CSV.File(joinpath(@__DIR__, "MandlaData/TEAData/Mandla_TEA_TPBs.csv"))); #TEA data for all Mandla total procedural blanks run
    NBHs = DataFrame(CSV.File(joinpath(@__DIR__, "MandlaData/TEAData/Mandla_TEA_NBHs.csv"))); #TEA data for all redbole zircons that were nobody home

    #Redux data
    Redbole_dates = DataFrame(CSV.File("MandlaData/All_Mandla_zircons.csv"));

    #List of accepted redbole zircon dates
    Accepted_Mandla_zircons = DataFrame(CSV.File("MandlaData/Included_Mandla_Zircons.csv"));

    #Sample names
    Mandla_samples = ("RBCT", "RBCU", "RBCV", "RBCZ", "RBCY", "RBCW", "RBFS", "RBFT", "RBEC", "RBFR", "RBEB", "RBFQ");

    #Colors for plotting
    colors = ColorSchemes.twelvebitrainbow;

## --- Plotting Zr/Hf

    #Individual samples
    Zr_sigma = 0.01 .* Redbole_TEA.Zr_uncert .* Redbole_TEA.Zr
    Hf_sigma = 0.01 .* Redbole_TEA.Hf_uncert .* Redbole_TEA.Hf
    Zr_Hf_sigma = sqrt.(((Zr_sigma ./ Redbole_TEA.Zr).^2 .+ (Hf_sigma ./ Redbole_TEA.Hf).^2) .* (Redbole_TEA.Zr ./ Redbole_TEA.Hf).^2)

    for i = 1:length(Mandla_samples)
        sample_name = Mandla_samples[i]
        t1 = plot(xlabel="206Pb/238U Th-corrected Date (Ma)", ylabel="Zr/Hf", title="", framestyle=:box, seriestype=:scatter, legend=:bottomright, xlims=(64.9,67.6), ylims=(0,70), margin=5mm, xtickfontsize=12, ytickfontsize=12, xlabelfontsize=12, ylabelfontsize=12, legendfontsize=12)
    
        for j = 1:length(Redbole_TEA[:,1])
            for k = 1:length(Accepted_Mandla_zircons[:,1])

                if (Redbole_TEA.FullSampleID[j] == Accepted_Mandla_zircons.FullSampleID[k]) && (Redbole_TEA.SampleName[j] == Mandla_samples[i]) && (Redbole_TEA.Fraction[j] == Accepted_Mandla_zircons.Fraction[k])
                    sample = Redbole_TEA.FullSampleID[j]
                    zircon = Redbole_TEA.Fraction[j]
                    plot!(t1, [Accepted_Mandla_zircons.Date_206Pb238U_Th[k]], [Redbole_TEA.Zr[j]./Redbole_TEA.Hf[j]], xerror=(Accepted_Mandla_zircons.Date_206Pb238U_Th_2sigma[k]), yerror=(Zr_Hf_sigma[j]), color=colors[i], markerstrokecolor=:black, seriestype=:scatter, label="", markersize=7)
                end
            
            end
        end

        plot!(t1, [NaN], [NaN], label="$sample_name", seriestype=:scatter, color=colors[i], markerstrokecolor=:black)
        display(t1)
        savefig(t1, "MandlaData/TEAFigures/ZrHf_$sample_name.pdf")
    end

    #All redboles
    Zr_sigma = 0.01 .* Redbole_TEA.Zr_uncert .* Redbole_TEA.Zr
    Hf_sigma = 0.01 .* Redbole_TEA.Hf_uncert .* Redbole_TEA.Hf
    Zr_Hf_sigma = sqrt.(((Zr_sigma ./ Redbole_TEA.Zr).^2 .+ (Hf_sigma ./ Redbole_TEA.Hf).^2) .* (Redbole_TEA.Zr ./ Redbole_TEA.Hf).^2)

    t2 = plot(xlabel="206Pb/238U Th-corrected Date (Ma)", ylabel="Zr/Hf", title="", framestyle=:box, seriestype=:scatter, legend=:bottomright, xlims=(64.9,67.6), ylims=(0,70), margin=5mm, xtickfontsize=10, ytickfontsize=10)
    
    for i = 1:length(Mandla_samples)
        sample_name = Mandla_samples[i]
        for j = 1:length(Redbole_TEA[:,1])
            for k = 1:length(Accepted_Mandla_zircons[:,1])

                if (Redbole_TEA.FullSampleID[j] == Accepted_Mandla_zircons.FullSampleID[k]) && (Redbole_TEA.SampleName[j] == Mandla_samples[i]) && (Redbole_TEA.Fraction[j] == Accepted_Mandla_zircons.Fraction[k])
                    sample = Redbole_TEA.FullSampleID[j]
                    zircon = Redbole_TEA.Fraction[j]
                    plot!(t2, [Accepted_Mandla_zircons.Date_206Pb238U_Th[k]], [Redbole_TEA.Zr[j]./Redbole_TEA.Hf[j]], xerror=(Accepted_Mandla_zircons.Date_206Pb238U_Th_2sigma[k]), yerror=(Zr_Hf_sigma[j]), color=colors[i], markerstrokecolor=:black, seriestype=:scatter, label="", markersize=5)
                end
            
            end
        end

        plot!(t2, [NaN], [NaN], label="$sample_name", seriestype=:scatter, color=colors[i], markerstrokecolor=:black)
    end

    display(t2)
    savefig(t2, "MandlaData/TEAFigures/ZrHf_all.pdf")
    savefig(t2, "MandlaData/TEAFigures/ZrHf_all.svg")

## --- Plotting Dy/Yb
    
    #Error propagation
    Dy_sigma = 0.01 .* Redbole_TEA.Dy_uncert .* Redbole_TEA.Dy
    Yb_sigma = 0.01 .* Redbole_TEA.Yb_uncert .* Redbole_TEA.Yb
    Dy_Yb_sigma = sqrt.(((Dy_sigma ./ Redbole_TEA.Dy).^2 .+ (Yb_sigma ./ Redbole_TEA.Yb).^2) .* (Redbole_TEA.Dy ./ Redbole_TEA.Yb).^2)

    for i = 1:length(Mandla_samples)
        sample_name = Mandla_samples[i]
        t3 = plot(xlabel="206Pb/238U Th-corrected Date (Ma)", ylabel="Dy/Yb", title="", legend=:bottomright, framestyle=:box, seriestype=:scatter, xlims=(64.9,67.6), ylims=(0,0.6), margin=5mm, xtickfontsize=12, ytickfontsize=12, xlabelfontsize=12, ylabelfontsize=12, legendfontsize=12)
        
        for j = 1:length(Redbole_TEA[:,1])
            for k = 1:length(Accepted_Mandla_zircons[:,1])

                if (Redbole_TEA.FullSampleID[j] == Accepted_Mandla_zircons.FullSampleID[k]) && (Redbole_TEA.SampleName[j] == Mandla_samples[i]) && (Redbole_TEA.Fraction[j] == Accepted_Mandla_zircons.Fraction[k])
                    sample = Redbole_TEA.FullSampleID[j]
                    zircon = Redbole_TEA.Fraction[j]  
                    plot!(t3, [Accepted_Mandla_zircons.Date_206Pb238U_Th[k]], [Redbole_TEA.Dy[j]./Redbole_TEA.Yb[j]], xerror=(Accepted_Mandla_zircons.Date_206Pb238U_Th_2sigma[k]), yerror=(Dy_Yb_sigma[j]), color=colors[i], markerstrokecolor=:black, seriestype=:scatter, label="", markersize=7)
                end
            
            end
        end

        plot!(t3, [NaN], [NaN], label="$sample_name", seriestype=:scatter, color=colors[i])

        display(t3)
        savefig(t3, "MandlaData/TEAFigures/DyYb_$sample_name.pdf")
    end

## --- Plotting Th/U

    for i = 1:length(Mandla_samples)
        sample_name = Mandla_samples[i]
        t4 = plot(xlabel="206Pb/238U Th-corrected Date (Ma)", ylabel="Th/U", title="", legend=:bottomright, framestyle=:box, seriestype=:scatter, xlims=(64.9,67.6), ylims=(0,2), margin=5mm, xtickfontsize=12, ytickfontsize=12, xlabelfontsize=12, ylabelfontsize=12, legendfontsize=12)
        for j = 1:length(Redbole_dates[:,1])
            for k = 1:length(Accepted_Mandla_zircons[:,1])
                if (Redbole_dates.SampleName[j] == Mandla_samples[i]) && (Accepted_Mandla_zircons.FullSampleID[k] == Redbole_dates.FullSampleID[j]) && (Redbole_dates.Fraction[j] == Accepted_Mandla_zircons.Fraction[k])
                    name = Redbole_dates.FullSampleID[j]
                    zircon = Redbole_dates.Fraction[j]  
                    plot!(t4, [Redbole_dates.Date_206Pb238U_Th[j]], [Redbole_dates.Th_U[j]], xerror=(Redbole_dates.Date_206Pb238U_Th_2sigma[j]), color=colors[i], markerstrokecolor=:black, seriestype=:scatter, label="", markersize=7)
                end
            end
        end

        plot!(t4, [NaN], [NaN], label="$sample_name", seriestype=:scatter, color=colors[i])

        display(t4)
        savefig(t4, "MandlaData/TEAFigures/ThU_$sample_name.pdf")
    end

## --- Zr/Hf vs Th/U

    #Individual samples
    Zr_sigma = 0.01 .* Redbole_TEA.Zr_uncert .* Redbole_TEA.Zr
    Hf_sigma = 0.01 .* Redbole_TEA.Hf_uncert .* Redbole_TEA.Hf
    Zr_Hf_sigma = sqrt.(((Zr_sigma ./ Redbole_TEA.Zr).^2 .+ (Hf_sigma ./ Redbole_TEA.Hf).^2) .* (Redbole_TEA.Zr ./ Redbole_TEA.Hf).^2)
 
    for i = 1:length(Mandla_samples)
        sample_name = Mandla_samples[i]
        t5 = plot(xlabel="Th/U", ylabel="Zr/Hf", title="", framestyle=:box, seriestype=:scatter, legend=:bottomright, margin=5mm, ylims=(0,70), xlims=(0.3,1.8), xtickfontsize=12, ytickfontsize=12, xlabelfontsize=12, ylabelfontsize=12, legendfontsize=12)
        
        for j = 1:length(Redbole_TEA[:,1])
            for k = 1:length(Redbole_dates[:,1])
                for l = 1:length(Accepted_Mandla_zircons[:,1])
                    if (Redbole_TEA.FullSampleID[j] == Redbole_dates.FullSampleID[k] == Accepted_Mandla_zircons.FullSampleID[l]) && (Redbole_TEA.SampleName[j] == Mandla_samples[i]) && (Redbole_TEA.Fraction[j] == Redbole_dates.Fraction[k] == Accepted_Mandla_zircons.Fraction[l])
                        sample = Redbole_TEA.FullSampleID[j]
                        zircon = Redbole_TEA.Fraction[j]  
                        plot!(t5, [Redbole_dates.Th_U[k]], [Redbole_TEA.Zr[j]./Redbole_TEA.Hf[j]], yerror=(Zr_Hf_sigma[j]), color=colors[i], markerstrokecolor=:black, seriestype=:scatter, label="", markersize=7)
                    end
                end
            
            end
        end

        plot!(t5, [NaN], [NaN], label="$sample_name", seriestype=:scatter, color=colors[i], markerstrokecolor=:black)

        display(t5)
        savefig(t5, "MandlaData/TEAFigures/ZrHf_ThU_$sample_name.pdf")
    end

    #All redboles
    Zr_sigma = 0.01 .* Redbole_TEA.Zr_uncert .* Redbole_TEA.Zr
    Hf_sigma = 0.01 .* Redbole_TEA.Hf_uncert .* Redbole_TEA.Hf
    Zr_Hf_sigma = sqrt.(((Zr_sigma ./ Redbole_TEA.Zr).^2 .+ (Hf_sigma ./ Redbole_TEA.Hf).^2) .* (Redbole_TEA.Zr ./ Redbole_TEA.Hf).^2)
 
    t6 = plot(xlabel="Th/U", ylabel="Zr/Hf", title="", framestyle=:box, seriestype=:scatter, legend=:bottomright, margin=5mm, ylims=(0,70), xlims=(0.3,1.8), xtickfontsize=10, ytickfontsize=10)
    for i = 1:length(Mandla_samples)
        sample_name = Mandla_samples[i]
        
        for j = 1:length(Redbole_TEA[:,1])
            for k = 1:length(Redbole_dates[:,1])
                for l = 1:length(Accepted_Mandla_zircons[:,1])
                    if (Redbole_TEA.FullSampleID[j] == Redbole_dates.FullSampleID[k] == Accepted_Mandla_zircons.FullSampleID[l]) && (Redbole_TEA.SampleName[j] == Mandla_samples[i]) && (Redbole_TEA.Fraction[j] == Redbole_dates.Fraction[k] == Accepted_Mandla_zircons.Fraction[l])
                        sample = Redbole_TEA.FullSampleID[j]
                        zircon = Redbole_TEA.Fraction[j]  
                        plot!(t6, [Redbole_dates.Th_U[k]], [Redbole_TEA.Zr[j]./Redbole_TEA.Hf[j]], yerror=(Zr_Hf_sigma[j]), color=colors[i], markerstrokecolor=:black, seriestype=:scatter, label="", markersize=5)
                    end
                end
            
            end
        end

        plot!(t6, [NaN], [NaN], label="$sample_name", seriestype=:scatter, color=colors[i], markerstrokecolor=:black)
    end
    display(t6)
    savefig(t6, "MandlaData/TEAFigures/ZrHf_ThU_all.pdf")
    savefig(t6, "MandlaData/TEAFigures/ZrHf_ThU_all.svg")

## --- Standards
    #Zr/Hf
        Zr_sigma = 0.01 .* ZrHf_STD.Zr_uncert .* ZrHf_STD.Zr
        Hf_sigma = 0.01 .* ZrHf_STD.Hf_uncert .* ZrHf_STD.Hf
        Zr_Hf_sigma = sqrt.(((Zr_sigma ./ ZrHf_STD.Zr).^2 .+ (Hf_sigma ./ ZrHf_STD.Hf).^2) .* (ZrHf_STD.Zr ./ ZrHf_STD.Hf).^2)
            
        t7 = plot(xlabel="Standard Measurement Number", ylabel="Zr/Hf", title="", framestyle=:box, seriestype=:scatter, legend=:bottomright, xlims=(1,length(ZrHf_STD.Zr)), ylims=(30,60), margin=5mm)
        for j = 1:length(ZrHf_STD[:,1])
             plot!(t7, [j], [ZrHf_STD.Zr[j]./ZrHf_STD.Hf[j]], yerror=(Zr_Hf_sigma[j]), color=:black, markerstrokecolor=:black, seriestype=:scatter, label="", markersize=7)
        end

        #Mean and standard deviation
            mean_ZrHf = mean(ZrHf_STD.Zr ./ ZrHf_STD.Hf)
            m1 = round(mean_ZrHf; digits=3)
            std_ZrHf = std(ZrHf_STD.Zr ./ ZrHf_STD.Hf)
            m2 = round(std_ZrHf; digits=3)
            hline!(t7, [mean_ZrHf], ribbon=(ZrHf_std), fillalpha=0.3, color=:blue, label="Zr/Hf Mean = $m1, STD = $m2", linewidth=2)

        #Weighted mean
            weights_ZrHf = 1 ./ ((Zr_Hf_sigma).^2)
            weighted_mean_ZrHf = mean(ZrHf_STD.Zr ./ ZrHf_STD.Hf, weights(weights_ZrHf))
            weighted_mean_ZrHf_sigma = 1 ./ sqrt(sum(weights_ZrHf))
            MSWD_ZrHf = (1 /(length(weights_ZrHf) - 1)) * sum((((ZrHf_STD.Zr ./ ZrHf_STD.Hf) .- weighted_mean_ZrHf).^2) ./ (Zr_Hf_sigma).^2)
            weighted_mean_ZrHf_sigma = weighted_mean_ZrHf_sigma * sqrt(MSWD_ZrHf)
            weighted_mean_ZrHf_2sigma = weighted_mean_ZrHf_sigma * 2

            w1 = round(weighted_mean_ZrHf; digits=3)
            w2 = round(weighted_mean_ZrHf_2sigma; digits=3)
            #hline!(t7, [weighted_mean_ZrHf], ribbon=(weighted_mean_ZrHf_2sigma), color=:cyan, label="Zr/Hf Weighted Mean = $w1 ± $w2", linewidth=2)

        #Standard value
            hline!(t7, [50], color=:red, label="Zr/Hf Standard Value", linewidth=2) 

            display(t7)
            savefig(t7, "MandlaData/TEAFigures/ZrHf_Standard.pdf")

    #GZ7
        #Weighted mean and uncertainties for the GZ7 standard, various REEs. Filter out any negative values to get rid of any slightly wonky calibration corrections (Sm, Tm, and Lu have slightly wonky GZ7 values on analysis day 4). 
                #If want to use the whole dataset, cut the if statements in the lines below
            Sm = []; Sm_uncert = [];
            for i = 1:length(GZ7_STD.Sm)
                if GZ7_STD.Sm[i] > 0
                    push!(Sm, GZ7_STD.Sm[i])
                    push!(Sm_uncert, GZ7_STD.Sm_uncert[i])
                end
            end
            Sm_sigma = 0.01 .* Sm_uncert .* Sm; weights_Sm = 1 ./ (Sm_sigma .^2); weighted_mean_Sm = mean(Sm, weights(weights_Sm)); weighted_mean_Sm_sigma = 1 ./ sqrt(sum(weights_Sm)); MSWD_Sm = (1 /(length(weights_Sm) - 1)) * sum(((Sm .- weighted_mean_Sm).^2) ./ (Sm_sigma).^2); weighted_mean_Sm_sigma = weighted_mean_Sm_sigma * sqrt(MSWD_Sm);
            
            Nd = []; Nd_uncert = [];
            for i = 1:length(GZ7_STD.Nd)
                if GZ7_STD.Nd[i] > 0
                    push!(Nd, GZ7_STD.Nd[i])
                    push!(Nd_uncert, GZ7_STD.Nd_uncert[i])
                end
            end
            Nd_sigma = 0.01 .* Nd_uncert .* Nd; weights_Nd = 1 ./ (Nd_sigma .^2); weighted_mean_Nd = mean(Nd, weights(weights_Nd)); weighted_mean_Nd_sigma = 1 ./ sqrt(sum(weights_Nd)); MSWD_Nd = (1 /(length(weights_Nd) - 1)) * sum(((Nd .- weighted_mean_Nd).^2) ./ (Nd_sigma).^2); weighted_mean_Nd_sigma = weighted_mean_Nd_sigma * sqrt(MSWD_Nd);
            
            Eu = []; Eu_uncert = [];
            for i = 1:length(GZ7_STD.Eu)
                if GZ7_STD.Eu[i] > 0
                    push!(Eu, GZ7_STD.Eu[i])
                    push!(Eu_uncert, GZ7_STD.Eu_uncert[i])
                end
            end
            Eu_sigma = 0.01 .* Eu_uncert .* Eu; weights_Eu = 1 ./ (Eu_sigma .^2); weighted_mean_Eu = mean(Eu, weights(weights_Eu)); weighted_mean_Eu_sigma = 1 ./ sqrt(sum(weights_Eu)); MSWD_Eu = (1 /(length(weights_Eu) - 1)) * sum(((Eu .- weighted_mean_Eu).^2) ./ (Eu_sigma).^2); weighted_mean_Eu_sigma = weighted_mean_Eu_sigma * sqrt(MSWD_Eu);
            
            Gd = []; Gd_uncert = [];
            for i = 1:length(GZ7_STD.Gd)
                if GZ7_STD.Gd[i] > 0
                    push!(Gd, GZ7_STD.Gd[i])
                    push!(Gd_uncert, GZ7_STD.Gd_uncert[i])
                end
            end
            Gd_sigma = 0.01 .* Gd_uncert .* Gd; weights_Gd = 1 ./ (Gd_sigma .^2); weighted_mean_Gd = mean(Gd, weights(weights_Gd)); weighted_mean_Gd_sigma = 1 ./ sqrt(sum(weights_Gd)); MSWD_Gd = (1 /(length(weights_Gd) - 1)) * sum(((Gd .- weighted_mean_Gd).^2) ./ (Gd_sigma).^2); weighted_mean_Gd_sigma = weighted_mean_Gd_sigma * sqrt(MSWD_Gd);
            
            Tb = []; Tb_uncert = [];
            for i = 1:length(GZ7_STD.Tb)
                if GZ7_STD.Tb[i] > 0
                    push!(Tb, GZ7_STD.Tb[i])
                    push!(Tb_uncert, GZ7_STD.Tb_uncert[i])
                end
            end
            Tb_sigma = 0.01 .* Tb_uncert .* Tb; weights_Tb = 1 ./ (Tb_sigma .^2); weighted_mean_Tb = mean(Tb, weights(weights_Tb)); weighted_mean_Tb_sigma = 1 ./ sqrt(sum(weights_Tb)); MSWD_Tb = (1 /(length(weights_Tb) - 1)) * sum(((Tb .- weighted_mean_Tb).^2) ./ (Tb_sigma).^2); weighted_mean_Tb_sigma = weighted_mean_Tb_sigma * sqrt(MSWD_Tb);
            
            Dy = []; Dy_uncert = [];
            for i = 1:length(GZ7_STD.Dy)
                if GZ7_STD.Dy[i] > 0
                    push!(Dy, GZ7_STD.Dy[i])
                    push!(Dy_uncert, GZ7_STD.Dy_uncert[i])
                end
            end
            Dy_sigma = 0.01 .* Dy_uncert .* Dy; weights_Dy = 1 ./ (Dy_sigma .^2); weighted_mean_Dy = mean(Dy, weights(weights_Dy)); weighted_mean_Dy_sigma = 1 ./ sqrt(sum(weights_Dy)); MSWD_Dy = (1 /(length(weights_Dy) - 1)) * sum(((Dy .- weighted_mean_Dy).^2) ./ (Dy_sigma).^2); weighted_mean_Dy_sigma = weighted_mean_Dy_sigma * sqrt(MSWD_Dy);
            
            Ho = []; Ho_uncert = [];
            for i = 1:length(GZ7_STD.Ho)
                if GZ7_STD.Ho[i] > 0
                    push!(Ho, GZ7_STD.Ho[i])
                    push!(Ho_uncert, GZ7_STD.Ho_uncert[i])
                end
            end
            Ho_sigma = 0.01 .* Ho_uncert .* Ho; weights_Ho = 1 ./ (Ho_sigma .^2); weighted_mean_Ho = mean(Ho, weights(weights_Ho)); weighted_mean_Ho_sigma = 1 ./ sqrt(sum(weights_Ho)); MSWD_Ho = (1 /(length(weights_Ho) - 1)) * sum(((Ho .- weighted_mean_Ho).^2) ./ (Ho_sigma).^2); weighted_mean_Ho_sigma = weighted_mean_Ho_sigma * sqrt(MSWD_Ho);
            
            Er = []; Er_uncert = [];
            for i = 1:length(GZ7_STD.Er)
                if GZ7_STD.Er[i] > 0
                    push!(Er, GZ7_STD.Er[i])
                    push!(Er_uncert, GZ7_STD.Er_uncert[i])
                end
            end
            Er_sigma = 0.01 .* Er_uncert .* Er; weights_Er = 1 ./ (Er_sigma .^2); weighted_mean_Er = mean(Er, weights(weights_Er)); weighted_mean_Er_sigma = 1 ./ sqrt(sum(weights_Er)); MSWD_Er = (1 /(length(weights_Er) - 1)) * sum(((Er .- weighted_mean_Er).^2) ./ (Er_sigma).^2); weighted_mean_Er_sigma = weighted_mean_Er_sigma * sqrt(MSWD_Er);
            
            Tm = []; Tm_uncert = [];
            for i = 1:length(GZ7_STD.Tm)
                if GZ7_STD.Tm[i] > 0
                    push!(Tm, GZ7_STD.Tm[i])
                    push!(Tm_uncert, GZ7_STD.Tm_uncert[i])
                end
            end
            Tm_sigma = 0.01 .* Tm_uncert .* Tm; weights_Tm = 1 ./ (Tm_sigma .^2); weighted_mean_Tm = mean(Tm, weights(weights_Tm)); weighted_mean_Tm_sigma = 1 ./ sqrt(sum(weights_Tm)); MSWD_Tm = (1 /(length(weights_Tm) - 1)) * sum(((Tm .- weighted_mean_Tm).^2) ./ (Tm_sigma).^2); weighted_mean_Tm_sigma = weighted_mean_Tm_sigma * sqrt(MSWD_Tm);
            
            Yb = []; Yb_uncert = [];
            for i = 1:length(GZ7_STD.Yb)
                if GZ7_STD.Yb[i] > 0
                    push!(Yb, GZ7_STD.Yb[i])
                    push!(Yb_uncert, GZ7_STD.Yb_uncert[i])
                end
            end
            Yb_sigma = 0.01 .* Yb_uncert .* Yb; weights_Yb = 1 ./ (Yb_sigma .^2); weighted_mean_Yb = mean(Yb, weights(weights_Yb)); weighted_mean_Yb_sigma = 1 ./ sqrt(sum(weights_Yb)); MSWD_Yb = (1 /(length(weights_Yb) - 1)) * sum(((Yb .- weighted_mean_Yb).^2) ./ (Yb_sigma).^2); weighted_mean_Yb_sigma = weighted_mean_Yb_sigma * sqrt(MSWD_Yb);
            
            Lu = []; Lu_uncert = [];
            for i = 1:length(GZ7_STD.Lu)
                if GZ7_STD.Lu[i] > 0
                    push!(Lu, GZ7_STD.Lu[i])
                    push!(Lu_uncert, GZ7_STD.Lu_uncert[i])
                end
            end
            Lu_sigma = 0.01 .* Lu_uncert .* Lu; weights_Lu = 1 ./ (Lu_sigma .^2); weighted_mean_Lu = mean(Lu, weights(weights_Lu)); weighted_mean_Lu_sigma = 1 ./ sqrt(sum(weights_Lu)); MSWD_Lu = (1 /(length(weights_Lu) - 1)) * sum(((Lu .- weighted_mean_Lu).^2) ./ (Lu_sigma).^2); weighted_mean_Lu_sigma = weighted_mean_Lu_sigma * sqrt(MSWD_Lu);
        
        #Plotting 
            t8 = plot(xlabel="REE", ylabel="Element Concentration", framestyle=:box, seriestype=:scatter, xlims=(0,12), margins=5mm, xticks = ([1,2,3,4,5,6,7,8,9,10,11], ["Sm","Nd","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu"]))
            plot!(t8, [1], [weighted_mean_Sm], yerror = [2 * weighted_mean_Sm_sigma], seriestype=:scatter, color=:black, label="Weighted Mean GZ7 Analyses", markersize=5)
            plot!(t8, [2], [weighted_mean_Nd], yerror = [2 * weighted_mean_Nd_sigma], seriestype=:scatter, color=:black, label="", markersize=5)
            plot!(t8, [3], [weighted_mean_Eu], yerror = [2 * weighted_mean_Eu_sigma], seriestype=:scatter, color=:black, label="", markersize=5)
            plot!(t8, [4], [weighted_mean_Gd], yerror = [2 * weighted_mean_Gd_sigma], seriestype=:scatter, color=:black, label="", markersize=5)
            plot!(t8, [5], [weighted_mean_Tb], yerror = [2 * weighted_mean_Tb_sigma], seriestype=:scatter, color=:black, label="", markersize=5)
            plot!(t8, [6], [weighted_mean_Dy], yerror = [2 * weighted_mean_Dy_sigma], seriestype=:scatter, color=:black, label="", markersize=5)
            plot!(t8, [7], [weighted_mean_Ho], yerror = [2 * weighted_mean_Ho_sigma], seriestype=:scatter, color=:black, label="", markersize=5)
            plot!(t8, [8], [weighted_mean_Er], yerror = [2 * weighted_mean_Er_sigma], seriestype=:scatter, color=:black, label="", markersize=5)
            plot!(t8, [9], [weighted_mean_Tm], yerror = [2 * weighted_mean_Tm_sigma], seriestype=:scatter, color=:black, label="", markersize=5)
            plot!(t8, [10], [weighted_mean_Yb], yerror = [2 * weighted_mean_Yb_sigma], seriestype=:scatter, color=:black, label="", markersize=5)
            plot!(t8, [11], [weighted_mean_Lu], yerror = [2 * weighted_mean_Lu_sigma], seriestype=:scatter, color=:black, label="", markersize=5)

            #Published GZ7 value
            plot!(t8, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], [4.14, 2.65, 0.513, 15, 4.71, 53.3, 18.5, 82.6, 18.0, 179, 29.8], yerror= [0.40, 0.370, 0.067, 0.90, 0.19, 2.3, 0.7, 2.8, 0.4, 12, 1.4], seriestype=:scatter, color=:red, label="GZ7 Standard Value", markersize=5)
        
            display(t8)
            savefig(t8, "MandlaData/TEAFigures/GZ7_Standard.pdf")

## --- End of File