## --- Load packages and move to appropriate directory
    using Plots, StatGeochem, Distributions, Chron
    using ColorSchemes, Measures, DataFrames, CSV
    using QuadGK, Dierckx
    cd(@__DIR__)

## --- Color palette for plotting
    colors = ColorSchemes.twelvebitrainbow;

## --- Duration for probability density calculations
    duration = 50.0:0.001:80.0;

## --- Bilinear Exponential Distribution function from Chron.jl
    function bilinear_exponential(x, p)
        i₀ = firstindex(p)
        d = BilinearExponential(p[i₀], p[i₀+1], abs(p[i₀+2]), abs(p[i₀+3]), abs(p[i₀+4]))
        return pdf.(d, x)
    end

## --- All Zircons Western Ghats, n = 308
    #Data from Schoene et al. 2019. 206Pb/238U <Th> ages
    western_ghats_zircons = importdataset("WesternGhatsData/WesternGhatsZircons.csv", ',', importas=:Tuple)
    western_ghats_dates = western_ghats_zircons.Date #206/238U <Th> dates
    western_ghats_2σ_uncertainty = western_ghats_zircons.Uncertainty #2σ uncertainty

    western_ghats_distribution = zeros(length(duration))
    for i =1:length(western_ghats_dates)
        western_ghats_distribution += normpdf(western_ghats_dates[i], western_ghats_2σ_uncertainty[i] ./ 2, duration)
    end
    spl_western_ghats_distribution = Spline1D(duration, western_ghats_distribution); western_ghats_distribution_area, _ = quadgk(x -> spl_western_ghats_distribution.(x), duration[1], duration[end])
    western_ghats_distribution = western_ghats_distribution ./ western_ghats_distribution_area #Normalize by area under the curve

    n_westernghats_all = length(western_ghats_dates)

## --- Basalt Eruption Ages, Western Ghats, n = 24
    #Data from Schoene et al. 2019, Bayesian MCMC Stratigraphic Age Model results, 95% CI with internal uncertainties

    WesternGhats_samples_distparameters = DataFrame(CSV.File("/Users/tg5176/Documents/GitHub/MandlaLobe_UPb/WesternGhatsData/DeccanUPbData/distparameters_Schoene2019.csv"))

    #Calculate probability density
    WesternGhats_eruptionage_distribution = zeros(length(duration))
    for i =1:length(WesternGhats_samples_distparameters.Sample) 
        WesternGhats_eruptionage_distribution += bilinear_exponential(duration, WesternGhats_samples_distparameters[i,2:end]) #Bayesian eruption age distribution parameters (calculates pdf)
    end

    spl_WesternGhats_eruptionage = Spline1D(duration, WesternGhats_eruptionage_distribution); WesternGhats_eruptionage_area, _ = quadgk(x -> spl_WesternGhats_eruptionage.(x), duration[1], duration[end])
    WesternGhats_eruptionage_distribution = WesternGhats_eruptionage_distribution ./ WesternGhats_eruptionage_area

    n_westernghats_eruptionages = length(WesternGhats_samples_distparameters.Sample)

## --- All Zircons Malwa Plateau, n = 57
    #Data from Eddy et al 2020. 206Pb/238U <Th> ages. All the zircons from the samples interpretted as tuffaceous and included in age model
    malwa_zircons = importdataset("MalwaData/Malwazircons.csv", ',', importas=:Tuple)
    malwa_dates = malwa_zircons.Date #206/238U <Th> dates
    malwa_2σ_uncertainty = malwa_zircons.Uncertainty #2σ uncertainty

    malwa_distribution = zeros(length(duration))
    for i =1:length(malwa_dates)
        malwa_distribution += normpdf(malwa_dates[i], malwa_2σ_uncertainty[i] ./ 2, duration)
    end
    spl_malwa_distribution = Spline1D(duration, malwa_distribution); malwa_distribution_area, _ = quadgk(x -> spl_malwa_distribution.(x), duration[1], duration[end])
    malwa_distribution = malwa_distribution ./ malwa_distribution_area

    n_malwa_all = length(malwa_dates)

## --- Basalt Eruption Ages, Malwa Plateau, n = 7 
    #Data from Eddy et al. 2020, Bayesian MCMC stratigraphic age model outputs for bootstrapped prior, 95% CI with internal uncertainties

    Malwa_samples_distparameters = DataFrame(CSV.File("/Users/tg5176/Documents/GitHub/MandlaLobe_UPb/MalwaData/MalwaUPbData/distparameters_bootstrapped.csv"))

    #Calculate probability density
    Malwa_eruptionage_distribution = zeros(length(duration))
    for i =1:length(Malwa_samples_distparameters.Sample) 
        Malwa_eruptionage_distribution += bilinear_exponential(duration, Malwa_samples_distparameters[i,2:end]) #Bayesian eruption age distribution parameters (calculates pdf)
        #Malwa_eruptionage_distribution += normpdf(Malwa_eruptionages[i], Malwa_eruptionages_sigma[i], duration)
    end
    spl_Malwa_eruptionage = Spline1D(duration, Malwa_eruptionage_distribution); Malwa_eruptionage_area, _ = quadgk(x -> spl_Malwa_eruptionage.(x), duration[1], duration[end])
    Malwa_eruptionage_distribution = Malwa_eruptionage_distribution ./ Malwa_eruptionage_area

    n_Malwa_eruptionages = length(Malwa_samples_distparameters.Sample)

## --- Combined Western Ghats and Malwa Plateau distribution
    combined_malwawesternghats_eruptionage_distribution = zeros(length(duration))
       for i =1:length(Malwa_samples_distparameters.Sample) 
        combined_malwawesternghats_eruptionage_distribution += bilinear_exponential(duration, Malwa_samples_distparameters[i,2:end])
    end
    for i =1:length(WesternGhats_samples_distparameters.Sample) 
        combined_malwawesternghats_eruptionage_distribution += bilinear_exponential(duration, WesternGhats_samples_distparameters[i,2:end])
    end
    spl_combined_malwawesternghats_eruptionage = Spline1D(duration, combined_malwawesternghats_eruptionage_distribution); combined_malwawesternghats_eruptionage_area, _ = quadgk(x -> spl_combined_malwawesternghats_eruptionage.(x), duration[1], duration[end])
    combined_malwawesternghats_eruptionage_distribution = combined_malwawesternghats_eruptionage_distribution ./ combined_malwawesternghats_eruptionage_area

    n_malwawesternghats_all = length(Malwa_samples_distparameters.Sample) + length(WesternGhats_samples_distparameters.Sample)

## --- All Zircons Mandla Lobe, n = 85
    #Data from this work, 206Pb/238U <Th> ages included in the age interpretations
    Mandla_zircons = importdataset("MandlaData/Included_Mandla_Zircons.csv",',', importas=:Tuple)
    Mandla_dates = Mandla_zircons.Date_206Pb238U_Th #206/238 <Th> dates
    Mandla_2σ_uncertainty = Mandla_zircons.Date_206Pb238U_Th_2sigma #2σ uncertainty 

    Mandla_distribution = zeros(length(duration))
    for i =1:length(Mandla_dates)
        Mandla_distribution += normpdf(Mandla_dates[i], Mandla_2σ_uncertainty[i] ./ 2, duration)
    end
    spl_Mandla_distribution = Spline1D(duration, Mandla_distribution); Mandla_distribution_area, _ = quadgk(x -> spl_Mandla_distribution.(x), duration[1], duration[end])
    Mandla_distribution = Mandla_distribution ./ Mandla_distribution_area

    n_mandla_all = length(Mandla_dates)

## --- Basalt Eruption Ages, Mandla Lobe, n = 12
    #Data from this work, Bayesian age model results for all except RBEC and RBFQ (max depositional ages, youngest single zircon)
    RBEC = DataFrame(CSV.File("MandlaData/UPbData/RBEC.csv",header=false)); sort!(RBEC, order(:Column1,rev=true));
    RBFQ = DataFrame(CSV.File("MandlaData/UPbData/RBFQ.csv",header=false)); sort!(RBFQ, order(:Column1,rev=true));
    Mandla_samples_distparameters = DataFrame(CSV.File("MandlaData/UPbData/BayesianModeling/distparameters_bootstrapped.csv"));

    Mandla_sample_data = DataFrame(CSV.File("MandlaData/UPbData/BayesianModeling/distresults_bootstrapped.csv"))
    Mandla_samples = ["RBCT", "RBCU", "RBCV", "RBCZ", "RBCY", "RBCW", "RBFS", "RBFT", "RBEC", "RBFR", "RBEB", "RBFQ"]
    Mandla_eruptionages = []
    Mandla_eruptionages_95CI_p = []
    Mandla_eruptionages_95CI_m = []
    Mandla_eruptionages_sigma = []

    for i = 1:length(Mandla_samples)
        if Mandla_samples[i] == "RBEC" #Max depositional age, not bayesian. Single youngest zircon
            push!(Mandla_eruptionages, RBEC[end,1])
            push!(Mandla_eruptionages_95CI_p, RBEC[end,2])
            push!(Mandla_eruptionages_95CI_m, RBEC[end,2])
            push!(Mandla_eruptionages_sigma, RBEC[end,2])

        elseif Mandla_samples[i] == "RBFQ" #Max depositional age, not bayesian. Single youngest zircon
            push!(Mandla_eruptionages, RBFQ[end,1])
            push!(Mandla_eruptionages_95CI_p, RBFQ[end,2])
            push!(Mandla_eruptionages_95CI_m, RBFQ[end,2])
            push!(Mandla_eruptionages_sigma, RBFQ[end,2])

        else
            push!(Mandla_eruptionages, Mandla_sample_data.Age[i])
            push!(Mandla_eruptionages_95CI_p, Mandla_sample_data."97.5% CI"[i] .- Mandla_sample_data.Age[i])
            push!(Mandla_eruptionages_95CI_m, Mandla_sample_data.Age[i] .- Mandla_sample_data."2.5% CI"[i])
            push!(Mandla_eruptionages_sigma, Mandla_sample_data.sigma[i])

        end
    end

    #Mandla Bayesian ages distribution
    Mandla_eruptionage_distribution = zeros(length(duration))
    for i =1:length(Mandla_samples_distparameters.Sample)
        if Mandla_samples_distparameters.Sample[i] == "RBEC"
            Mandla_eruptionage_distribution += normpdf(RBEC[end,1], RBEC[end,2], duration) #Gaussian distribution for this MDA sample age
        elseif Mandla_samples_distparameters.Sample[i] == "RBFQ"
            Mandla_eruptionage_distribution += normpdf(RBFQ[end,1], RBFQ[end,2], duration) #Gaussian distribution for this MDA sample age
        else
            Mandla_eruptionage_distribution += bilinear_exponential(duration, Mandla_samples_distparameters[i,2:end]) #Bayesian eruption age distribution parameters (calculates pdf)
            #Mandla_eruptionage_distribution += normpdf(Mandla_eruptionages[i], Mandla_eruptionages_sigma[i], duration) #Assumes that the eruption ages are Guassian distributions
    
        end
    end
    spl_Mandla_eruptionage_distribution = Spline1D(duration, Mandla_eruptionage_distribution); Mandla_eruptionage_distribution_area, _ = quadgk(x -> spl_Mandla_eruptionage_distribution.(x), duration[1], duration[end])
    Mandla_eruptionage_distribution = Mandla_eruptionage_distribution ./ Mandla_eruptionage_distribution_area

    n_mandla_eruptionages = length(Mandla_samples_distparameters.Sample)

## --- Scaled probability density vs age
    #All individual zircons
        #Uncertainties are all the 2sigma uncertainties on individual zircon dates
        h1 = plot(xlabel="Age (Ma)", ylabel="Normalized probability", xlims=(65.0, 67.0), ylims=(0,3), framestyle=:box, seriestype=:line, legend=:topleft, margin=5mm)
        plot!(h1, duration, western_ghats_distribution, color=:black, label="Western Ghats \nAll Zircon Dates, n=$n_westernghats_all", linewidth=1.5)
        plot!(h1, duration, malwa_distribution, color=:purple, label="Malwa Plateau \nAll Zircon Dates, n=$n_malwa_all", linewsideth=1.5)
        plot!(h1, duration, Mandla_distribution, color=ColorSchemes.YlGnBu[6], label="Mandla Lobe \nAll Zircon Dates, n=$n_mandla_all", linewidth=1.5)

        display(h1)    
        savefig(h1, "Zircon_PDF.pdf")

    #Eruption ages
        #These use the eruption age distributions (95% CI for all the Bayesian ones, 2sigma for Mandla RBEC and RBFQ since maximum depositional age)
        h2 = plot(xlabel="Age (Ma)", ylabel="Normalized probability", xlims=(65.0, 67.0), ylims=(0,12), framestyle=:box, seriestype=:line, legend=:topright)
        plot!(h2, duration, WesternGhats_eruptionage_distribution, color=:black, label="Basalt Eruption Ages \nWestern Ghats, n=$n_westernghats_eruptionages", linewidth=1.5)
        plot!(h2, duration, Malwa_eruptionage_distribution, color=:purple, label="Basalt Eruption Ages \nMalwa Plateau, n=$n_Malwa_eruptionages", linewidth=1.5)
        plot!(h2, duration, Mandla_eruptionage_distribution, color=ColorSchemes.YlGnBu[6], label="Basalt Eruption Ages \n Mandla Lobe, n =$n_mandla_eruptionages", linewidth=1.5)
        display(h2)    
        savefig(h2, "Zircon_PDF_EruptionAges.pdf")

## --- Scaled probability vs age for Pulsed eruption model and Mandla Bayesian ages
    h3 = plot(xlabel="Age (Ma)", ylabel="Normalized probability", xlims=(65.2, 66.8), ylims=(0,11.5), framestyle=:box, seriestype=:line, legend=:right, size=(900,600), margin=5mm, xtickfontsize=9, ytickfontsize=9)

    #Plot scaled probability vs age for Western Ghats eruption ages only, like pulsed eruption model plot in Schoene et al. 2019
        #plot!(h3, duration, WesternGhats_eruptionage_distribution, color=:blue, label="Western Ghats \nEruption Ages")
        #plot!(h3, duration, WesternGhats_eruptionage_distribution, label="", fillrange=zeros(length(duration)), fillalpha=0.6, fillcolor=:blue, linealpha=0, ylims=ylims())
    
    #Plot scaled probability vs age for Malwa Plateau and Western Ghats eruption ages. Updates pulsed eruption model plot in Schoene et al. 2019 to include all high precision U-Pb dates
        plot!(h3, duration, combined_malwawesternghats_eruptionage_distribution, color=:black, label="Western Ghats & Malwa \nPlateau Eruption Ages")
        plot!(h3, duration, combined_malwawesternghats_eruptionage_distribution, label="", fillrange=zeros(length(duration)), fillalpha=0.6, fillcolor=:black, linealpha=0, ylims=ylims())

    #KPg boundary
        rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
        plot!(h3, rectangle((2*0.050), 10, (66.016-0.050), 0), opacity=0.3, label="", color=:firebrick2, linecolor=:firebrick2)
        vline!(h3,[66.016], color=:firebrick2, markerstrokecolor=:firebrick2, label="K-Pg boundary", linealpha=1)

    #Sort Mandla eruption ages into rank order so easier to read figure
        t_Mandla = sortperm(Mandla_eruptionages, rev=true)
        Mandla_eruptionages_sorted = Mandla_eruptionages[t_Mandla]
        Mandla_eruptionages_95CI_p_sorted = Mandla_eruptionages_95CI_p[t_Mandla]
        Mandla_eruptionages_95CI_m_sorted = Mandla_eruptionages_95CI_m[t_Mandla]
        Mandla_samples_sorted = Mandla_samples[t_Mandla]
        Mandla_colors_sorted = colors[t_Mandla]

    #Plot Mandla Bayesian eruption ages
        position = 0.8
       
        for i in 1:length(Mandla_eruptionages)
            sample_name = Mandla_samples_sorted[i]
            plot!(h3, [Mandla_eruptionages_sorted[i]], [position], xerror=([Mandla_eruptionages_95CI_m_sorted[i]],[Mandla_eruptionages_95CI_p_sorted[i]]), seriestype=:scatter, color=Mandla_colors_sorted[i], markerstrokecolor=Mandla_colors_sorted[i], label="", markersize=5)
            annotate!(h3, [Mandla_eruptionages_sorted[i]], [position], text.(["$sample_name"], 9, :top))
            position += 0.8
        end

    #(Optional) Add Western Ghats formatiions at the top of the figure
        #Starting and ending point of formation shapes on y axis
            y1 = 10
            y2 = 10.5

        #Each shape is a trapezoid with edges based on the Bayesian age model outputs with triangular prior from Schoene et al. 2019 for the samples within that unit.
        #Slope of trapezoid reflects the 95% confidence interval of the age model output for the sample
        
        #Lonavala and Kalsubai subgroups (individual formations)
            ##No dates given for Igatpuri, Neral, or Bhimashankar formations. Bracketed by the ones with dates
            #plot!(h3, Shape(reverse!([66.296+0.037, 66.296,  66.296, 66.296-0.030], dims=1), [y1, y2, y2, y1]), opacity=1, color=:mediumpurple2, linecolor=:black, label="") #Jawhar
            #annotate!(h3, [66.296+0.037], [y2], text.(["Jawhar"], 7, :top))
            ##Igatpuri: No dates
            ##Neral: No dates
            #plot!(h3, Shape(reverse!([66.225+0.077, 66.225,  66.185, 66.185-0.056], dims=1), [y1, y2, y2, y1]), opacity=1, color=:lavender, linecolor=:black, label="") #Thakurvadi
            #annotate!(h3, [66.225+0.077], [y2], text.(["Thakurvadi"], 7, :top))
            ##Bhimashankar: No dates
            #plot!(h3, Shape(reverse!([66.161+0.066, 66.161,  66.161, 66.161-0.069], dims=1), [y1, y2, y2, y1]), opacity=1, color=:palegreen, linecolor=:black, label="") #Khandala
            #annotate!(h3, [66.161+0.066], [y2], text.(["Khandala"], 7, :top))
            #plot!(h3, Shape(reverse!([66.132+0.069, 66.132,  66.132, 66.132-0.058], dims=1), [y1, y2, y2, y1]), opacity=1, color=:olivedrab1, linecolor=:black, label="") #Bushe
            #annotate!(h3, [66.132+0.069], [y2], text.(["Bushe"], 7, :top))
        
        #Lonavala and Kalsubai subgroups (combined)
            plot!(h3, Shape(reverse!([66.296+0.037, 66.296,  66.132, 66.132-0.058], dims=1), [y1, y2, y2, y1]), opacity=1, color=:gray, linecolor=:black, label="")
            annotate!(h3, [66.296+0.037], [y2], text.(["Lonavala and \n Kalsubai Subgroups"], 7, :top))
       
        #Wai subgroup
            plot!(h3, Shape(reverse!([66.088+0.032, 66.088,  66.039, 66.039-0.030], dims=1), [y1, y2, y2, y1]), opacity=1, color=:yellow, linecolor=:black, label="") #Polapur
            annotate!(h3, [66.088+0.032], [y2], text.(["Poladpur"], 7, :top))
            plot!(h3, Shape(reverse!([65.926+0.035, 65.926,  65.875, 65.875-0.042], dims=1), [y1, y2, y2, y1]), opacity=1, color=:orange, linecolor=:black, label="") #Ambenali
            annotate!(h3, [65.926+0.035], [y2], text.(["Ambenali"], 7, :top))
            plot!(h3, Shape(reverse!([65.631+0.053, 65.631,  65.590, 65.590-0.027], dims=1), [y1, y2, y2, y1]), opacity=1, color=:coral, linecolor=:black, label="") #Mahabaleswhar
            annotate!(h3, [65.631+0.053], [y2], text.(["Mahab."], 7, :top))

    #Save and display figure
        display(h3)
        savefig(h3, "Pulsedmodel_RankOrder.pdf")
        savefig(h3, "Pulsedmodel_RankOrder.svg")

## --- Alternative plotting of above figure, pulsed eruptions and Mandla Bayesian zircon ages as pdfs (much busier)

    #Plot scaled probability density vs age, like Basu et al. 2021 Figure 5
        h4 = plot(xlabel="Age (Ma)", ylabel="Normalized probability", xlims=(65.2, 66.8), ylims=(0,15), framestyle=:box, seriestype=:line, legend=:topright, size=(900,600), margin=5mm, xtickfontsize=9, ytickfontsize=9)
        plot!(h4, duration, combined_malwawesternghats_eruptionage_distribution, color=:black, label="Malwa and Western Ghats \nEruption Ages")
        plot!(h4, duration, combined_malwawesternghats_eruptionage_distribution, label="", fillrange=zeros(length(duration)), fillalpha=0.6, fillcolor=:black, linealpha=0, ylims=ylims())
        #plot!(h4, duration, Mandla_eruptionage_distribution, color=ColorSchemes.YlGnBu[6], label="Mandla Lobe Eruption Ages", linewidth=1.5)

   #Plot each Mandla red bole's eruption age as pdf distribution (bootstrapped Bayesian ages (bilinear exponential distribution) for all but the MDA ages of RBFQ and RBEC)
        for i = 1:length(Mandla_eruptionages)
            sample_name = Mandla_samples[i]
            if sample_name == "RBFQ" || sample_name =="RBEC"
                dist = normpdf(Mandla_eruptionages[i], Mandla_eruptionages_sigma[i], duration)
                spl_dist = Spline1D(duration, dist); dist_area, _ = quadgk(x -> spl_dist.(x), duration[1], duration[end])
                dist = normpdf(Mandla_eruptionages[i], Mandla_eruptionages_sigma[i], duration) ./ dist_area #Normalize by area under the curve
                plot!(h4, duration, dist, color=colors[i], label="$sample_name", linewidth=2)
            else
                dist = bilinear_exponential(duration, Mandla_samples_distparameters[i,2:end])
                spl_dist = Spline1D(duration, dist); dist_area, _ = quadgk(x -> spl_dist.(x), duration[1], duration[end])
                dist = bilinear_exponential(duration, Mandla_samples_distparameters[i,2:end]) ./ dist_area #Normalize by area under the curve
                plot!(h4, duration, dist, color=colors[i], label="$sample_name", linewidth=2)
            end
        end
    
    #KPg boundary (Schoene et al. 2019 recalculation of Clyde et al. 2016 U-Pb date, internal uncertainties only)
        rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
        plot!(h4, rectangle((2*0.050), 20, (66.016-0.050), 0), opacity=0.3, label="", color=:firebrick2, linecolor=:firebrick2)
        vline!(h4,[66.016], color=:firebrick2, markerstrokecolor=:firebrick2, label="K-Pg boundary", linealpha=1)
       
    display(h4)
    savefig(h4, "Pulsedmodel_PDF.pdf")

## --- End of file