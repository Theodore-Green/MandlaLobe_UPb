## --- Load required packages
using Statistics, StatsBase, DelimitedFiles, SpecialFunctions
using Chron
using Plots; gr();
using ColorSchemes
using CSV, DataFrames
using Distributions, Measures
using QuadGK, Dierckx, DSP, LsqFit

## --- Required function for fitting Bayesian eruption age distribution from Chron.jl
    function bilinear_exponential(x, p)
        i₀ = firstindex(p)
        d = BilinearExponential(p[i₀], p[i₀+1], abs(p[i₀+2]), abs(p[i₀+3]), abs(p[i₀+4]))
        return pdf.(d, x)
    end

## --- Duration for probability density calculations
    duration = 50:0.001:80

## --- Import datasets
    cd(@__DIR__)
    Mandla_samples_distparameters = DataFrame(CSV.File("MandlaData/UPbData/BayesianModeling/distparameters_bootstrapped.csv"))
    RBEC = DataFrame(CSV.File("MandlaData/UPbData/RBEC.csv",header=false)); sort!(RBEC, order(:Column1,rev=true));
    RBFQ = DataFrame(CSV.File("MandlaData/UPbData/RBFQ.csv",header=false)); sort!(RBFQ, order(:Column1,rev=true));
    WesternGhats_samples_distparameters = DataFrame(CSV.File("WesternGhatsData/DeccanUPbData/distparameters_Schoene2019.csv"))
    Malwa_samples_distparamters = DataFrame(CSV.File("MalwaData/MalwaUPbData/distparameters_bootstrapped.csv"))

## --- Poladpur
    #Data from Schoene et al. 2019, Bayesian MCMC Stratigraphic Age Model results, 95% CI with internal uncertainties
    Poladpur_samples = ["RBBI", "RBBM", "RBBR", "RBBJ", "RBBS", "RBAW", "RBX"]
    Poladpur_samples_distparameters = filter(:Sample => x -> x in Poladpur_samples, WesternGhats_samples_distparameters)

    h1 = plot(xlabel="Age (Ma)", ylabel="Probability (unscaled)", xlims=(65.3, 66.2), framestyle=:box, seriestype=:line, legend=:topright, size=(900,600), margin=5mm)
    Poladpur_eruptionage_distribution = zeros(length(duration))
    for i = 1:length(Poladpur_samples)
        Poladpur_eruptionage_distribution += bilinear_exponential(duration, Poladpur_samples_distparameters[i,2:end])
    end

    plot!(h1, duration, Poladpur_eruptionage_distribution, color=:black, label="Poladpur Eruption Ages")
    
## --- Ambenali
    #Data from Schoene et al. 2019, Bayesian MCMC Stratigraphic Age Model results, 95% CI with internal uncertainties
    Ambenali_samples = ["RBBH", "RBBF", "RBB2", "RBAY", "RBAO", "RBAN", "RBP", "RBO"]
    Ambenali_samples_distparameters = filter(:Sample => x -> x in Ambenali_samples, WesternGhats_samples_distparameters)

    Ambenali_eruptionage_distribution = zeros(length(duration))
    for i =1:length(Ambenali_samples) 
        Ambenali_eruptionage_distribution += bilinear_exponential(duration, Ambenali_samples_distparameters[i,2:end])
    end

    plot!(h1, duration, Ambenali_eruptionage_distribution, color=:blue, label="Ambenali Eruption Ages")

## --- Mahabaleshwar
    #Data from Schoene et al. 2019, Bayesian MCMC Stratigraphic Age Model results, 95% CI with internal uncertainties
    Mahabaleshwar_samples = ["RBE", "RBF", "DEC13-08", "RBG"]
    Mahabaleshwar_samples_distparameters = filter(:Sample => x -> x in Mahabaleshwar_samples, WesternGhats_samples_distparameters)

    Mahabaleshwar_eruptionage_distribution = zeros(length(duration))
    for i =1:length(Mahabaleshwar_samples) 
        Mahabaleshwar_eruptionage_distribution += bilinear_exponential(duration, Mahabaleshwar_samples_distparameters[i,2:end])
    end

    plot!(h1, duration, Mahabaleshwar_eruptionage_distribution, color=:green, label="Mahabaleshwar Eruption Ages")

## --- Mandla Lobe eruption ages for each red bole. Bayesian age model results for all (not stratigraphic age-depth) except RBEC and RBFQ (maximum depositional ages, youngest single zircon)
    #Mandla Lobe data (this work)
    Mandla_samples = ["RBCT", "RBCU", "RBCV", "RBCZ", "RBCY", "RBCW", "RBFS", "RBFT", "RBEC", "RBFR", "RBEB", "RBFQ"]

## --- Probability that Mandla eruptions happen during the Western Ghats hiatuses (area under the curve where the distributions overlap)
    #Poladpur-Ambenali hiatus
        PoladpurAmbenali_hiatus = DataFrame(Sample = [], HiatusFraction = [], HiatusPercent = [], TotalArea = [], SumCalculatedAreas=[])

        for i = 1:length(Mandla_samples_distparameters.Sample)
            sample = Mandla_samples_distparameters.Sample[i]

            #Pdf distribution and total area under curve for chosen red bole sample. Total area should be ~1 because they're close to gaussian distributions
            if sample == "RBEC"
                redbole_distribution = normpdf(RBEC[end,1], RBEC[end,2], duration)
                totalarea_redbole = normcdf(RBEC[end,1], RBEC[end,2], duration[end]) - normcdf(RBEC[end,1], RBEC[end,2], duration[1])
            elseif sample == "RBFQ"
                redbole_distribution = normpdf(RBFQ[end,1], RBFQ[end,2], duration)
                totalarea_redbole = normcdf(RBFQ[end,1], RBFQ[end,2], duration[end]) - normcdf(RBFQ[end,1], RBFQ[end,2], duration[1])
            else
                redbole_distribution = bilinear_exponential(duration, Mandla_samples_distparameters[i,2:end])
                dist_parameters = BilinearExponential(Mandla_samples_distparameters[i,2], Mandla_samples_distparameters[i,3], Mandla_samples_distparameters[i,4], Mandla_samples_distparameters[i,5], Mandla_samples_distparameters[i,6])
                totalarea_redbole = cdf(dist_parameters, duration[end]) - cdf(dist_parameters, duration[1])
            end

            #Distributions of the overlap areas. Minimum values at each point in the overlap area, allows for shading in the overlap area on plots and calculating curves to define that area
            min_pdf_redbole_Poladpur = min.(Poladpur_eruptionage_distribution, redbole_distribution) #Distribution of the overlap area between Poladpur and chosen red bole
            min_pdf_redbole_Ambenali = min.(Ambenali_eruptionage_distribution, redbole_distribution) #Distribution of the overlap area between Ambenali and chosen red bole

            #Fit curves to the overlap areas
            spl_Poladpur = Spline1D(duration, min_pdf_redbole_Poladpur) #Spline fit. Need to define the curve where the red bole and Poladpur distributions overlap. Fits a curve to the points in the overlap area
            Poladpur_redbole_area, _ = quadgk(x -> spl_Poladpur.(x), duration[1], duration[end], order=100)

            spl_Ambenali = Spline1D(duration, min_pdf_redbole_Ambenali) #Spline fit. Need to define the curve where the red bole and Ambenali distributions overlap. Fits a curve to the points in the overlap area
            Ambenali_redbole_area, _ = quadgk(x -> spl_Ambenali.(x), duration[1], duration[end], order=100)

            #Find the distribution and calculate the area under the other parts of the curve (those not overlapping the Ambenali or Poladpur pulses) 
            notpulses_pdf = redbole_distribution .- min_pdf_redbole_Poladpur .- min_pdf_redbole_Ambenali
            
            spl_notpulses = Spline1D(duration, notpulses_pdf)
            #Start and stop points for area calculations are roughly the center of the pulse's distribution, just to easily split up the sections and prevent any double counting. Shouldn't be anything in the overlap zones anyway, so just makes the calculation easy
            afterAmbenali_area, _ = quadgk(x -> spl_notpulses.(x), duration[1], 65.9, order=100)
            beforePoladpur_area, _ = quadgk(x -> spl_notpulses.(x), 66.075, duration[end], order=100)
            areabetweenpulses_redbole, _ = quadgk(x -> spl_notpulses.(x), 65.9, 66.075, order=100)

            sum_area = afterAmbenali_area + Ambenali_redbole_area + Poladpur_redbole_area + beforePoladpur_area + areabetweenpulses_redbole #Should be same as total area, good to double check

            #Record the amount of the chosen red bole distribution that falls during the hiatus (fraction of total area)
            hiatus_fraction = areabetweenpulses_redbole ./ totalarea_redbole
            push!(PoladpurAmbenali_hiatus, [sample, hiatus_fraction, hiatus_fraction*100, totalarea_redbole, sum_area])

            #Plot to visualize the intersections (optional)
            h2 = plot(duration, Poladpur_eruptionage_distribution, label="Poladpur", lw=2, fillalpha=0.3, fillrange=0, xlabel="Age (Ma)", ylabel="Probability density (unscaled)", color=:yellow, xlims=(65.75,66.15))
            plot!(h2, duration, redbole_distribution, label="$sample", lw=2, fillalpha=0.3, fillrange=0, color=:blue)
            plot!(h2, duration, min_pdf_redbole_Poladpur, fillrange=0, seriesalpha=0.8, color=:green, label="Overlap Area")
            plot!(h2, duration, Ambenali_eruptionage_distribution, label="Ambenali", lw=2, fillalpha=0.3, fillrange=0, color=:red)
            plot!(h2, duration, min_pdf_redbole_Ambenali, fillrange=0, seriesalpha=0.8, color=:purple, label="Overlap Area")
            display(h2)
        end

        CSV.write("Hiatus_calculation_PoladpurAmbenali.csv", PoladpurAmbenali_hiatus)
    
    #Ambenali-Mahabaleshwar hiatus
        AmbenaliMahabaleshwar_hiatus = DataFrame(Sample = [], HiatusFraction = [], HiatusPercent = [], TotalArea = [], SumCalculatedAreas = [])

        for i = 1:length(Mandla_samples_distparameters.Sample)
            sample = Mandla_samples_distparameters.Sample[i]

            #Pdf distribution and total area under curve for chosen red bole sample. Total area should be ~1 because they're close to gaussian distributions
            if sample == "RBEC"
                redbole_distribution = normpdf(RBEC[end,1], RBEC[end,2], duration)
                totalarea_redbole = normcdf(RBEC[end,1], RBEC[end,2], duration[end]) - normcdf(RBEC[end,1], RBEC[end,2], duration[1])
            elseif sample == "RBFQ"
                redbole_distribution = normpdf(RBFQ[end,1], RBFQ[end,2], duration)
                totalarea_redbole = normcdf(RBFQ[end,1], RBFQ[end,2], duration[end]) - normcdf(RBFQ[end,1], RBFQ[end,2], duration[1])
            else
                redbole_distribution = bilinear_exponential(duration, Mandla_samples_distparameters[i,2:end])
                dist_parameters = BilinearExponential(Mandla_samples_distparameters[i,2], Mandla_samples_distparameters[i,3], Mandla_samples_distparameters[i,4], Mandla_samples_distparameters[i,5], Mandla_samples_distparameters[i,6])
                totalarea_redbole = cdf(dist_parameters, duration[end]) - cdf(dist_parameters, duration[1])
            end

            #Distributions of the overlap areas. Minimum values at each point in the overlap area, allows for shading in the overlap area on plots and calculating curves to define that area
            min_pdf_redbole_Ambenali = min.(Ambenali_eruptionage_distribution, redbole_distribution) #Distribution of the overlap area between Ambenali and chosen red bole
            min_pdf_redbole_Mahabaleshwar = min.(Mahabaleshwar_eruptionage_distribution, redbole_distribution) #Distribution of the overlap area between Mahabaleshwar and chosen red bole

            #Fit curves to the overlap areas
            spl_Ambenali = Spline1D(duration, min_pdf_redbole_Ambenali) #Spline fit. Need to define the curve where the red bole and Ambenali distributions overlap. Fits a curve to the points in the overlap area
            Ambenali_redbole_area, _ = quadgk(x -> spl_Ambenali.(x), duration[1], duration[end], order=100)

            spl_Mahabaleshwar = Spline1D(duration, min_pdf_redbole_Mahabaleshwar) #Spline fit. Need to define the curve where the red bole and Mahabaleshwar distributions overlap. Fits a curve to the points in the overlap area
            Mahabaleshwar_redbole_area, _ = quadgk(x -> spl_Mahabaleshwar.(x), duration[1], duration[end], order=100)

            #Find the distribution and calculate the area under the other parts of the curve (those not overlapping the Mahabaleshwar or Ambenali pulses) 
            notpulses_pdf = redbole_distribution .- min_pdf_redbole_Ambenali .- min_pdf_redbole_Mahabaleshwar
            
            spl_notpulses = Spline1D(duration, notpulses_pdf)
            #Start and stop points for area calculations are roughly the center of the pulse's distribution, just to easily split up the sections and prevent any double counting. Shouldn't be anything in the overlap zones anyway, so just makes the calculation easy
            afterMahabaleshwar_area, _ = quadgk(x -> spl_notpulses.(x), duration[1], 65.65, order=100)
            beforeAmbenali_area, _ = quadgk(x -> spl_notpulses.(x), 65.9, duration[end], order=100)
            areabetweenpulses_redbole, _ = quadgk(x -> spl_notpulses.(x), 65.65, 65.9, order=100)

            sum_area = afterMahabaleshwar_area + Mahabaleshwar_redbole_area + Ambenali_redbole_area + beforeAmbenali_area + areabetweenpulses_redbole #Should be same as total area, good to double check

            #Record the amount of the chosen red bole distribution that falls during the hiatus (fraction of total area)
            hiatus_fraction = areabetweenpulses_redbole ./ totalarea_redbole
            push!(AmbenaliMahabaleshwar_hiatus, [sample, hiatus_fraction, hiatus_fraction*100, totalarea_redbole, sum_area])

            #Plot to visualize the intersections (optional)
            h3 = plot(duration, Ambenali_eruptionage_distribution, label="Ambenali", lw=2, fillalpha=0.3, fillrange=0, xlabel="Age (Ma)", ylabel="Probability density (unscaled)", color=:yellow, xlims=(65.5,66.05))
            plot!(h3, duration, redbole_distribution, label="$sample", lw=2, fillalpha=0.3, fillrange=0, color=:blue)
            plot!(h3, duration, min_pdf_redbole_Ambenali, fillrange=0, seriesalpha=0.8, color=:green, label="Overlap Area")
            plot!(h3, duration, Mahabaleshwar_eruptionage_distribution, label="Mahabaleshwar", lw=2, fillalpha=0.3, fillrange=0, color=:red)
            plot!(h3, duration, min_pdf_redbole_Mahabaleshwar, fillrange=0, seriesalpha=0.8, color=:purple, label="Overlap Area")
            display(h3)
        end

        CSV.write("Hiatus_calculations_AmbenaliMahabaleshwar.csv", AmbenaliMahabaleshwar_hiatus)

## --- End of file