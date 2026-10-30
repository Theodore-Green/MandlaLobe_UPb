## --- Load required packages
using Plots; gr();
using Isoplot
using CSV, DataFrames, Measures

## --- Load dataset
cd(@__DIR__)
ET100_ratios = DataFrame(CSV.File(joinpath(@__DIR__, "ET100_ratios.csv"))) #Data is stored in csv in the format 207Pb/235U ratio	207Pb/235U 1sigma (abs)	    206Pb/238U ratio    206Pb/238U 1sigma (abs) 	Correlation coefficient	
ET100_dates = DataFrame(CSV.File(joinpath(@__DIR__, "ET100.csv"))); sort!(ET100_dates, :Date_206Pb238U)

## --- Turn ratios into UPbAnalysis objects (Isoplot.jl function) for Wetherhill concordia plot
    analyses = UPbAnalysis.(eachcol(ET100_ratios)...,)

    h1 = plot(xlabel="²⁰⁷Pb/²³⁵U", ylabel="²⁰⁶Pb/²³⁸U", framestyle=:box, title="ET100", margins=5mm)
    plot!(h1, analyses, color=:darkblue, fillalpha=0.4, linecolor=:black, label="")
    concordiacurve!(h1) # Add concordia curve
    savefig(h1, "ET100_concordia.pdf")
    display(h1)

## --- Rank-order plot of 206/238 dates
    h2 = plot(framestyle=:box, legend=:topleft, title="ET100", ylabel="²⁰⁶Pb/²³⁸U Age [Ma]", xlabel="", xticks=false)
    plot!(h2, 1:length(ET100_dates.Date_206Pb238U), ET100_dates.Date_206Pb238U, yerror=(ET100_dates.Uncertainty_2sigma), seriestype=:scatter, color=:darkblue, mscolor=:darkblue, label="", margins=5mm)
    hline!(h2, [100.206], color=:black, ribbon=(0.005), fillalpha=0.6, label="100.206 ± 0.005/0.025 Ma \n n = 39; MSWD = 2.76")
    hline!(h2, [100.206], color=:black, ribbon=(0.025), fillalpha=0.3, label="")
    savefig(h2, "ET100_rankorder.pdf")
    display(h2)

## --- End of file
