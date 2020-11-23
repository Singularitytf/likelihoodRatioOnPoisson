module Tools
export scanBKGGif, getMuBkgPlot

push!(LOAD_PATH, ".")
using PoissonLikelihoodRatio
using Plots

function scanBKGGif(start_point, end_point)
    #=
    Generate a GIF animation for scanned bkg.
    start_point: of scanned bkg.
    end_point: of scanned bkg.
    =#
    plot_array = Any[]
    bkg_list = start_point:0.001:end_point
    mu_list = 0:0.005:20
    anim = @animate for i in bkg_list
        temp = constructBelt(i, 0.9, 60, 20.0) # import from PoissonLikelihoodRatio
        plot(temp, mu_list, 
            legend = false, 
            xlim = [0,15], 
            ylim = [0,15], 
            title="bkg = $(i)")
    end
gif(anim, "variousBKG$(start_point)_$(end_point).gif", fps=60)
end

function showSortedResult(mu, bkg, n0_upper::Int = 30)
    # Input teturn from getSortedRatio function.
    hypo_mu_list, hypo_mu_best_list, ratio_list, rank = getSortedRatio(mu, bkg, n0_upper, true)
    println("n0", "\t", "hypo_mu", "\t", "hypo_mu_best", "\t", "ratio", "\t", "rank")
    println("===="*10)
    for i in 1:length(rank)
        println(i-1, "\t", hypo_mu_list[i], "\t", 
        hypo_mu_best_list[i], "\t", 
        ratio_list[i], "\t",
        rank[i])
    end
end

function getMuBkgPlot(start_point::Float64, end_point::Float64, cfdent_level)
    # Return a 2 by x array to get mu_2 against bkg plot.
    # after
    mu_list = 0:0.005:50
    bkg_scan = start_point:0.1:end_point
    upper_mu_for_ = Dict{Int, Array}()

    # initialize data structure.
    for i in 0:10
        upper_mu_for_[i] = ones(length(bkg_scan))
    end
    # end of initialization.
    for i in 1:length(bkg_scan)
        n0_limit_list = constructBelt(bkg_scan[i], cfdent_level) # get constructBelt for a fixed bkg.
        upper_mu_dict = selectMuRegion(mu_list, n0_limit_list) # get mu_2 for each of n0.
        # fill each upper mu corresponding to a n0 with a fixed b.
        for j in 0:10
            upper_mu_for_[j][i] = upper_mu_dict[j]
        end
    end
    return bkg_scan, upper_mu_for_
end


end