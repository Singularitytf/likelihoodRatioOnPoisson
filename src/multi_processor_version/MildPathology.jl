push!(LOAD_PATH, ".")
include("HighPerformenceLHR.jl")
using Distributed


function getMu2ofBkg(bkg_scan, mu_list, cfdent_level; nrange = (0, 10))
    # Return a 2 by x array to get mu_2 against bkg plot.
    # after
    # mu_list = 0:0.005:50
    upper_mu_for_ = Dict{Int, Array}()
    nlower, nupper = (nrange)

    # initialize data structure.
    for i in nlower:nupper
        upper_mu_for_[i] = ones(length(bkg_scan))
    end
    # end of initialization.
    for i in 1:length(bkg_scan)
        n0_limit_list = constructBelt(bkg_scan[i], mu_list, cfdent_level) # get constructBelt for a fixed bkg.
        upper_mu_dict = selectMuRegion(mu_list, n0_limit_list) # get mu_2 for each of n0.
        # fill each upper mu corresponding to a n0 with a fixed b.
        for j in nlower:nupper
            upper_mu_for_[j][i] = upper_mu_dict[j]
        end
    end
    return upper_mu_for_
end

function getStepIndex(list, bkg_scan_precision)
    # input an dict mapping n0 to mu_2(b). output its step's index
    
    step_index_list = []

    len = length(list)
    for i in 1:(len-1)
        if (abs(list[i+1] - list[i])/bkg_scan_precision) > 5 # slightly larger than mu precision. If mu_2 is non-continous, it should larger than mu_precision.
            push!(step_index_list, i+1) # push peak index.
        end
    end
    return step_index_list
end