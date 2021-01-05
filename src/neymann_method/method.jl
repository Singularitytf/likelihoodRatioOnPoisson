module method
export getMU2Dict, selectMuRegion, constructBelt, getMuofBkg

using poissonNeymann

function getMU2Dict(mu_list, n0_limit_list)
    #=
    This function gives the mu interval with dictionary structure.
    NOTICE: This function naturely select the topmost mu_2 of a fixed n_0, since it would overwrite the keys(n_0).
    =#
    length_of_muList::Int = length(mu_list)
    n0_upper_mu_dic = Dict{Int, Float64}()

    for i in 1:(length_of_muList-1)
        if n0_limit_list[i+1, 1] != n0_limit_list[i, 1]
            n0_upper_mu_dic[n0_limit_list[i, 1]] = mu_list[i]
        end
    end
    return n0_upper_mu_dic
end

function selectMuRegion(mu_list, n0_limit_list)
    #=
    This function gives the mu interval with dictionary structure.
    NOTICE: This function naturely select the topmost mu_2 of a fixed n_0, since it would overwrite the keys(n_0).
    =#
    length_of_muList::Int = length(mu_list)
    n0_upper_mu_dic = Dict{Int, Float64}()
    n0_lower_mu_dic = Dict{Int, Float64}()

    for i in 1:(length_of_muList-1)
        if n0_limit_list[i+1, 1] != n0_limit_list[i, 1]
            n0_upper_mu_dic[n0_limit_list[i, 1]] = mu_list[i]
        end
        if n0_limit_list[i+1, 2] != n0_limit_list[i, 2]
            n0_lower_mu_dic[n0_limit_list[i, 2]+1] = mu_list[i] # +1 for bottommost.
        end
    end
    return n0_upper_mu_dic, n0_lower_mu_dic
end


function constructBelt(bkg::Float64, mu_list, cfdent_level::Float64=0.9, n0_upper::Int=100)
    #=
    This function construct the confident belt with a certainty CL by scanning the parameter(mu) space.
    n0_upper -> scan limit for n0.
    mu_upper_limit -> scan limit for mu.
    RETURN: a 2 by len(mu_list) ARRAY.
    =#
    # Scan for the whole mu parameter space.
    # mu_list = 0:0.005:mu_upper_limit
    length_of_list = length(mu_list)
    n0_limit_list = ones(Int, length_of_list, 2)
    for i in 1:length_of_list
        n0_limit_list[i, :] = searchIntervalEdge(mu_list[i], bkg, cfdent_level, n0_upper)
    end
    return n0_limit_list
end

function getMuofBkg(bkg_scan, mu_list, cfdent_level; nrange = (0, 10))
    # Return a 2 by x array to get mu_2 against bkg plot.
    # after
    # mu_list = 0:0.005:50
    nlower, nupper = nrange
    upper_mu_for_ = Dict{Int, Array}()
    lower_mu_for_ = Dict{Int, Array}()

    # initialize data structure.
    for i in nlower:nupper
        upper_mu_for_[i] = ones(length(bkg_scan))
        lower_mu_for_[i] = ones(length(bkg_scan))
    end
    # end of initialization.
    for i in 1:length(bkg_scan)
        n0_limit_list = constructBelt(bkg_scan[i], mu_list, cfdent_level) # get constructBelt for a fixed bkg.
        upper_mu_dict, lower_mu_dict = selectMuRegion(mu_list, n0_limit_list) # get mu_2 for each of n0.
        # fill each upper mu corresponding to a n0 with a fixed b.
        for j in nlower:nupper
            try
                upper_mu_for_[j][i] = upper_mu_dict[j]
            catch e
                if isa(e, KeyError)
                    # println("WARNING: KeyMissing in bkg ", bkg_scan[i], " at n0 ", j, "!!!")
                    upper_mu_for_[j][i] = NaN
                    continue 
                end
            end

            try
                lower_mu_for_[j][i] = lower_mu_dict[j]
            catch e
                if isa(e ,KeyError)
                    lower_mu_for_[j][i] = 0.0
                end
            end
        end
    end
    return upper_mu_for_, lower_mu_for_
end

end # The end of module.