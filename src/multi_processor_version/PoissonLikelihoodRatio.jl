#=
The algorthim is form PHYSICAL REVIEW D Unified approach to the classical statistical analysis of small signals.
=#


module PoissonLikelihoodRatio
export likelihoodRatio, selectMuRegion, getSortedRatio, showSortedResult, constructBelt
using Distributions


function log_poisson(mu, n0)
    sum_ = sum(log.(1:n0))
    return n0 * log(mu) - mu - sum_
end

function poisson_pdf(mu, n0)
    return exp(log_poisson(mu, n0))
end


function getMuBest(n0::Int, bkg::Float64)
    #=
    This function will not be exported. No need to input a object
    =#
    if (n0 - bkg) >= 0
        mu_best = n0 - bkg
    else
        mu_best = 0.0
    end
    return mu_best
end

# function likelihoodRatio(mu::Int, bkg::Float64, n0::Int)
#     #=
#     This function will not be exported. No need to input a object
#     =#
#     mu_best = getMuBest(n0, bkg)
#     return ((mu+bkg) / (mu_best+bkg))^n0 * exp(mu_best - mu)
# end

function getSortedRatio(mu, bkg, n0_upper, check::Bool = false)
    # A sub-purcess of likelihoodRatio.
    #=
    for each n0:
    1. Get mu_best for each n0.
    2. Compute p(n0|mu) and p(n0|mu_best).
    3. Compute likehood ratio, and save it as a list.
    4. Arrange the list, and take first # n for building confident interval.
    =#
    # For convinience of checking data, we could return p(n0|mu) and p(n0|mu_best) list.
    if check == false
        ratio_list = ones(n0_upper+1) # Since n0 counts form 0 to n0_upper.
        for n0 in 0:n0_upper
            mu_best = getMuBest(n0, bkg)
            ratio_list[n0+1] = pdf(Poisson(mu + bkg), n0)/pdf(Poisson(mu_best + bkg), n0)
            # ratio_list[n0+1] = likelihoodRatio(mu, bkg, n0)
        end
        sorted_index_rank = sortperm(ratio_list, rev = true) # Due to index start form 1, while the n0 start form 0, there exists one difference between index and n0, i.e. n0 = index - 1.
        return sorted_index_rank
    else
        ratio_list = ones(n0_upper+1)
        hypo_mu_list = ones(n0_upper+1)
        hypo_mu_best_list = ones(n0_upper+1)
        for n0 in 0:n0_upper
            mu_best = getMuBest(n0, bkg)
            hypo_mu_list[n0+1] = poisson_pdf(mu + bkg, n0)
            hypo_mu_best_list[n0+1] = poisson_pdf(mu_best + bkg, n0)
            ratio_list[n0+1] = hypo_mu_list[n0+1]/hypo_mu_best_list[n0+1]
        end
        sorted_index_rank = sortperm(ratio_list, rev = true) # Due to index start form 1, while the n0 start form 0, there exists one difference between index and n0, i.e. n0 = index - 1.
        rank = sortperm(sorted_index_rank, rev = false)
        return [hypo_mu_list, hypo_mu_best_list, ratio_list, rank]
    end
end


function selectN0(mu, bkg, index_sorted_list, cfdent_level)
    # A sub-purcess of likelihoodRatio.
    # Select upper and lower n0 for each mu with A confident level.
    sum_of_hypo::Float64 = 0.0
    selected_n0 = []
    n0_upper = length(index_sorted_list) - 1
    i::Int8 = 0
    while sum_of_hypo < cfdent_level
        sum_of_hypo += pdf(Poisson(mu + bkg), index_sorted_list[i+1]-1)
        push!(selected_n0, index_sorted_list[i+1]-1)
        i += 1
    end
    return selected_n0
end

# Likelihood ratio method
function likelihoodRatio(mu::Float64, bkg::Float64, cfdent_level::Float64, n0_upper::Integer = 100, probability::Bool = false)
    index_sorted_list = getSortedRatio(mu, bkg, n0_upper)
    # println(index_sorted_list)
    # Adding all the p(n0|mu) at arrange with sorted index, until it reached confident level we need.
    # The following twe parameters are created for outputing the information we may want to see later.
    selected_n0 = selectN0(mu, bkg, index_sorted_list, cfdent_level)
    # Select upper and lower limit.
    n0_upper_limit = maximum(selected_n0)
    n0_lower_limit = minimum(selected_n0)
    # Return parameters we need.
    if probability == false
        return [n0_lower_limit, n0_upper_limit]
    else
        return [[n0_lower_limit, n0_upper_limit], sum_of_hypo]
    end
end

function selectMuRegion(mu_list, n0_limit_list)
    #=
    This function gives the mu interval with dictionary structure.
    NOTICE: This function naturely select the topmost mu_2 of a fixed n_0, since it would overwrite the keys(n_0).
    =#
    temp::Int = 0
    length_of_length::Int = length(mu_list)
    n0_upper_mu_dic = Dict{Int32, Float64}()
    for i in 1:(length_of_length-1)
        if n0_limit_list[i+1, 1] != n0_limit_list[i, 1]
            n0_upper_mu_dic[n0_limit_list[i, 1]] = mu_list[i]
        end
    end
    return n0_upper_mu_dic
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
    n0_limit_list = ones(length_of_list, 2)
    for i in 1:length_of_list
        n0_limit_list[i, :] = likelihoodRatio(mu_list[i], bkg, cfdent_level, n0_upper)
    end
    return n0_limit_list
end


    # length_of_list = 10000 # parameter used in [S0556-2821(98)00109-X] is 10,000.
# mu_upper_limit = 60
# bkg = 5.0
# cfdent_level = 0.9
# n0_upper = 60

# mu_list = range(0, mu_upper_limit, length = length_of_list)
# n0_limit_list = ones(length_of_list, 2)
# for i in 1:length_of_list
#     temp = likelihoodRatio(mu_list[i], bkg, cfdent_level, n0_upper)
#     n0_limit_list[i, :] = temp
# end
# for i in 100:200
#     println(n0_limit_list[i,1], "\t",n0_limit_list[i,2], "\t",mu_list[i])
# end
# plot(n0_limit_list[:, 1], mu_list, "r",
#      n0_limit_list[:, 2], mu_list, "r")
# dict_test = selectMuRegion(mu_list, n0_limit_list)
# for i in 0:100
#     println(i, "\t",dict_test[i])
# end
end # End of module.