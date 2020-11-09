module PoissonLikelihoodRatio
export likelihoodRatio, selectMuRegion, showSortedResult

using Distributions
using GR


function poissonbkg(mu::Float64, bkg::Float64)
    #=
    Return poisson distribution(a distribution object) with background.
    lambda is poisson parameter.
    bkg is background. And lambda + bkg as the parameter of in-build poisson distribution by Julia Distributions.
    =#
    poisson_para = mu + bkg
    return Poisson(poisson_para)
end

function searchIntervalEdge(func, cfdentLevel::Float64)
    """
    Our search way is find interval confident level.
    """
    # Initialize Parameters
    alpha = 1 - cfdentLevel
    upper_limit::UInt16 = 0
    lower_limit::UInt16 = 0
    upper_n0 = 25 # Change this parameter!!!
    is_upper = true # If we find upper limit, these two parameters would help iteration not update the upper(or lower)_limit value and help us check whether it update or not.
    is_lower = true
    # Search for upper and lower limit.
    for n0 in 0:upper_n0
        cdf_temp = cdf(func, n0)
        if (cdf_temp >= alpha/2) & is_upper
            upper_limit = n0
            is_upper = false
        end
        if (cdf_temp >= (1 - alpha/2)) & is_lower
            lower_limit = n0
            is_lower = false
        end
    end
    # Note here!!!
    return [upper_limit, lower_limit]
end

function searchUpperEdge(func, cfdent_level::Float64)
    alpha = 1 - cfdent_level
    upper_limit::UInt16 = 0
    upper_n0 = 20 # Change this parameter!!!
    for n0 in 0:upper_n0
        cdf_temp = cdf(func, n0)
        if cdf_temp >= alpha
            upper_limit = n0
            break
        end
    end
    return upper_limit
end

function getMuBest(n0::Int, bkg::Float64)
    if (n0 - bkg) >= 0
        mu_best = n0 - bkg
    else
        mu_best = 0.0
    end
    return mu_best
end

function getSortedRatio(mu, bkg, n0_upper, check::Bool = false)
    #=
    for each n0:
    1. Scan for mu_best.
    2. Compute p(n0|mu) and p(n0|mu_best).
    3. Compute likehood ratio, and save it as a list.
    4. Arrange the list, and take first # n for building confident interval.
    =#
    # For convinience of checking data, we could return p(n0|mu) and p(n0|mu_best) list.
    if check == false
        ratio_list = ones(n0_upper+1)
        for n0 in 0:n0_upper
            mu_best = getMuBest(n0, bkg)
            ratio_list[n0+1] = pdf(Poisson(mu + bkg), n0)/pdf(Poisson(mu_best + bkg), n0)
        end
        sorted_index_rank = sortperm(ratio_list, rev = true) # Due to index start form 1, while the n0 start form 0, there exists one difference between index and n0, i.e. n0 = index - 1.
        return sorted_index_rank
    else
        ratio_list = ones(n0_upper+1)
        hypo_mu_list = ones(n0_upper+1)
        hypo_mu_best_list = ones(n0_upper+1)
        for n0 in 0:n0_upper
            mu_best = getMuBest(n0, bkg)
            hypo_mu_list[n0+1] = pdf(Poisson(mu + bkg), n0)
            hypo_mu_best_list[n0+1] = pdf(Poisson(mu_best + bkg), n0)
            ratio_list[n0+1] = hypo_mu_list[n0+1]/hypo_mu_best_list[n0+1]
        end
        sorted_index_rank = sortperm(ratio_list, rev = true) # Due to index start form 1, while the n0 start form 0, there exists one difference between index and n0, i.e. n0 = index - 1.
        rank = sortperm(sorted_index_rank, rev = false)
        return [hypo_mu_list, hypo_mu_best_list, ratio_list, rank]
    end
end

function showSortedResult(mu, bkg, n0_upper)
    # Input teturn from getSortedRatio function.
    hypo_mu_list, hypo_mu_best_list, ratio_list, rank = getSortedRatio(mu, bkg, n0_upper, true)
    println("n0", "\t", "hypo_mu", "\t", "hypo_mu_best", "\t", "ratio", "\t", "rank")
    for i in 1:length(rank)
        println(i-1, "\t",hypo_mu_list[i], "\t",hypo_mu_best_list[i], "\t",ratio_list[i], "\t",rank[i])
    end
end

function selectN0(mu, bkg, index_sorted_list, cfdent_level)
    # Select upper and lower n0 for each mu with A confident level.
    sum_of_hypo::Float64 = 0.0
    selected_n0 = []
    for i in 0:n0_upper
       if sum_of_hypo < cfdent_level
        sum_of_hypo = pdf(Poisson(mu + bkg), index_sorted_list[i+1] - 1) + sum_of_hypo
        push!(selected_n0, index_sorted_list[i+1] - 1)
        # print(sum_of_hypo, "\n") # This is for testing the part is work well or not.
       else
        break
       end
    end
    # Select upper and lower limit.
    n0_upper_limit = maximum(selected_n0)
    n0_lower_limit = minimum(selected_n0)
end

# Likelihood ratio method
function likelihoodRatio(mu, bkg, cfdent_level::Float64 = 0.99, n0_upper::Integer = 60, probability = false)
    index_sorted_list = getSortedRatio(mu, bkg, n0_upper)
    # println(index_sorted_list)
    # Adding all the p(n0|mu) at arrange with sorted index, until it reached confident level we need.
    # The following twe parameters are created for outputing the information we may want to see later.
    sum_of_hypo = 0.0
    selected_n0 = []
    for i in 0:n0_upper
       if sum_of_hypo < cfdent_level
        sum_of_hypo = pdf(Poisson(mu + bkg), index_sorted_list[i+1] - 1) + sum_of_hypo
        push!(selected_n0, index_sorted_list[i+1] - 1)
        # print(sum_of_hypo, "\n") # This is for testing the part is work well or not.
       else
        break
       end
    end
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