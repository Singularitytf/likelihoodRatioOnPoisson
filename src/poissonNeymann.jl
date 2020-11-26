module poissonNeymann
export poissonbkg, searchIntervalEdge, searchUpperEdge

using Distributions



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
end # The end of moudle.