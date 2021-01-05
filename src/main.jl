push!(LOAD_PATH, ".")
using lookItUp, Distributions, JLD, Plots, StatsBase, Random

# Loading the init data.
loadData(0.9)


list_ = rand(Poisson(1), 30000)
upper_limit_list = Float64[]
lower_limit_list = Float64[]
for i in list_
    temp = findInterval(i, 3.0)
    append!(upper_limit_list, temp[1])
    append!(lower_limit_list, temp[2])
end
# interval_list = findInterval.(list_, 3.0)

hist = append!(Histogram(0:0.1:2), upper_limit_list)
plot(hist)
# scatter(lower_limit_list, upper_limit_list)
# histogram2d(upper_limit_list, lower_limit_list)


# function fit_function(p::NamedTuple{(:a, :mu, :sigma)}, x::Real)
#     p.a[1] * pdf(Normal(p.mu[1], p.sigma), x) +
#     p.a[2] * pdf(Normal(p.mu[2], p.sigma), x)
# end

