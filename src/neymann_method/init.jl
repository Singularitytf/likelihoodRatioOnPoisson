push!(LOAD_PATH, ".")
using method, Plots

mu_list = 0:0.005:50
n0_limit_list = constructBelt(3.0, mu_list)
plot(n0_limit_list, mu_list, label = ["upper" "lower"], xlim = [0,15], ylim = [0,15])