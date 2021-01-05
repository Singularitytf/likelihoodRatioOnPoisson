using JLD, Plots
a = load("/home/fang/Documents/likelihoodRatioOnPoisson/src/neymann_mu_bkg_data/mu1_0_20_0.9.jld")
fig = plot(a["bkg_scan"], a["10"], title = "bkg vs. n0=10", xlabel = "bkg", ylabel = "mu_2")
savefig(fig, "n0_10_bkg.png")