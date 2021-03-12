push!(LOAD_PATH, ".")
# using MildPathology # Top level packages
using Plots, JLD
# using GR

# mu_list = 0:0.005:50
# n0_limit_list = constructBelt(2.5, mu_list)


# scanBKGGif(2,4)
# bkgScanGif(4,6)
# bkgScanGif(20,22)

# dict_test = selectMuRegion(mu_list, n0_limit_list)

# for i in 0:30
#     println(i, "\t",dict_test[i])
# end
# plot(n0_limit_list, mu_list, label = ["upper" "lower"], xlim = [0,15], ylim = [0,15]) # for plot
# plot(n0_limit_list[:,1], mu_list, n0_limit_list[:,2], mu_list, xlim = [0,15], ylim = [0,15]) # for GR

# plot uppper limit for n = 0
upper_mu_dict = load("/Users/fangchanghao/Documents/rareEventsAnalysis/analysisCode/src/likelihood_mu_bkg_dat/mu2_0_20_0.9.jld")
bkg_scan = upper_mu_dict["bkg_scan"]
# upper_mu_for_0 = ones(length(bkg_scan))
# upper_mu_for_1 = ones(length(bkg_scan))
# upper_mu_for_2 = ones(length(bkg_scan))
# upper_mu_for_3 = ones(length(bkg_scan))
# upper_mu_for_4 = ones(length(bkg_scan))
# for i in 1:length(bkg_scan)
#     # n0_limit_list = constructBelt(bkg_scan[i])
#     # upper_mu_dict = selectMuRegion(mu_list, n0_limit_list)
#     upper_mu_for_0[i] = upper_mu_dict["0"]
#     upper_mu_for_1[i] = upper_mu_dict["1"]
#     upper_mu_for_2[i] = upper_mu_dict["2"]
#     upper_mu_for_3[i] = upper_mu_dict["3"]
#     upper_mu_for_4[i] = upper_mu_dict["4"]
# end
fig = plot(bkg_scan, upper_mu_dict["0"])
plot!(bkg_scan, upper_mu_dict["1"])
plot!(bkg_scan, upper_mu_dict["2"])
plot!(bkg_scan, upper_mu_dict["3"])
plot!(bkg_scan, upper_mu_dict["4"])
# savefig(fig, "mu2_bkg.png")
# bkg_scan = 0:0.1:20
# mu_list = 0:0.005:50
# dic = Dict()
# a = getMu2ofBkg(bkg_scan, mu_list, 0.9)
# for i in 0:10
#     dic["$(i)"] = a[i]
# end
# save("mu_bkg_dict_.jld", dic)


# showSortedResult(0.9800, 5.0, 60)

# for i in 1:200
#     println(n0_limit_list[i,1], "\t",n0_limit_list[i,2], "\t",mu_list[i])
# end

# Caculate mu_2 as function of bkg.