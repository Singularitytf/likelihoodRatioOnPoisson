push!(LOAD_PATH, ".")
using MildPathology, PoissonLikelihoodRatio# Top level packages
using Plots, JLD
# using GR

mu_list = 0:0.005:50
n0_limit_list = constructBelt(3.0, mu_list)


# # # scanBKGGif(2,4)
# # # bkgScanGif(4,6)
# # # bkgScanGif(20,22)

dict_test1 = getMU2Dict(mu_list, n0_limit_list)

# for i in 0:20
#     println(i, "\t", dict_test1[i], "\t")
# end
# # plot(n0_limit_list, mu_list, label = ["upper" "lower"], xlim = [0,15], ylim = [0,15]) # for plot
# plot(n0_limit_list[:,1], mu_list, n0_limit_list[:,2], mu_list, xlim = [0,15], ylim = [0,15]) # for GR

# plot uppper limit for n = 0
# bkg_scan = 0:0.01:25
# upper_mu_for_0 = ones(length(bkg_scan))
# upper_mu_for_1 = ones(length(bkg_scan))
# upper_mu_for_2 = ones(length(bkg_scan))
# upper_mu_for_3 = ones(length(bkg_scan))
# upper_mu_for_4 = ones(length(bkg_scan))
# for i in 1:length(bkg_scan)
#     n0_limit_list = constructBelt(bkg_scan[i], mu_list, 0.9)
#     upper_mu_dict = getMu2ofBkg(mu_list, n0_limit_list, 0.9)
#     upper_mu_for_0[i] = upper_mu_dict[0]
#     upper_mu_for_1[i] = upper_mu_dict[1]
#     upper_mu_for_2[i] = upper_mu_dict[2]
#     upper_mu_for_3[i] = upper_mu_dict[3]
#     upper_mu_for_4[i] = upper_mu_dict[4]
# end
# fig = plot(bkg_scan, upper_mu_for_0)
# plot!(bkg_scan, upper_mu_for_1)
# plot!(bkg_scan, upper_mu_for_2)
# plot!(bkg_scan, upper_mu_for_3)
# plot!(bkg_scan, upper_mu_for_4)
# savefig(fig, "mu2_bkg.png")
# bkg_scan = 0:0.001:20
# mu_list = 0:0.005:50
# dic = Dict()
# a = getMu2ofBkg(bkg_scan, mu_list, 0.9)
# for i in 0:10
#     dic["$(i)"] = a[i]
# end
# save("mu_bkg_dict_.jld", dic)
# plot(bkg_scan ,a[0])


# showSortedResult(0.9800, 5.0, 60)

# for i in 1:200
#     println(n0_limit_list[i,1], "\t",n0_limit_list[i,2], "\t",mu_list[i])
# end

# Caculate mu_2 as function of bkg.
bkg_scan = 0:0.001:25
mu_2_dic = load("/Users/fangchanghao/Documents/rareEventsAnalysis/analysisCode/meta_data/mu_bkg_dict_0.9.jld")
b = getStepIndex(mu_2_dic["10"], 0.001)
println(mu_2_dic["10"][b])
plot(bkg_scan, mu_2_dic["10"])
scatter!(bkg_scan[b], mu_2_dic["10"][b])