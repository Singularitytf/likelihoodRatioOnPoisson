# Initialize this Porgram, RUN THIS FIRST!!!
# Compute some meta_data


push!(LOAD_PATH, ".")
using MildPathology # Top level packages
using Plots, JLD
print("Please wait for a while, it may take a long time...")
bkg_upper_limit = 20
mu_upper_limit = 50
bkg_scan = 0:0.001:20
mu_list = 0:0.005:50
dic = Dict{String, Array}()
temp = getMu2ofBkg(bkg_scan, mu_list, 0.9)
# re-assign key for storing!! Since JLD need keys to be String.
for i in 0:10
    dic["$(i)"] = temp[i]
end
save("./mu2_bkg_data/mu_bkg_dict_$(mu_upper_limit)_$(bkg_upper_limit)_0.001.jld", dic)