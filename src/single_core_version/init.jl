# Initialize this Porgram, RUN THIS FIRST!!!
# Compute some meta_data


push!(LOAD_PATH, ".")
using MildPathology # Top level packages
using JLD
println("Please wait for a while, it may take several hours...")
bkg_upper_limit = 20
mu_upper_limit = 50
bkg_precision = 0.1
bkg_scan = 0:bkg_precision:25
mu_list = 0:0.001:50
CL = 0.9
nrange = (0, 10)
mu2 = Dict{String, Array}()
mu1 = Dict{String, Array}()
upper_dic, lower_dic = getMuofBkg(bkg_scan, mu_list, CL, nrange = (0, 10))
# re-assign key for storing!! Since JLD need keys to be String.
for i in nrange[1]:nrange[2]
    mu2["$(i)"] = upper_dic[i]
    mu1["$(i)"] = lower_dic[i]
end
mu2["bkg_scan"] = bkg_scan
mu1["bkg_scan"] = bkg_scan
save("../mu_bkg_data/mu2_bkg_dict_$(mu_upper_limit)_$(bkg_upper_limit)_$(bkg_precision).jld", mu2)
save("../mu_bkg_data/mu1_bkg_dict_$(mu_upper_limit)_$(bkg_upper_limit)_$(bkg_precision).jld", mu1)
print("Done!")