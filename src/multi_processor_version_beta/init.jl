# Initialize this Porgram, RUN THIS FIRST!!!
# Compute some meta_data


push!(LOAD_PATH, ".")
include("MildPathology.jl") # Top level packages

using Plots, JLD
println("Please wait for a while, it may take a long time...")

bkg_upper_limit = 20
mu_upper_limit = 50
cfdent_level = 0.9
nrange = (0, 10)
bkg_scan = 0:0.001:25
mu_list = 0:0.005:50
dic = Dict{String, Array}()
temp = getMu2ofBkg(bkg_scan, mu_list, cfdent_level, nrange = (0, 10))


# re-assign key for storing!! Since JLD need keys to be String.
for i in nrange[1]:nrange[2]
    dic["$(i)"] = temp[i]
end
save("../mu2_bkg_data/mu_bkg_dict_$(mu_upper_limit)_$(bkg_upper_limit)_0.001.jld", dic)
print("Done!")