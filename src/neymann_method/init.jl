push!(LOAD_PATH, ".")
using method, Plots, JLD

println("Please wait for a while, it may take several hours...")
bkg_upper_limit = 20
mu_upper_limit = 50
bkg_precision = 0.01
bkg_scan = 0:bkg_precision:25
mu_list = 0:0.01:50
CL = 0.9
nrange = (0, 20)
mu2 = Dict{String, Array}()
mu1 = Dict{String, Array}()
upper_dic, lower_dic = getMuofBkg(bkg_scan, mu_list, CL, nrange = nrange)
# re-assign key for storing!! Since JLD need keys to be String.
for i in nrange[1]:nrange[2]
    mu2["$(i)"] = upper_dic[i]
    mu1["$(i)"] = lower_dic[i]
end
mu2["bkg_scan"] = bkg_scan
mu1["bkg_scan"] = bkg_scan
try
    mkdir("../NMdata/")
    save("../NMdata/mu2_$(nrange[1])_$(nrange[2])_$(CL).jld", mu2)
    save("../NMdata/mu1_$(nrange[1])_$(nrange[2])_$(CL).jld", mu1)
catch e
    save("../../NMdata/mu2_$(nrange[1])_$(nrange[2])_$(CL).jld", mu2)
    save("../../NMdata/mu1_$(nrange[1])_$(nrange[2])_$(CL).jld", mu1)
end
print("Done!")