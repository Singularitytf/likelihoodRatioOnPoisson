using HDF5, JLD
nrange = (0, 20)
CL = 0.99
function loadData(CL, nrange = (0,20))
    global mu1_dic = load("../FCdata/mu1_$(nrange[1])_$(nrange[2])_$(CL).jld")
    global mu2_dic = load("../FCdata/mu2_$(nrange[1])_$(nrange[2])_$(CL).jld")
end

loadData(CL)
fid = h5open("../FCdata/$(nrange[1])_$(nrange[2])_$(CL).h5", "w")
# for i in nrange[1]:nrange[end]
#     println("wirte n0 = ", i)
#     fid["$(i)"] = mu2_dic["$(1)"]
# end
for i in nrange[1]:nrange[end]
    println("wirte n0 = ", i)
    fid["mu1_$(i)"] = mu1_dic["$(i)"]
    fid["mu2_$(i)"] = mu2_dic["$(i)"]
end
fid["bkg"] = mu2_dic["bkg_scan"]

println("Done!!!")

close(fid)