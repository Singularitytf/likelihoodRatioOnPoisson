# # testing of Poisson classical interval
# length_of_list = 1000000
# list_test = ones(length_of_list, 2)
# mulist = range(0, 15, length = length_of_list)
# for i in 1:length_of_list
#     func = poissonBG(Float64(mulist[i]), 3.0)
#     list_test[i, 1] = searchIntervalEdge(func, 0.9)[1]
#     list_test[i, 2] = searchIntervalEdge(func, 0.9)[2]
# end
# pic1 = plot(list_test[:,1], mulist, show = true)
# plot!(list_test[:,2], mulist)
# savefig(pic1, "pic1.png")
# upper_list = ones(length_of_list)
# for i in 1:length_of_list
#     func = poissonBG(Float64(mulist[i]), 3.0)
#     upper_list[i] = searchUpperEdge(func, 0.9)
#     # print(searchUpperEdge(func, 0.9))
# end
# pic2 = plot(upper_list, mulist, show = true)
# savefig(pic2, "pic2.png")
# func = poissonBG(0.5, 3.0)
# for i in 0:10
#     print(pdf(func, i), "\n")
# end
push!(LOAD_PATH, ".")
using PoissonLikelihoodRatio
using GR

# set parameters
length_of_list = 10000 # parameter used in [S0556-2821(98)00109-X] is 10,000.
mu_upper_limit = 50
bg = 3.0
cfdent_level = 0.97
n0_upper = 100

mu_list = range(0, mu_upper_limit, length = length_of_list)
n0_limit_list = ones(length_of_list, 2)
for i in 1:length_of_list
    n0_limit_list[i, :] = likelihoodRatio(mu_list[i], bg, cfdent_level, n0_upper)
end
pic1 = plot(n0_limit_list[:,1], mu_list, n0_limit_list[:,2], mu_list)

dict_test = selectMuRegion(mu_list, n0_limit_list)
for i in 0:20
    println(i, "\t",dict_test[i])
end

# for i in 1:200
#     println(n0_limit_list[i,1], "\t",n0_limit_list[i,2], "\t",mu_list[i])
# end