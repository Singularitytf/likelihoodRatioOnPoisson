using JLD
using Plots

mu2_dic = load("../meta_data/mu_bkg_dict_0.9.jld")
mu1_dic = 0

function getStepIndex(list, bkg_scan_precision)
    # input an dict mapping n0 to mu_2(b). output its step's index
    
    step_index_list = []

    len = length(list)
    for i in 1:len-1
        if (list[i+1] - list[i])/bkg_scan_precision > 20 # slightly larger than mu precision. If mu_2 is non-continous, it should larger than mu_precision.
            push!(step_index_list, i+1) # push peak index.
        end
    end
    return step_index_list
end

function findIndex(bkg, bkg_scan)
    percision = bkg_scan[2] - bkg_scan[1]
    for i in 1:length(bkg_scan)
        if abs(bkg - bkg_scan[i]) < percision/2 # bkg_scan precision
            index = i
            break
        end
    end
    return index
end


function findModifiedUpperLimit(n0::Int, bkg::Float64, mu2_dic, step_index)
    # Step index is index when mu2 is jumpping.
    bkg_scan = mu2_dic["bkg_scan"]
    mu2_list = mu2_dic[string(n0)]
    bkg_step_points = bkg_scan[step_index]
    mu2_step_points = mu2_list[step_index]
    index::Int = findIndex(bkg, bkg_scan)
    print(index)
    println("The nearest bkg is ", bkg_scan[index], ", whose mu2 is ", mu2_list[index], ".")
    for i in 1:length(step_index)-1
        # find nearest bkg in the dictionary.
        if bkg > bkg_step_points[i] && bkg < bkg_step_points[i+1]
            println(bkg_step_points[i], "\t", bkg_step_points[i+1])
            # if mu is in the dip reigon, modified it!
            if mu2_list[index] <= mu2_step_points[i+1]
                println(mu2_list[index], "\t", mu2_step_points[i+1])
                return mu2_step_points[i+1]
            else
                return mu2_list[index]
            end
        end
    end
end

function findLowerLimit(n0, bkg, mu1_dic)


end

bkg_scan = 0:0.001:25
step_index = getStepIndex(mu2_dic["0"], 0.001)
plot(bkg_scan, mu2_dic["0"])
scatter!(bkg_scan[mu2_step_index], mu2_dic["0"][mu2_step_index])
findModifiedUpperLimit(0, 5.0, bkg_scan, mu2_dic, step_index)