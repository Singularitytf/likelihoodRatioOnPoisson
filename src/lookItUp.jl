using JLD
using Plots

function getStepIndex(list)
    # input an dict mapping n0 to mu_2(b). output its step's index
    
    step_index_list = []
    bkg_scan_precision = abs(list[2] - list[1])
    # println(bkg_scan_precision)

    len = length(list)
    for i in 1:len-1
        if (list[i+1] - list[i])/bkg_scan_precision > 15 # slightly larger than mu precision. If mu_2 is non-continous, it should larger than mu_precision.
            push!(step_index_list, i+1) # push peak index.
        end
    end
    return step_index_list
end

function findIndex(bkg, bkg_scan)
    percision = bkg_scan[2] - bkg_scan[1]
    index::Int = 0
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
    println("--------------------------")
    println("The nearest bkg is ", bkg_scan[index], ", whose mu2 is ", mu2_list[index], ".")
    println("--------------------------")
    for i in 1:length(step_index)-1
        # find nearest bkg in the dictionary.
        if bkg > bkg_step_points[i] && bkg < bkg_step_points[i+1] # bkg in the abnormal part.
            println("The bkg interval is:")
            println(bkg_step_points[i], "\t", bkg_step_points[i+1])
            # if mu is in the dip reigon, modified it!
            if mu2_list[index] <= mu2_step_points[i+1]
                println("The mu interval is:")
                println(mu2_list[index], "\t", mu2_step_points[i+1])
                return mu2_step_points[i+1]
            else
                return mu2_list[index]
            end
        elseif bkg < bkg_step_points[1] # bkg in the left region.
            return mu2_list[index]
        elseif bkg > last(bkg_step_points)
            println("Warning: Need to know the next step point!!! Hint: enlarge bkg_scan region")
        end
    end
end

function findLowerLimit(n0, bkg, mu1_dic)
    bkg_scan = mu1_dic["bkg_scan"]
    mu1_list = mu1_dic[string(n0)]
    index::Int = findIndex(bkg, bkg_scan)
    println("--------------------------")
    println("The nearest bkg is ", bkg_scan[index], ", whose μ_1 is ", mu1_list[index], ".")
    println("--------------------------")
    return mu1_list[index]
end

function findInterval(n0, bkg)
    step_index = getStepIndex(mu2_dic[string(n0)])
    mu2 = findModifiedUpperLimit(n0, bkg, mu2_dic, step_index)
    mu1 = findLowerLimit(n0, bkg, mu1_dic)
    return [mu1, mu2]
end

# Loading the init data.
global mu1_dic = load("mu_bkg_data/mu1_bkg_dict_50_20_0.001.jld")
global mu2_dic = load("mu_bkg_data/mu2_bkg_dict_50_20_0.001.jld")
bkg_scan = mu2_dic["bkg_scan"]
step_index = getStepIndex(mu2_dic["2"])
plot(bkg_scan, mu2_dic["2"])
scatter!(bkg_scan[step_index], mu2_dic["2"][step_index])
print("μ interval is:", findInterval(3, 10.0))