dic1 = Dict{Int, Int}()
for i in 1:10
    dic1[i] = 10*i
end

for i in 10:15
    try
        dic1[i]
    catch e
        if isa(e, KeyError)
            dic1[i] = i
        end
        continue
    end
end
