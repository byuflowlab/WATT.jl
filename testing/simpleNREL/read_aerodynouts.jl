using DelimitedFiles

fullouts = readdlm("./sn5_addriver.1.out", skipstart=6)

names = fullouts[1,:]

# data = readdlm("./simpleNREL/sn5_ADdriver.1.out", skipstart=8)
# data = readdlm("./simpleNREL/sn5_input.out", skipstart=8)

data = Float64.(fullouts[3:end,:])
outs = Dict(names[i] => data[:,i] for i in eachindex(names))


@show outs["AB1N200Phi"].*(pi/180)