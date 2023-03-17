using DelimitedFiles, MatPrint, OpenFASTsr

of = OpenFASTsr

path = dirname(@__FILE__)
cd(path)

adblade = of.read_adblade("sn5_ADblade.dat", "./")

fullouts = readdlm("./sn5_input.out", skipstart=6)

names = fullouts[1,:]


data = Float64.(fullouts[3:end,:])

outs = Dict(names[i] => data[:,i] for i in eachindex(names))

tvec = outs["Time"]
nt = length(tvec)
n = Int(adblade["NumBlNds"])

mat = zeros(nt, 3*n+1)
# mat = zeros(nt+1, 3*n+1)

# mat[1,1] = nt
# mat[1,2] = n

mat[:,1] = tvec
# mat[2:end,1] = tvec

for i = 1:n
    if i<10
        number = "00$i"
    elseif i<100
        number = "0$i"
    else
        number = "$i"
    end
    mat[:,1+i] = outs["AB1N"*number*"Fx"]
    mat[:,1+n+i] = outs["AB1N"*number*"Fy"]
    mat[:,1+n+n+i] = outs["AB1N"*number*"Mm"]
    # mat[2:end,1+i] = outs["AB1N"*number*"Fx"]
    # mat[2:end,1+n+i] = outs["AB1N"*number*"Fy"]
    # mat[2:end,1+n+n+i] = outs["AB1N"*number*"Mm"]
end

# writedlm("./sn5_output.dat", mat, '\t')
writemat("./sn5_output.dat", [nt 3*n+1]; format="%i", delimiter=" ")
writemat("./sn5_output.dat", mat; append=true, delimiter = " ")

#TODO: Don't forget to add the size to the file after writing to file. 