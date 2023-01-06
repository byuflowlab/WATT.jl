using DelimitedFiles


fpcfppc = Float64.(readdlm("/Users/adamcardoza/.julia/dev/Rotors/testing/openfast_fpc_fppc_output.txt")[:,1])

nfile = Int(length(fpcfppc)/2)

fpcvec = zeros(nfile)
fppcvec = zeros(nfile)

let
    IDX = 1
    JDX = 1
    for i = 1:2*nfile
        if iseven(i)
            fppcvec[IDX] = fpcfppc[i]
            IDX += 1
        else
            fpcvec[IDX] = fpcfppc[i]
            JDX += 1
        end
    end
end

fpc = zeros(Int(nfile/19), 19)
fppc = zeros(Int(nfile/19), 19)

for i = 1:Int(nfile/19)
    idx = 19*(i-1)+1:19*i
    fpc[i, :] = fpcvec[idx] 
    fppc[i, :] = fppcvec[idx] 
end
