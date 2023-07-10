function nearestto(xvec, x)
    mins = abs.(xvec.-x)
    minval, minidx = findmin(mins)
    minval = xvec[minidx]
    return minval, minidx
end

function getfieldnames(x)
    return fieldnames(typeof(x))
end

function rotate2d(x, y, theta)
    return x*cos(theta) - y*sin(theta), x*sin(theta) + y*cos(theta)
end