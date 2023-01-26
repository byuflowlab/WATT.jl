using Plots, LaTeXStrings, Statistics, DelimitedFiles

include("utils.jl")

x = zeros(3)
y = zeros(3)
u = [1, 0, cos(pi/4)*0.66]
v = [0, 1, sin(pi/4)*0.66]

axplt = quiver(x, y, gradient=(u, v), seriescolor=:black, lw=5, leg=false, axis=false, ticks=false)
# annotate!([0.95], [-0.1], [text(L"z_s", 14)]) #Far right label
# annotate!([0.0725], [0.95], [text(L"Y_s", 14)]) #Top label
# annotate!([-0.35], [-0.44], [text(L"X_a", 14)]) #Bottom label
annotate!([1.0], [0.075], [text(L"Z_a", 14)]) #Far right label
annotate!([0.06], [0.95], [text(L"Y_a", 14)]) #Top label
annotate!([0.51], [0.44], [text(L"X_a", 14)]) #Middle label
# display(axplt)
# savefig("/Users/adamcardoza/Library/CloudStorage/Box-Box/research/FLOW/projects/bladeopt/figures/coupling/mycoupling/documentation/referenceframes/Rotors_aerodynamic.png")




###### Read in data
s809cor = readdlm("/Users/adamcardoza/Library/CloudStorage/Box-Box/research/FLOW/experimentaldata/UAE/s809.cor", skipstart=2)

x = s809cor[:,1]
y = s809cor[:,2]

xavg = mean(x)
yavg = mean(y)

xle, xleidx = nearestto(x, 0.0)

xtop = reverse(x[1:xleidx])
ytop = reverse(y[1:xleidx])

ybot = y[xleidx:end]
xbot = x[xleidx:end]

xrot = zero(x)
yrot = zero(y)
theta = 15*pi/180

for i in eachindex(x)
    xrot[i], yrot[i] = rotate2d(x[i], y[i], theta)
end

xshift = -1.0 #xavg
yshift = 0.0 #-yavg
xnew = xrot .+ xshift
ynew = yrot .+ yshift

ynewavg = mean(ynew)

xq = ones(2).*-0.1
yq = ones(2).*-0.1
scale = 0.25
uq = [0, 1].*scale
vq = [-1, 0].*scale

airfoilplt = plot(xnew, ynew, aspectratio=:equal, leg=false, seriescolor=:black, axis=false, ticks=true)
quiver!(xq, yq, gradient=(uq,vq), seriescolor=:black, lw=5)
plot!([-1.2, -1], [-.3, -.15], arrow=:filledtriangle, seriescolor=:navyblue)
hline!([ynewavg], linestyle=:dash, seriescolor=:black)
annotate!([-0.06], [-0.3], [text(L"u", 14)])
annotate!([0.1], [-0.15], [text(L"v", 14)])
display(airfoilplt)