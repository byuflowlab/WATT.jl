#=
Determine why the dynamic stall states at the tip of the simpleNREL turbine in a decoupled aerostructural solution are not changing. 
Adam Cardoza 5/9/23
=#

using DynamicStallModels, OpenFASTsr, Plots

of = OpenFASTsr
dsm = DynamicStallModels

localpath = @__DIR__
cd(localpath)

aft = of.read_airfoilinput("./Airfoils/NACA64_A17.dat") 

c = 1.4190122

airfoils = Vector{dsm.Airfoil}(undef, 1)
airfoils[1] = of.make_dsairfoil(aft, c)


U = 76.11175070589553
alpha = 0.10148806382354443

tvec = collect(0:0.001:0.01)
nt = length(tvec)
Uvec = ones(nt).*U
aoavec = ones(nt).*alpha


states, loads = dsm.solve_indicial(airfoils, tvec, Uvec, aoavec)


#=
5/9/23
Well... I'm fairly confident in the ability of the solve_indicial function to solve correctly, as it compares quite well against OpenFAST when fed the same inflow velocities and angles that OpenFAST calculates. Thus somehow, this angle of attack is somehow at steady state. I don't understand how that works. I mean I could pop open the aero-only solve and see what the loading across the blade looks like (both to compare smoothness and what the tip does across time). 
=#