

#=
The states for this is just the aero states (bem-riso) concatonated with the gxbeam states. Really, all of the inputs that we're going to pass to Differential Equations are just the aero and the structural inputs concatonated together. 

X_as = [x_aero, x_gxbeam]

Y_as = [y_aero, y_gxbeam]

P_as = [p_aero, p_gxbeam]

You now that I'm thinking about it... there might be some repeated parameters that we can easily get rid of. We'll see. 


=#



function create_aerostructuralfun(riso::Riso, bem::BEM, gxmodel::gxbeam, blade::Blade, env::Environment)

    ### Aero preamble - Initialize structs and vectors that don't need to be recreated every iteration. 
    na = length(blade.airfoils)

    rarray = Array{CCBlade.Rotor}(undef, 1)
    afarray = Array{CCBlade.AFType}(undef,1)
    secarray = Array{CCBlade.Section}(undef,1)
    oparray = Array{CCBlade.OperatingPoint}(undef, 1)

    alphavec = collect(-pi:.1:pi)

    ### Structural preamble - Initialize structs and vectors that don't need to be recreated every iteration. 
    ## Create prescribed conditions
    prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0))

    ## Create point masses
    point_masses = Dict(1 => PointMass(0.0, SVector(0.0, 0.0, 0.0), @SMatrix zeros(3,3)))
end

