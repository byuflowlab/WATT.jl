#=
Code used to interact with GXBeam.jl. 
Also some code to run GXBeam in the way that I'm doing it in the aerostructural simulation. 

=#
function retrieve_eulerangles(R)
    return SVector{3}([atan(R[3,2], R[3,3]), asin(-R[3,1]), atan(R[2,1], R[1,1])])
end

function WMPtoangle(c)
    scaling = GXBeam.rotation_parameter_scaling(c)
    R = GXBeam.wiener_milenkovic(scaling*c)'
    return retrieve_eulerangles(R)
end

function get_bladelength_vector(assembly::GXBeam.Assembly)
    ns = length(assembly.points)
    rgx = zeros(ns)
    rgx[1] = norm(assembly.points[1])
    for i = 2:ns
        delr = assembly.points[i] - assembly.points[i-1]
        rgx[i] = rgx[i-1] + norm(delr)
    end
    return rgx
end

function update_assembly(assembly; compliance=nothing, stiffness=nothing)

    points = assembly.points

    start = assembly.start
    stop = assembly.stop

    ### Update Elements
    if !isnothing(compliance) || !isnothing(stiffness)
        if !isnothing(stiffness)
            if isa(stiffness, Vector)
                if length(stiffness) != length(assembly.elements)
                    error("The number of elements in the stiffness vector does not match the number of elements in the assembly.")
                end
                compliance = [SMatrix{6,6}(inv(stiffness[i])) for i = eachindex(assembly.elements)]
            else
                compliance = [inv(stiffness) for i = eachindex(assembly.elements)]
            end
        else
            if isa(compliance, Vector)
                if length(compliance) != length(assembly.elements)
                    error("The number of elements in the compliance vector does not match the number of elements in the assembly.")
                end
            else
                compliance = [SMatrix{6,6}(compliance) for i = eachindex(assembly.elements)]
            end
        end
    end

    elements = [GXBeam.Element(assembly.elements[i].L, assembly.elements[i].x, compliance[i], assembly.elements[i].mass, assembly.elements[i].Cab, assembly.elements[i].mu) for i = eachindex(assembly.elements)]

    

    return GXBeam.Assembly(points, start, stop, elements)

end

"""
    interpolate_matrix_symmetric(f1, f2, Kmat)
Interpolate a set of symmetric matrices from one vector of fractions to any number of new fractions.

### Inputs
- f1::Vector{Float64}: The fractional locations of the original stiffness matrices.
- f2::Vector{Float64}: The fractional locations of the new stiffness matrices.
- Kmat::Array{Float64, 3}: The stiffness matrices at the original locations.

### Outputs
- Kfit::Array{Float64, 3}: The stiffness matrices at the new locations.
"""
function interpolate_matrix_symmetric(f1, f2, Kmat; fit=Linear)
    if length(f1)!=size(Kmat, 3)
        error("interpolate_stiffness(): The number of f1 points and stiffness matrices must be the same.")
    end

    K11fit = fit(f1, Kmat[1,1,:])
    K12fit = fit(f1, Kmat[1,2,:])
    K13fit = fit(f1, Kmat[1,3,:])
    K14fit = fit(f1, Kmat[1,4,:])
    K15fit = fit(f1, Kmat[1,5,:])
    K16fit = fit(f1, Kmat[1,6,:])

    K22fit = fit(f1, Kmat[2,2,:])
    K23fit = fit(f1, Kmat[2,3,:])
    K24fit = fit(f1, Kmat[2,4,:])
    K25fit = fit(f1, Kmat[2,5,:])
    K26fit = fit(f1, Kmat[2,6,:])

    K33fit = fit(f1, Kmat[3,3,:])
    K34fit = fit(f1, Kmat[3,4,:])
    K35fit = fit(f1, Kmat[3,5,:])
    K36fit = fit(f1, Kmat[3,6,:])

    K44fit = fit(f1, Kmat[4,4,:])
    K45fit = fit(f1, Kmat[4,5,:])
    K46fit = fit(f1, Kmat[4,6,:])

    K55fit = fit(f1, Kmat[5,5,:])
    K56fit = fit(f1, Kmat[5,6,:])

    K66fit = fit(f1, Kmat[6,6,:])

    Kfit = zeros(6,6,length(f2))
    for i in eachindex(f2)
        Kfit[1,1,i] = K11fit(f2[i])
        Kfit[1,2,i] = K12fit(f2[i])
        Kfit[1,3,i] = K13fit(f2[i])
        Kfit[1,4,i] = K14fit(f2[i])
        Kfit[1,5,i] = K15fit(f2[i])
        Kfit[1,6,i] = K16fit(f2[i])

        Kfit[2,1,i] = K12fit(f2[i])
        Kfit[2,2,i] = K22fit(f2[i])
        Kfit[2,3,i] = K23fit(f2[i])
        Kfit[2,4,i] = K24fit(f2[i])
        Kfit[2,5,i] = K25fit(f2[i])
        Kfit[2,6,i] = K26fit(f2[i])

        Kfit[3,1,i] = K13fit(f2[i])
        Kfit[3,2,i] = K23fit(f2[i])
        Kfit[3,3,i] = K33fit(f2[i])
        Kfit[3,4,i] = K34fit(f2[i])
        Kfit[3,5,i] = K35fit(f2[i])
        Kfit[3,6,i] = K36fit(f2[i])

        Kfit[4,1,i] = K14fit(f2[i])
        Kfit[4,2,i] = K24fit(f2[i])
        Kfit[4,3,i] = K34fit(f2[i])
        Kfit[4,4,i] = K44fit(f2[i])
        Kfit[4,5,i] = K45fit(f2[i])
        Kfit[4,6,i] = K46fit(f2[i])

        Kfit[5,1,i] = K15fit(f2[i])
        Kfit[5,2,i] = K25fit(f2[i])
        Kfit[5,3,i] = K35fit(f2[i])
        Kfit[5,4,i] = K45fit(f2[i])
        Kfit[5,5,i] = K55fit(f2[i])
        Kfit[5,6,i] = K56fit(f2[i])

        Kfit[6,1,i] = K16fit(f2[i])
        Kfit[6,2,i] = K26fit(f2[i])
        Kfit[6,3,i] = K36fit(f2[i])
        Kfit[6,4,i] = K46fit(f2[i])
        Kfit[6,5,i] = K56fit(f2[i])
        Kfit[6,6,i] = K66fit(f2[i])
    end

    return Kfit
end

function pane_assembly(assembly; ne=nothing, verbose::Bool=false, fit=Linear)
    if isnothing(ne)
        return assembly
    elseif ne<1
        if verbose
            @warn("The number of elements must be at least one. Increasing number of elements to 1.")
        end
        ne = 1
    end


    ### Extract the current assembly points
    np_cur = length(assembly.points)
    ne_cur = np_cur - 1

    points = zeros(np_cur, 3)
    for i = 1:np_cur
        points[i,:] = collect(assembly.points[i].data)
    end

    ### Calculate the element lengths (so I can calculate the assembly length)
    lvec = zeros(ne_cur)
    for i = 2:np_cur
        p1 = view(points, i, :)
        p2 = view(points, i-1, :)
        delp = p2.-p1
        
        lvec[i-1] = sqrt(sum(delp.^2))
    end

    L = sum(lvec) #Calculate the assembly length

    ### Calculate non-dimensional position of points along beam
    svec = zeros(np_cur)
    for i = 2:np_cur
        svec[i] = sum(lvec[1:i-1])/L
    end

    Xfit = fit(svec, points[:,1])
    Yfit = fit(svec, points[:,2])
    Zfit = fit(svec, points[:,3])

    # elvec = [assembly.elements[i].L for i = 1:ne_cur]
    sevec = [(sum(lvec[1:i-1])+(lvec[i]/2))/L for i in 1:ne_cur]






    ### New Beam
    np = ne + 1
    svec_new = collect(range(0.0, 1.0, np))
    points_new = zeros(np, 3)

    points_new[:,1] = Xfit.(svec_new)
    points_new[:,2] = Yfit.(svec_new)
    points_new[:,3] = Zfit.(svec_new)

    pointsvec = [SVector{3}(points_new[i,:]) for i = 1:np]
    newlvec = (L/ne).*ones(ne)

    x_elements = [SVector{3}((points_new[i,:].+points_new[i+1,:])./2) for i in 1:ne]

    elvec_new = [(L/(2*ne)) + (i-1)*(L/ne) for i = 1:ne]

    sevec_new = [elvec_new[i]/L for i in 1:ne]

    return sevec_new
end

function plotpoints(points; xdim = true, ydim = true, zdim=true)
    n = length(points)
    x = [points[i][1] for i in 1:n]
    y = [points[i][2] for i in 1:n]
    z = [points[i][3] for i in 1:n]

    if iszero(x)
        xdim = false
        xplt = y
        yplt = z
        xlab = "Y distance"
        ylab = "Z distance"
    elseif iszero(y)
        ydim = false
        xplt = x
        yplt = z
        xlab = "X distance"
        ylab = "Z distance"
    elseif iszero(z)
        zdim = false
        xplt = x
        yplt = y
        xlab = "X distance"
        ylab = "Y distance"
    end

    if xdim&&ydim&&zdim #Check if the plot is 3 dimensional
        plt = scatter(x, y, z, xaxis="X distance", yaxis="Y distance", zaxis="Z distance", legend=false, aspectratio=:equal)
    else
        plt = scatter(xplt, yplt, xaxis=xlab, yaxis=ylab, legend=false, aspectratio=:equal)
    end
    return plt
end

function plotassembly(assembly; xdim = true, ydim = true, zdim=true)
    points = cat(assembly.points...,dims=2)'

    ### Get element points
    elements = zeros(length(assembly.elements), 3)
    for i = 1:length(assembly.elements)
        elements[i,:] = assembly.elements[i].x
    end
    
    if iszero(points[:,1]) #X direction is zero
        xdim = false
        xplt = points[:,2]
        yplt = points[:,3]
        xlab = "Y distance"
        ylab = "Z distance"
        xmplt = elements[:,2]
        ymplt = elements[:,3]
    elseif iszero(points[:,2]) #Y direction is zero
        ydim = false
        xplt = points[:,1]
        yplt = points[:,3]
        xlab = "X distance"
        ylab = "Z distance"
        xmplt = elements[:,1]
        ymplt = elements[:,3]
    elseif iszero(points[:,3]) #Z direction is zero
        zdim = false
        xplt = points[:,1]
        yplt = points[:,2]
        xlab = "X distance"
        ylab = "Y distance"
        xmplt = elements[:,1]
        ymplt = elements[:,2]
    end



    if xdim&&ydim&&zdim #Check if the plot is 3 dimensional
        plt = plot(points[:,1], points[:,2], points[:,3], linewidth=4, linecolor=:black, xaxis="X distance", yaxis="Y distance", zaxis="Z distance", legend=false, aspectratio=:equal)
        scatter!(points[:,1], points[:,2], points[:,3])
        scatter!(elements[:,1], elements[:,2], elments[:,3], markershape=:cross)
    else
        plt = plot(xplt, yplt, linewidth=4, linecolor=:black, xaxis=xlab, yaxis=ylab, legend=false, aspectratio=:equal)
        scatter!(xplt, yplt, markersize=6)
        scatter!(xmplt, ymplt, markershape=:diamond, markersize=4)
    end
    return plt

end

function gxbeam_initial_conditions(env::Environment, assembly, prescribed_conditions, distributed_loads, t0, azimuth0, g, structural_damping, linear, flag)

    Omega0 = SVector(0.0, 0.0, -env.RS(t0))
    gravity0 = SVector(-g*cos(azimuth0), -g*sin(azimuth0), 0.0)

    if flag==:steady
        @warn("Steady initialization not yet prepared, starting from no load no deflection. ")

        system, history0, converged = GXBeam.time_domain_analysis(assembly, [t0]; prescribed_conditions, distributed_loads, angular_velocity = Omega0, gravity=gravity0, steady_state=false, structural_damping, linear) 

    elseif flag==:spinning
        system, history0, converged = GXBeam.time_domain_analysis(assembly, [t0]; prescribed_conditions = prescribed_conditions, angular_velocity = Omega0, gravity=gravity0, steady_state=true, structural_damping, linear)

    else #No load, no deflection initialization. 
        system, history0, converged = GXBeam.time_domain_analysis(assembly, [t0]; prescribed_conditions = prescribed_conditions, distributed_loads = distributed_loads, angular_velocity = Omega0, gravity=gravity0, steady_state=false, structural_damping, linear) 
    end


    gxstate = history0[end]
    return gxstate, system
end

function update_forces!(distributed_loads, Fx, Fy, Mx, blade, assembly; fit=DS.Linear)

    #Todo: I think that this is a bit of a problem, because what if the rvec already includes rhub? (problem from before that needs to be resolved, see next Todo statement. ) #Todo: I need to nail down behavior outside of the aero node regions.
    Fzfit = fit(blade.r, -Fx)  
    Fyfit = fit(blade.r, Fy) 
    Mxfit = fit(blade.r, Mx)


    for ielem = eachindex(assembly.elements)
        # r1 = assembly.points[ielem][1] #Todo,: I want a vector of just lengths, of the points. Not just the X distance. 
        # r2 = assembly.points[ielem+1][1]
        r1 = norm(assembly.points[ielem])
        r2 = norm(assembly.points[ielem+1])
        distributed_loads[ielem] = GXBeam.DistributedLoads(assembly, ielem; fy_follower = (s) -> Fyfit(s), fz_follower = (s) -> Fzfit(s), s1=r1, s2=r2) #, mx = (s) -> Mxfit(s)
        # distributed_loads[ielem] = GXBeam.DistributedLoads(assembly, ielem; fy = (s) -> Fyfit(s), fz = (s) -> Fzfit(s), s1=r1, s2=r2) #, mx = (s) -> Mxfit(s)
        #Todo: There is a slight problem here, if changing from follower loads to dead loads does absolutely nothing... then I'm not sure that what Taylor says they are doing is what they are actually doing. I need to look into that behavior. -> He applies the rotation matrix to the follower loads... And it looks like he does it correctly, or rather 
    end

end



function simulate_gxbeam(rvec, rhub, rtip, tvec, azimuth, Fx, Fy, Mx, env::Environment, assembly::GXBeam.Assembly; verbose::Bool=false, speakiter=100, structural_damping::Bool=true, linear::Bool=false, g=9.81)
    # println("Why isn't it printing inside the function.")
    if verbose
        println("Rotors.jl initializing solution...")
    end


    ### Initialization information
    nt = length(tvec)

    t0 = tvec[1]
    azimuth0 = azimuth[1]

    gxhistory = Array{GXBeam.AssemblyState{eltype(rvec), Vector{GXBeam.PointState{eltype(rvec)}}, Vector{GXBeam.ElementState{eltype(rvec)}}}}(undef, nt)



    ### Extract CCBlade Loads and create a distributed load 
    # Fzfit = Linear(vcat(rhub, rvec, rtip), vcat(0, -Fx[1,:], 0)) 
    # Fyfit = Linear(vcat(rhub, rvec, rtip), vcat(0, Fy[1,:], 0))  #Todo: I think that this is a bit of a problem, because what if the rvec already includes rhub? 
    # Mxfit = Linear(vcat(rhub, rvec, rtip), vcat(0, Mx[1,:], 0))

    # Fzfit = Linear(rvec, -Fx[1,:]) 
    # Fyfit = Linear(rvec, Fy[1,:])  
    # Mxfit = Linear(rvec, -Mx[1,:])

    nelem = length(assembly.elements)
    rgx = [assembly.elements[i].x[1] for i in 1:nelem]
    rgxp = [assembly.points[i][1] for i in eachindex(assembly.points)]

    distributed_loads = Dict{Int64, GXBeam.DistributedLoads{eltype(rvec)}}()
    update_forces!(distributed_loads, Fx[1,:], Fy[1,:], Mx[1,:], rvec, assembly)
    # f_follower = @SVector zeros(3) #Todo: Should I be applying a follower force instead of a dead load?
    # # m = @SVector zeros(3)
    # m_follower = @SVector zeros(3)

    # #Todo. I'm not including gravitational loads here.
    # #Todo. I need to add moment loads.  
    # for ielem = 1:nelem #Iterate through the elements and apply the distributed load at every element.  
    #     f = SVector(0.0, Fyfit(rgx[ielem]), Fzfit(rgx[ielem]))
    #     m = SVector(Mxfit(rgx[ielem]), 0.0, 0.0) 
    #     distributed_loads[ielem] = GXBeam.DistributedLoads(f, f, m, m, f_follower, f_follower, m_follower, m_follower)
    # end 



    ### Prepare GXBeam inputs  
    Omega = SVector(0.0, 0.0, -env.RS(t0)) #Todo: This might be spinning the wrong way. 
    gravity0 = SVector(g*sin(azimuth0), -g*cos(azimuth0), 0.0)
    prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)) # root section is fixed
    # V0 = [SVector(0.0, -rgxp[i]*env.RS(t0), 0.0) for i in eachindex(assembly.points)]


    ### GXBeam initial solution #TODO: I might make it an option to initialize from rest, or from steady state at the intial conditions (as current). 
    system, history0, converged = GXBeam.time_domain_analysis(assembly, tvec[1:1]; prescribed_conditions = prescribed_conditions, distributed_loads = distributed_loads, angular_velocity = Omega, gravity=gravity0, steady=true) #Initialize from undeflected, everything moving. -> Taylor fixed GXBema so that adding the angular_velocity will update the velocity of everything in the global frame when initializing. 

    gxhistory[1] = history0[end]





    ### Iterate through time 
    for i = 1:nt-1 #Todo: I need to change how I'm doing this for this to work with what I'm doing. 
        t = tvec[i]



        ### Update GXBeam loads #Todo. This is also not accounting for gravitational loads. 
        ##Todo: Does GXBeam account for actual angular movement when keeping the states and what not? Does it integrate and find the new position? -> i.e. changing the gravity. 
        update_forces!(distributed_loads, Fx[i,:], Fy[i,:], Mx[i,:], rvec, assembly) #Todo: I'm not confident that the time step that I'm passing this will be correct. 

        Omega = SVector(0.0, 0.0, -env.RS(tvec[i]))
        gravity = SVector(g*sin(azimuth[i]), -g*cos(azimuth[i]), 0.0)


        ### Solve GXBeam for time step
        system, localhistory, converged = GXBeam.time_domain_analysis!(system, assembly, tvec[i:i+1]; prescribed_conditions=prescribed_conditions, distributed_loads, linear, angular_velocity = Omega, reset_state=false, initialize=false, structural_damping, gravity) #TODO: Add gravity. 



        ### Extract GXBeam outputs
        gxhistory[i] = localhistory[end] 



        if verbose & (mod(i, speakiter)==0)
            println("")
            println("Simulation time: ", t)
            # @show Vxvec[end], Vyvec[end] 
            # @show aeroV[end]
            # @show azimuth[i], localhistory[end].elements[end].V
            # @show maximum(N[i,:]), maximum(T[i,:])
        end
    end

    gxhistory[end] = gxhistory[end-1] #Hack: 

    return gxhistory
end

function steady_simulate_gxbeam(rvec, azimuth, Fx, Fy, Mx, env::Environment, assembly::GXBeam.Assembly; verbose::Bool=false, speakiter=100, structural_damping::Bool=true, linear::Bool=false, g=9.81)
    # println("Why isn't it printing inside the function.")
    if verbose
        println("Rotors.jl initializing solution...")
    end


    ### Initialization information

    t0 = 0.0
    azimuth0 = azimuth[1]


    ### Extract CCBlade Loads and create a distributed load 
    nelem = length(assembly.elements)
    rgx = [assembly.elements[i].x[1] for i in 1:nelem]
    rgxp = [assembly.points[i][1] for i in eachindex(assembly.points)]

    distributed_loads = Dict{Int64, GXBeam.DistributedLoads{eltype(rvec)}}()
    update_forces!(distributed_loads, Fx[1,:], Fy[1,:], Mx[1,:], rvec, assembly)
    



    ### Prepare GXBeam inputs  
    Omega = SVector(0.0, 0.0, -env.RS(t0)) 
    gravity0 = SVector(g*sin(azimuth0), -g*cos(azimuth0), 0.0)
    prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)) # root section is fixed
    


    ### GXBeam initial solution  
    system, converged = GXBeam.steady_state_analysis(assembly; prescribed_conditions = prescribed_conditions, distributed_loads = distributed_loads, angular_velocity = Omega, gravity=gravity0) 

    if !converged
        println("Yo, GXBeam didn't converge.")
    end

    state = AssemblyState(system, assembly; prescribed_conditions)



    return state
end


function get_blade_weight(assembly::GXBeam.Assembly)

    n = length(assembly.elements)

    rvec = zeros(n)
    masses = zeros(n)

    for i = 1:n
        rvec[i] = norm(assembly.elements[i].x)
        masses[i] = assembly.elements[i].mass[1,1]
    end

    return trapz(rvec, masses)
end