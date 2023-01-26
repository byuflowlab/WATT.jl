#=
Code used to interact with GXBeam.jl. 
Also some code to run GXBeam in the way that I'm doing it in the aerostructural simulation. 

=#
function retrieve_eulerangles(R)
    return [atan(R[3,2], R[3,3]), asin(-R[3,1]), atan(R[2,1], R[1,1])]
end

function WMPtoangle(c)
    R = GXBeam.wiener_milenkovic(c)'
    return retrieve_eulerangles(R)
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


function update_forces!(distributed_loads, Fx, Fy, Mx, rvec, rgx, rhub, rtip)
    # Fzfit = Linear(vcat(rhub, rvec, rtip), vcat(0, -Fx, 0) ) #I use the loads of the previous time. We want everything to be calculated based off of the previous time step to update the next time step. Then we'll move to that step.  
    # Fyfit = Linear(vcat(rhub, rvec, rtip), vcat(0, Fy, 0))
    # Mxfit = Linear(vcat(rhub, rvec, rtip), vcat(0, Mx, 0))

    Fzfit = Linear(rvec, -Fx) 
    Fyfit = Linear(rvec, Fy)  
    Mxfit = Linear(rvec, Mx)

    f_follower = @SVector zeros(3)
    # m = @SVector zeros(3)
    m_follower = @SVector zeros(3)
    
    for ielem = eachindex(rgx)
        
        f = SVector(0.0, Fyfit(rgx[ielem]), Fzfit(rgx[ielem])) #Using GXBeam's internal damping
        m = SVector(0.0, 0.0, Mxfit(rgx[ielem]))
        distributed_loads[ielem] = GXBeam.DistributedLoads(f, f, m, m, f_follower, f_follower, m_follower, m_follower)
    end
end



function simulate_gxbeam(rvec, rhub, rtip, tvec, azimuth, Fx, Fy, Mx, env::Environment, assembly::GXBeam.Assembly; verbose::Bool=false, speakiter=100, structural_damping::Bool=true, linear::Bool=false, g=9.81)
    # println("Why isn't it printing inside the function.")
    if verbose
        println("Rotors.jl initializing solution...")
    end


    ### Initialization information
    na = length(rvec)
    nt = length(tvec)

    t0 = tvec[1]
    azimuth0 = azimuth[1]

    gxhistory = Array{GXBeam.AssemblyState{eltype(rvec), Vector{GXBeam.PointState{eltype(rvec)}}, Vector{GXBeam.ElementState{eltype(rvec)}}}}(undef, nt)



    ### Extract CCBlade Loads and create a distributed load 
    # Fzfit = Linear(vcat(rhub, rvec, rtip), vcat(0, -Fx[1,:], 0)) 
    # Fyfit = Linear(vcat(rhub, rvec, rtip), vcat(0, Fy[1,:], 0))  #Todo: I think that this is a bit of a problem, because what if the rvec already includes rhub? 
    # Mxfit = Linear(vcat(rhub, rvec, rtip), vcat(0, Mx[1,:], 0))

    Fzfit = Linear(rvec, -Fx[1,:]) 
    Fyfit = Linear(rvec, Fy[1,:])  
    Mxfit = Linear(rvec, Mx[1,:])

    nelem = length(assembly.elements)
    rgx = [assembly.elements[i].x[1] for i in 1:nelem]

    distributed_loads = Dict{Int64, GXBeam.DistributedLoads{eltype(rvec)}}()
    f_follower = @SVector zeros(3)
    # m = @SVector zeros(3)
    m_follower = @SVector zeros(3)

    #Todo. I'm not including gravitational loads here.
    #Todo. I need to add moment loads.  
    for ielem = 1:nelem #Iterate through the elements and apply the distributed load at every element.  
        f = SVector(0.0, Fyfit(rgx[ielem]), Fzfit(rgx[ielem]))
        m = SVector(0.0, 0.0, Mxfit(rgx[ielem])) 
        distributed_loads[ielem] = GXBeam.DistributedLoads(f, f, m, m, f_follower, f_follower, m_follower, m_follower)
    end 



    ### Prepare GXBeam inputs  
    Omega = SVector(0.0, 0.0, -env.RS(t0))
    gravity0 = SVector(g*sin(azimuth0), -g*cos(azimuth0), 0.0)
    prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)) # root section is fixed




    ### GXBeam initial solution #TODO: I might make it an option to initialize from rest, or from steady state at the intial conditions (as current). #Todo: This might be making a large difference, might need to simulate from rest. 
    system, history0, converged = GXBeam.time_domain_analysis(assembly, tvec[1:2]; prescribed_conditions = prescribed_conditions, distributed_loads = distributed_loads, linear = false, angular_velocity = Omega, gravity=gravity0) 

    # gxhistory[1] = AssemblyState(system, assembly; prescribed_conditions = prescribed_conditions)
    gxhistory[1] = history0[end]





    ### Iterate through time 
    for i = 2:nt-1 #Todo: I need to change how I'm doing this for this to work with what I'm doing. 
        t = tvec[i]



        ### Update GXBeam loads #Todo. This is also not accounting for gravitational loads. 
        ##Todo: Does GXBeam account for actual angular movement when keeping the states and what not? Does it integrate and find the new position? -> i.e. changing the gravity. 
        update_forces!(distributed_loads, Fx[i,:], Fy[i,:], Mx[i,:], rvec, rgx, rhub, rtip) #Todo: I'm not confident that the time step that I'm passing this will be correct. 

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
