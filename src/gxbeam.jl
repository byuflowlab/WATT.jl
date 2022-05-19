
#=
############## States ###################

    x = [x_p, ..., x_e, ...] 

where x_p is a vector of states for a given point, and x_e is a vector states for a given gxbeam element. There is a x_p for every point in our beam, and a x_e for every gxbeam element. 

x_p has two options, depending if we are constraining displacement or not. Regardless of if we are constraining displacement, there are six elements of x_p; which implies that the total number of elements in x from the point states are 6*n_p (where n_p is the number of points).

    If we're constraining displacement, then the vector of states is:
        x_p = [F_x, F_y, F_z, M_x, M_y, M_z]. 
    These are the forces and moments required to keep the point stationary. 

    If we're not constraining displacement (it has that degree of freedom), then the vector of states is 
    x_p = [delta_x, delta_y, delta_z, theta_x, theta_y, theta_z], 
    where delta and theta are linear and angular displacements respectively. 


x_e has only one option with 18 elements that are the same regardless of displacement constraints. These states are: 

        x_e = [delta_x, delta_y, delta_z, theta_x, theta_y, theta_z, F_x, F_y, F_z, M_x, M_y, M_z, V_x, V_y, V_z, omega_x, omega_y, omega_z], 

    where delta and theta are linear and angular displacements of the gxbeam element, F and M are the forces and moments on the gxbeam element, V and omega are the linear and angular velocities. I'm pretty sure that I never alter any of these, but I'm not 100%. I think that the applied forces come in through the inputs. There are 18 elements in x_e, so there are 18*n_e elements in x (where n_e is the number of gxbeam elements). 
    I know that the magnitudes of the forces and moments are defined at the center of the beam element, I'm not 100% certain if the same is true for the linear displacements and velocities. I'm pretty sure that it doesn't matter for the angular displacements and velocities, but not 100%. 

- Note that the states are all actually just vectors (columns). 

########### Inputs ##############

    y = [y_pl, ..., y_dl, ..., y_pm, ..., y_s]

where y_pl is the point loadings, y_dl is the distributed loads, y_pm are any point masses, and y_s is general system inputs. There is a y_pl for every point, a y_dl for every element, a y_pm for every element, and only one y_s. 




#Todo: Shouldn't point masses be a parameter, not an input? -> I'm not sure that it will matter


=#
# """
#    gxbeam(displacement)
# A struct to hold information about the gxbeam model.

# ### Inputs
# - displacement::Array{Float64, 2} - a matrix that holds fixity/DoF information for all the points in the beam. The matrix will be nx6 where n is the number of points, and the row represents the x, y, z, theta_x, theta_y, theta_z degree of freedom (DoF). If the value is true, then the point has that degree of freedom. To iterate, if the value is false, then the point is fixed in that direction. For example, if the first row of the matrix read [false, false, true, false, false, false], then the start of the beam would be fixed in the x and y directions, but free to translate in the z. It would not be allowed to rotate at all. 
# """
# struct gxbeam{TF}
#     displacement::Array{TF, 2} 
# end

"""
    gxbeam(displacement)
A struct to hold information about the gxbeam model. 

### Inputs
- L::Float64 - The length of the beam
- n::Int64 - Number of beam elements
"""
struct gxbeam{TI}
    n::TI #Do I even need this? 
end

### residual function
#=
- I don't want to be forming an assembly every iteration... Plus... that should remain constant across the entire simulation. 

- Do I want the assembly passed in? 
    - Well... I guess it depends what goes into the assembly. Is it dependent on any parameters? Very probably. 
        -> Yes, it is dependent on parameters, which implies that I need to create the assembly every iteration. 

- What goes in the prescribed conditions? 
    - Fixity conditions
    - "..., prescribed conditions are forces and/or displacement boundary conditions applied to points." - GXBeam documentation. 
    - "One instance of `PrescribedConditions` must be created for every point on which prescribed conditions are applied." - GXBeam documentation. 
    - It appears to hold applied loads (force and moment) and follower loads at each point. 
    - Not point masses, or distributed loads

    - > At this point, I don't see any other need for prescribed conditions other than the fixity conditions. And in reality... this doesn't need to be customizable... will there be any other fixity constraints on a wind turbine? I don't think so. -> TODO: I can remove displacement from the struct, point loadings from inputs, and potentially from the states... although... I might still need the loadings for the residual function. Definitely don't need them in parameters if they show up there. 

- What goes in distributed loads? 
    - directional and follower forces and moments -> Todo: Should I put the aerodynamic forces as follower or directional forces. I think they should be follower. -> Taylor says that they should directional (or "dead loads") because at every iteration you're going to recalculate them anyway. -> I want to see what DG thinks. 
    - Note: All loads are defined in the body frame (not the local frame). Follower loads will be changed as the beam deflects, but they are still defined in the body frame. 
    - They are integrated using 4 point Gauss-Legendre quadrature. 
    - The load is a function of s, and the start and end of the beam can be specified using keywords s1 and s2. #TODO: I could have a varying force load on the beam.. or I could do sections. I wonder which would be better. I should probably start with a constant load on each element, for simplicity. 

    - Outcome is, I think I'm going to pass in the distributed loads. So at this point, I'm not sure that I really need a whole new function to pass stuff into GXBeam... it seems like I'd just be adding a little overhead. 

- What are these irow_point and friends?
    - irow_point: row index of first equilibrium equation for each point
    - irow_elem: row index of first equation for just this beam element
    - irow_elem1: row index of first equation for the left side of each beam
    - irow_elem2: row index of first equation for the right side of each beam
    - icol_point: column index of first state variable for each point
    - icol_elem: column index of first state variable for each beam element

    - All of these bad boys come from the function GXBeam.system_indices. They should be constant. 

- What are these x0 and friends? 
    - x0 - body frame linear origin
    - v0 - body frame linear velocity
    - omega0 - body frame angular velocity
    - a0 - body frame linear acceleration
    - alpha0 - body frame angular acceleration
=#
function gxbeam_residual(residuals, dx, x, assembly, prescribed_conditions) #Not going to use this function, going to use the GXBeam function. #TODO: Clean this function away. 


    GXBeam.dynamic_system_residual!(residuals, dx, x, assembly, prescribed_conditions, distributed_loads, point_masses, gvec,
    force_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem,
    x0, v0, omega0, a0, alpha0)
end

create_gxbeam_point(p) = SVector(p[1], p[2], p[3])

function create_gxbeam_element(p, points, start, stop)
    # separate element parameters
    e1x, e1y, e1z, e2x, e2y, e2z, e3x, e3y, e3z, C11, C12, C13, C14, C15, C16, C22, C23, C24, C25, C26, C33, C34, C35, C36, C44, C45, C46, C55, C56, C66, mu, xm2, xm3, i22, i33, i23 = p  
        # e - reference frame orientation vector (1,2,3 correspond to x,y,z respectively)
    
    # element length
    DeltaL = norm(points[stop] - points[start])
    
    # element location
    x = (points[start] + points[stop])/2

    # element compliance matrix
    C = @SMatrix [
        C11 C12 C13 C14 C15 C16;
        C12 C22 C23 C24 C25 C26;
        C13 C23 C33 C34 C35 C36;
        C14 C24 C34 C44 C45 C46;
        C15 C25 C35 C45 C55 C56;
        C16 C26 C36 C46 C56 C66]

    # element mass matrix
    mass = @SMatrix [
            mu       0     0       0 mu*xm3 -mu*xm2;
            0       mu     0  -mu*xm3     0      0;
            0       0     mu   mu*xm2     0      0;
            0  -mu*xm3 mu*xm2 i22+i33     0      0;
         mu*xm3      0     0       0   i22   -i23;
        -mu*xm2      0     0       0  -i23    i33;
    ]
    # element triad
    Cab = @SMatrix [
        e1x e2x e3x;
        e1y e2y e3y;
        e1z e2z e3z;
    ]
    
    return GXBeam.Element(DeltaL, x, C, mass, Cab)
    end

function create_gxbeam_assembly(gxmodel::gxbeam, p, start, stop)
    ne = gxmodel.n
    np = gxmodel.n+1
    points = [create_gxbeam_point(view(p, 3*(ip-1) + 1 : 3*(ip-1) + 3)) for ip = 1:np] 
    elements = [create_gxbeam_element(view(p, 3*np + 36*(ie-1) + 1 : 3*np + 36*(ie-1) + 36),
    points, start[ie], stop[ie]) for ie = 1:ne]

    return GXBeam.Assembly(points, start, stop, elements)
end

function create_gxbeam_assembly(gxmodel::gxbeam, p) #For ease of use outside of the package. 
    ne = gxmodel.n
    np = gxmodel.n+1

    start = 1:gxmodel.n
    stop = 2:gxmodel.n+1

    points = [create_gxbeam_point(view(p, 3*(ip-1) + 1 : 3*(ip-1) + 3)) for ip = 1:np] 
    elements = [create_gxbeam_element(view(p, 3*np + 36*(ie-1) + 1 : 3*np + 36*(ie-1) + 36),
    points, start[ie], stop[ie]) for ie = 1:ne]

    return GXBeam.Assembly(points, start, stop, elements)
end



function get_element_velocity(states, ielem; allstates=false) #TODO: Apply multiple dispatch
    idx = 6 + 18*(ielem-1)
    if allstates
    end
    # @show idx
    return states[idx+13:idx+15]
end

function create_gxbeamfun(gxmodel::gxbeam, env::Environment, distributedload::Function; g = 9.817, damping=true, b=0.01)

    ### Todo. I might accidently be "globalizing" these variables. Allocating them outside of the function might actually be worse. -> I don't think that I am, because I tried to declare force_scaling a constant, and it balked saying that it is a local variable. 
    ### Create prescribed conditions
    prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)) # root section is fixed #0.000019 seconds 0.000003 seconds (4 allocations: 2.156 KiB)

    ### Create point masses
    #TODO: Do I want point masses to be a differentiable thing? If so then this should probably go inside gxbeamfun!(). 
    # create dictionary of point masses
    point_masses = Dict(1 => PointMass(0.0, SVector(0.0, 0.0, 0.0), @SMatrix zeros(3,3)))

    # force_scaling = 1.0 #"unsupported 'const' declaration on local variable". I guess this is a local variable.

    ### Create GXBeam pass ins.
    start = 1:gxmodel.n
    stop = 2:gxmodel.n+1
    N, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem = GXBeam.system_indices(start, stop, false; prescribed_points=1:1) #The third argument is whether or not the system is static. #prescribed_points=[1, gxmodel.n+1]
    # @show N
    # @show irow_point
    # @show irow_elem
    # @show icol_point
    # @show icol_elem
    # @show irow_elem1 
    # @show irow_elem2

    gxbeamfun! = function(outs, dx, x, p, t) #Todo: I'm not sure that this is type stable. Adding types to the inputs didn't decrease the number of allocations. 
        ### Create assembly
        assembly = create_gxbeam_assembly(gxmodel, p, start, stop) #0.000004 seconds (2 allocations: 7.172 KiB)

        ### Create distributed load 
        elements = view(assembly.elements, :) #0.000000 seconds 

        f = distributedload(t) #0.000007 seconds (7 allocations: 240 bytes) ... 0.000005 seconds (3 allocations: 80 bytes)
        m = @SVector zeros(3)
        # distributed_loads = Dict(ielem => DistributedLoads(assembly, ielem;fx = (s) -> f[1], fy = (s) -> f[2], fz = (s) -> f[3]) for ielem in 1:nelem)

        ### follower distributed loads
        f1_follower = f2_follower = @SVector zeros(3) 
        m1_follower = m2_follower = @SVector zeros(3)
        
        f1 = f2 = elements[1].L*f/2 #0.000002 seconds (3 allocations: 80 bytes)
        m1 = m2 = elements[1].L*m/2
        if damping
            ## Write a function that gets the element velocities. 
            ue = get_element_velocity(x, 1)

            ## Calculate the damping force
            Fd = -b.*ue

            ## Apply the damping force
            f1_follower += Fd./elements[1].L
            f2_follower += Fd./elements[1].L
        end
        distributed_loads = Dict(1 => GXBeam.DistributedLoads(f1, f2, m1, m2, f1_follower, f2_follower, m1_follower, m2_follower)) #0.000007 seconds (6 allocations: 3.922 KiB)
        # distributed_loads = Dict() #0.000002 seconds (4 allocations: 608 bytes)

        for ielem in 2:gxmodel.n #0.000007 seconds (50 allocations: 4.844 KiB)... 0.000004 seconds (36 allocations: 2.531 KiB)
            f1 = f2 = elements[ielem].L*f/2 #TODO: Why am I multiplying by the length of the element? Isn't the force integrated? 
            m1 = m2 = elements[ielem].L*m/2
            if damping
                ## Write a function that gets the element velocities. 
                ue = get_element_velocity(x, ielem)

                ## Calculate the damping force
                Fd = -b.*ue
                # @show ue

                ## Apply the damping force
                f1_follower += Fd./elements[ielem].L
                f2_follower += Fd./elements[ielem].L
                # @show f1_follower
                # @show typeof(f1_follower)
            end
            # if maximum(f1_follower)>0.0 #It does eventually go positive. 
            #     @show f1_follower
            # end
            distributed_loads[ielem] = GXBeam.DistributedLoads(f1, f2, m1, m2, f1_follower, f2_follower, m1_follower, m2_follower)
        end

        # println(distributed_loads[5])

        # for ielem in 1:gxmodel.n
        #     f1_follower = f2_follower = elements[ielem].L*f/2
        #     m1_follower = m2_follower = elements[ielem].L*m/2
        #     distributed_loads[ielem] = GXBeam.DistributedLoads(f1, f2, m1, m2, f1_follower, f2_follower, m1_follower, m2_follower)
        # end


        ### System velocities and accelerations
            ### create the origin
        x0 = @SVector zeros(3)
        v0 = @SVector zeros(3) #System linear velocity #0.000000 seconds
        omega0 = SVector(0.0, 0.0, env.Omega(t)) #System angular velocity #0.000000 seconds
        a0 = @SVector zeros(3) #System linear acceleration #0.000000 seconds
        alpha0 = SVector(0.0, 0.0, env.Omegadot(t)) #System angular acceleration #0.000000 seconds

        ### Create gravity vector #TODO: Convert this from steady to dependent on actual position? 
        alpha = pi/2 - env.Omega(t)*t #Assuming the blade starts horizontal. # 0.000004 seconds (5 allocations: 80 bytes) #Todo. Why are there allocations here? This is scalar math.... evaluating env.Omega is two allocations. 0.000000 seconds The Environment object wasn't type stable.  
        gx = -g*cos(alpha)
        gy = -g*sin(alpha)
        gz = 0.0
        gvec = SVector(gx, gy, gz) #[gx gy gz]

        force_scaling = GXBeam.default_force_scaling(assembly) #"unsupported 'const' declaration on local variable". I guess this is a local variable.

        # if any(isnan, dx)
        #     println("dx: ", dx)
        # end

        # if any(isnan, x)
        #     println("x: ", x)
        # end
        # println("")
        # println(outs[end-18:end])

        #Todo. Are the states actually the way that I think that they are? Like... Is GXBeam expecting a different set of states? Now they are. I needed to cut the excess states (the zero states). 

        # println(length(outs)) #This is the length that I expect it to be. 


        # println(typeof(outs)) #Vector{Float64} I'm not sure that this can be a static vector. 
        # println(typeof(dx)) #Depends on what is passed in. 
        # println(typeof(irow_point)) #Vector{Int64}

        # println(typeof(prescribed_conditions))
        # println(typeof(distributed_loads))
        # println(typeof(point_masses))

        # @show distributed_loads

        GXBeam.dynamic_system_residual!(outs, dx, x, assembly, prescribed_conditions, distributed_loads, point_masses, gvec, force_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem, x0, v0, omega0, a0, alpha0) # 0.000020 seconds (17 allocations: 5.266 KiB)
        
        return outs
    end
    return SciMLBase.DAEFunction{true, true}(gxbeamfun!)
end


function differentialvars(gxmodel; include_extra_states=false) #Todo. Make sure that these actually line up with the states that I think they do. It appears Taylor uses fewer states than what he said in AerostructuralDynamics. He is eliminating states of points that aren't needed. He keeps the states of any constrained point, and the last point. I'm going to assume that he'd also keep the states of the first point if it wasn't constrained.
    np = gxmodel.n + 1
    ne = gxmodel.n

    if include_extra_states
        differential_vars = fill(false, 6*np + 18*ne)

        for i in 1:ne
            idx = 6*np + 18*(i-1)
            differential_vars[idx+1:idx+3] .= true # u (for the beam element)
            differential_vars[idx+4:idx+6] .= true # θ (for the beam element)
            differential_vars[idx+13:idx+15] .= true # V (for the beam element)
            differential_vars[idx+16:idx+18] .= true # Ω (for the beam element)
        end
    else
        differential_vars = fill(false, 12 + 18*ne)

        for i in 1:ne
            idx = 6 + 18*(i-1) #Skip past the points states. There are 6 at the beginning and 6 at the end. 
            differential_vars[idx+1:idx+3] .= true # u (for the beam element)
            differential_vars[idx+4:idx+6] .= true # θ (for the beam element)
            differential_vars[idx+13:idx+15] .= true # V (for the beam element)
            differential_vars[idx+16:idx+18] .= true # Ω (for the beam element)
        end
    end
    return differential_vars
end

function rotate_x(alpha_x)
    return [
        1.0     0.0             0.0;
        0.0     cos(alpha_x)   -sin(alpha_x);
        0.0     sin(alpha_x)    cos(alpha_x)]
end

function rotate_y(alpha_y)
    return [cos(alpha_y) 0 sin(alpha_y);
            0.0 1.0 0.0;
            -sin(alpha_y) 0 cos(alpha_y)]
end

function rotate_z(alpha_z)
    return [cos(alpha_z) -sin(alpha_z) 0.0;
            sin(alpha_z) cos(alpha_z) 0.0;
            0.0 0.0 1.0]
end

function create_simplebeam(radii, chords, twists, rhub, rtip, thicknesses; density=2.69e3, E= 6.83e10, nu=0.34, ky = 1.2000001839588001, kz = 14.625127919304001, kt = 65.85255016982444) #Aluminum default, AZO materials
    # shear and torsion correction factors: ky, kz, kt
    # material properties: E, nu

    N = length(radii) #Number of elements

    ### Create the point location parameters
    interior_points = zeros(N-1)
    for i=1:N-1
        interior_points[i] = (radii[i] + radii[i+1])/2
    end

    points = vcat(rhub, interior_points, rtip) #All of the points together

    pp = zeros(3*(N+1))

    for i = 1:(N+1) #Beam extends in x direction -> points go into every third position of pp, starting with the first position. . 
        idx = 3*(i-1)
        pp[idx+1] = points[i] 
    end

    ### Create the element parameters
    pe = zeros(36*N)

    ## Create the orientation vectors of the local reference frame (Cab matrix)

    Cab = [1.0 0.0 0.0; #e1x e1y e1z
           0.0 1.0 0.0  #e2x e2y e2z
           0.0 0.0 1.0] #e3x e3y e3z

    for i = 1:N
        idx = 36*(i-1)

        # Local Cab matrix (The local reference frame should be rotated by the local twist. The beam cross section will initially be in the YZ plane with the horizontal in the Z. For convinience, I'm going to let that be and rotate the airfoil here. So the airfoil needs.... actually... that's mildly complicated. It depends on how the airfoil is defined. So.... let's deal with that problem when we get there. For now we'll just rotate the beam so the width is vertical. 
        # Cab_local = rotate_x(pi/2 - twists[i])*Cab #Todo: Wait... If my beam is extending in the x direction, and I have no rotation... then .... I shouldn't be multiplying by a shift. right? 
        # Cab_local = rotate_x(pi/2 - twists[i])'*Cab*rotate_x(pi/2 - twists[i]) #todo:

        Cab_local = Cab
        # if i==1
        #     @show rotate_x(pi/2 - twists[i])
        #     @show twists[i]
        #     @show Cab_local
        # end

        pe[idx+1:idx+3] = Cab_local[1,:] # Store e1
        pe[idx+4:idx+6] = Cab_local[2,:] # Store e2
        pe[idx+7:idx+9] = Cab_local[3,:] # Store e3

        ## Create compliance matrix elements 
        # cross section
        w = chords[i] # inch
        h = thicknesses[i] #0.063 # inch

        A = h*w
        Iz = w*(h^3)/12 #Second moment of area
        Iy = h*(w^3)/12
        Izy = 0.0

        # if i==1
        #     println("Inside CreateSimpleBeam")
        #     @show Iz, Iy
        # end

        Izz, Iyy, Izy = rotate_smoa(Iz, Iy, Izy, pi/2 - twists[i]) #Todo: Izy might be non-zero now.... I wonder if that'll play with the compliance matrix. 

        # Izz = Iz
        # Iyy = Iy
        # if i==1
        #     @show Izz, Iyy, Izy
        # end

        J = Iyy + Izz #Second polar moment of area about axis through centroid

        # apply corrections
        # Ay = A/ky
        # Ax = A/kz
        # Jx = J/kt

        G = E/(2*(1+nu)) ##Todo. Is this to correct relationship? -> For an isotropic material, yes. 

        #= 
            - Why does he have the area moment of inertia in the compliance matrix? Should I have that? -> Yes, he uses the integrated unit length compliance matrix, so yes, it has the area moment of inertia in the bottom.

            - What about the other terms. I did just the main diagonal here. That's fairly problematic. -> Added entries for an isotripic material. The isotropic compliance matrix looks like the orthotropic compliance matrix. But simpler. -> Yes, but this is the integrated unit length compliance matrix, so it won't look like the compliance matrix that I was previously looking at. It will in fact be just the main diagonal. 

            - Who determines what axis these mass moment of inertia should be. For instance, if I'm extending in the z direction, then should I have Iyy and Ixx? I'm going to guess so. -> NO. So, GXBeam assumes that your beam extends in the x direction. So all properties here should be Y and Z properties. 
        =#


        ### Trying Taylor's isotropic version
        pe[idx + 10] = 1/(E*A) #c11 
        pe[idx + 16] = 1/(G*A) #c22
        pe[idx + 21] = 1/(G*A) #c33
        pe[idx + 25] = (1)/(G*J) #c44 
        pe[idx + 28] = (1)/(E*Iyy) #c55 
        pe[idx + 30] = (1)/(E*Izz) #c66 

        

        xm2 = 0.0 #For an isotropic beam, the center of mass will be in the center of the beam. 
        xm3 = 0.0

        ## Create mass matrix elements
        pe[idx + 31] = density*A # mu (distributed weight) #Todo. This is going to be in kilograms, not Newtons. Should this be in Newtons? -> No, this should be kg/m, so mass per length. This will be used to calculated accelerations based on forces. 
        pe[idx + 32] = xm2 #xm2 distance to center of mass 
        pe[idx + 33] = xm3 #xm3 distance to center of mass
        pe[idx + 34] = density*Iyy #i22 # Mass moment of inertia #For an isotropic beam, the unit mass moment of inertia is equal to the density times the second moment of area. 
        pe[idx + 35] = density*Izz #i33
        pe[idx + 36] = density*Izy #i23 #Todo. This says i23. Should it be J, or should it be Izy? He says on the getting started page that it should be the product of inertia. 
    end

    ### Create p
    p = vcat(pp, pe)

    return N, p
end

function convert_assemblystate(state;include_extra_states=false) #Convert a the assembly state of a beam with a fixed end into a vector of states. 
    points_states = view(state.points, :)
    elements_states = view(state.elements, :)

    np = length(points_states)
    ne = length(elements_states)

    ### Extract states
    if include_extra_states #Extract all the states, including the ones that are likely not used. 
        xp = zeros(6*np) #Initialize array 

        #The first states are forces and moments since the first state is constrained not to displace. 
        xp[1:3] = points_states[1].F
        xp[4:6] = points_states[1].M

        #Iterate through all the other points, their states are displacement states. 
        for i = 2:np
            idx = 6*(i-1)
            xp[idx+1:idx+3] = points_states[i].u #Linear displacement states
            xp[idx+4:idx+6] = points_states[i].theta #Angular displacement states
        end

        xe = zeros(18*ne)

        #Iterate through all the elements, each element uses all 18 states. 
        for i = 1:ne
            idx = 18*(i-1)
            xe[idx+1:idx+3] = elements_states[i].u #Linear displacement states
            xe[idx+4:idx+6] = elements_states[i].theta #Angular displacement states
            xe[idx+7:idx+9] = elements_states[i].F #Force states
            xe[idx+10:idx+12] = elements_states[i].M #Moment states
            xe[idx+13:idx+15] = elements_states[i].V #Linear velocity states
            xe[idx+16:idx+18] = elements_states[i].Omega #Angular velocity States
        end
    else #Extract the used states. Since this is a cantilevered beam, GXBeam only uses the states of the first and last points. 
        ## Extract point states
        xp = zeros(12)
        # Extract states for the first point
        xp[1:3] = points_states[1].F
        xp[4:6] = points_states[1].M

        # Extract states for the last point
        xp[end-5:end-3] = points_states[end].u
        xp[end-2:end] = points_states[end].theta


        ## Extract Element states
        xe = zeros(18*ne)
        for i = 1:ne
            idx = 18*(i-1)
            xe[idx+1:idx+3] = elements_states[i].u
            xe[idx+4:idx+6] = elements_states[i].theta
            xe[idx+7:idx+9] = elements_states[i].F
            xe[idx+10:idx+12] = elements_states[i].M
            xe[idx+13:idx+15] = elements_states[i].V
            xe[idx+16:idx+18] = elements_states[i].Omega
        end
    end
    
    ## Combine the point and element states and return them.
    return vcat(xp, xe)
end

function convert_history(history; include_extra_states=false)
    np = length(history[1].points)
    ne = length(history[1].elements)
    nt = length(history)
    
    if include_extra_states
        n = 6*np + 18*ne
    else
        n = 12 + 18*ne
    end

    x = zeros(nt, n)
    for i = 1:nt
        x[i,:] = convert_assemblystate(history[i])
    end
    
    return x
end

function initialize_gxbeam(gxmodel, p, distributedload)

    start = 1:gxmodel.n
    stop = 2:gxmodel.n+1

    assembly = create_gxbeam_assembly(gxmodel, p, start, stop)
    elements = view(assembly.elements, :)

    ### Create prescribed conditions
    prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)) # root section is fixed

    ### Create distributed loads
    f1_follower = f2_follower = @SVector zeros(3) 
    m1_follower = m2_follower = @SVector zeros(3)
    
    m = @SVector zeros(3)

    f = distributedload(0.0)

    distributed_loads = Dict()
    for ielem in 1:gxmodel.n
        f1 = f2 = elements[ielem].L*f/2
        m1 = m2 = elements[ielem].L*m/2
        distributed_loads[ielem] = GXBeam.DistributedLoads(f1, f2, m1, m2, f1_follower, f2_follower, m1_follower, m2_follower)
    end

    system, converged = static_analysis(assembly; prescribed_conditions = prescribed_conditions, distributed_loads = distributed_loads, linear = false) #Todo. Should this be linear? Or nonlinear? Nonlinear, this is a nonlinear problem. 

    state = AssemblyState(system, assembly; prescribed_conditions = prescribed_conditions)

    return convert_assemblystate(state)
end

function initialize_gxbeam2(gxmodel, p, distributedload) #Todo: Need to add kwargs that go to GXBeam's initalize function so that I can have the same initial conditions. 

    start = 1:gxmodel.n
    stop = 2:gxmodel.n+1

    assembly = create_gxbeam_assembly(gxmodel, p, start, stop)
    elements = view(assembly.elements, :)

    ### Create prescribed conditions
    prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)) # root section is fixed

    ### Create distributed loads
    # f1_follower = f2_follower = @SVector zeros(3) 
    # m1_follower = m2_follower = @SVector zeros(3)
    
    # m = @SVector zeros(3)

    # f = distributedload(0.0)

    # distributed_loads = Dict()
    # for ielem in 1:gxmodel.n
    #     f1 = f2 = elements[ielem].L*f/2
    #     m1 = m2 = elements[ielem].L*m/2
    #     distributed_loads[ielem] = GXBeam.DistributedLoads(f1, f2, m1, m2, f1_follower, f2_follower, m1_follower, m2_follower)
    # end

    f1 = f2 = @SVector zeros(3) 
    m1 = m2 = @SVector zeros(3)
    
    m = @SVector zeros(3)

    f = distributedload(0.0)

    distributed_loads = Dict()
    for ielem in 1:gxmodel.n
        f1_follower = f2_follower = elements[ielem].L*f/2
        m1_follower = m2_follower = elements[ielem].L*m/2
        distributed_loads[ielem] = GXBeam.DistributedLoads(f1, f2, m1, m2, f1_follower, f2_follower, m1_follower, m2_follower)
    end 

    system, converged = GXBeam.initial_condition_analysis(assembly, 0.0; prescribed_conditions, distributed_loads)

    state = AssemblyState(system, assembly; prescribed_conditions = prescribed_conditions)

    return convert_assemblystate(state)
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

# function parsesolution(sol, gxmodel::gxbeam) #Todo: This isn't complete/working
#     t = sol.t
#     x = Array(sol)'
#     n = length(t)

#     # diffvars = differentialvars(gxmodel)

#     points = Array{Vector}(undef, n)
#     elements = Array{Vector}(undef, n)

#     pnts = Array{Vector}(undef, gxmodel.n+1)
#     els = Array{Vector}(undef, gxmodel.n)

    
#     for i = 1:n

#         xl = x[i,:]

#         #Extract point states
#         for j = 1:gxmodel.n+1
#             pidx = 6*(j-1)
#             pnts[j] = xl[pidx+1:pidx+6]
            
#         end

#         #Extract Element states
#         for j = 1:gxmodel.n
#             elidx = 18*(j-1) + 6*(gxmodel.n+1)
#             els[j] = xl[elidx+1:elidx+18]
#         end
#         points[i] = pnts
#         elements[i] = els
#     end
#     return t, x, points, elements
# end


function secondmomentofarea(x, y)
    nx = length(x)
    ny = length(y)

    if nx!=ny
        error("Second moment of area function: X and Y vectors must be the same length")
    end

    ### Todo: I should probably do some sort of check to see if the points make a shape

    ### Initialize the moments of area
    Ix = 0.0
    Iy = 0.0
    Ixy = 0.0

    ### Iterate through all the points but the last
    for i = 1:nx-1
        v = (x[i]*y[i+1] - x[i+1]*y[i])
        Ix += v*(y[i]^2 + y[i]*y[i+1] + y[i+1]^2)
        Iy += v*(x[i]^2 + x[i]*x[i+1] + x[i+1]^2)
        Ixy += v*(x[i]*y[i+1] + 2*x[i]*y[i] + 2*x[i+1]*y[i+1] + x[i+1]*y[i])
    end

    ### Calculate the last point
    Ix += (x[end]*y[1] - x[1]*y[end])*(y[end]^2 + y[end]*y[1] + y[1]^2)
    Iy += (x[end]*y[1] - x[1]*y[end])*(x[end]^2 + x[end]*x[1] + x[1]^2)
    Ixy += (x[end]*y[1] - x[1]*y[end])*(x[end]*y[1] + 2*x[end]*y[end] + 2*x[1]*y[1] + x[1]*y[end])

    ### Divide by the constant and return the value. 
    return Ix/12, Iy/12, Ixy/24
end

function secondmomentofarea(pts)
    return secondmomentofarea(pts[:,1], pts[:,2])
end

"""
    rotate_smoa(Ix, Iy, Ixy, phi)

Calculates the second moments of area (Area moment of inertia) about an axis that is phi radians counterclockwise from the original axes. 

### Inputs
- Ix - The second moment of area about the x axis
- Iy - The second moment of area about the Y axis
- Ixy - The product moment of area a bout the XY axis

### Outputs
- Iu - 
- Iv - 
- Iuv - 
"""
function rotate_smoa(Ix, Iy, Ixy, phi)
    p = (Ix + Iy)/2
    n = (Ix - Iy)/2

    Iu = p + n*cos(2*phi) - Ixy*sin(2*phi)
    Iv = p - n*cos(2*phi) + Ixy*sin(2*phi)
    Iuv = n*sin(2*phi) + Ixy*cos(2*phi)

    return Iu, Iv, Iuv
end

struct assemblystate{Tp, Te} #Stand in struct for GXBeam.AssemblyState because it doesn't have the constructor I want. 
    points::Tp
    elements::Te
end

function parsesolution(sol, gxmodel, p) #Todo: Broken. 

    start = 1:gxmodel.n
    stop = 2:gxmodel.n+1
    N, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem = GXBeam.system_indices(start, stop, false; prescribed_points=1:1)

    # @show icol_elem

    assembly = create_gxbeam_assembly(gxmodel, p, start, stop)

    ### Unpack
    states = Array(sol)'
    elements = assembly.elements

    ### Find the number of time steps, states, elements, and points
    nt, ns = size(states)
    ne = Int((ns-12)/18)
    np = ne + 1

    history = Vector{assemblystate}(undef, nt)
    
    for i = 1:nt
        pointstates = Vector{GXBeam.PointState}(undef, np)
        elementstates = Vector{GXBeam.ElementState}(undef, ne)

        ### Save the initial point and element states #Checked, I'm arriving at the correct states. At least the correct index. 
        idx1 = 6 
        ue1 = SVector(states[idx1+1:idx1+3]...)
        te1 = SVector(states[idx1+4:idx1+6]...)
        fe1 = SVector(states[idx1+7:idx1+9]...)
        me1 = SVector(states[idx1+10:idx1+12]...)
        ve1 = SVector(states[idx1+13:idx1+15]...)
        oe1 = SVector(states[idx1+16:idx1+18]...)
        elementstates[1] = GXBeam.ElementState(ue1, te1, fe1, me1, ve1, oe1)
        fp1 = SVector(states[1:3]...)
        mp1 = SVector(states[4:6]...)
        pointstates[1] = GXBeam.PointState(@SVector(zeros(3)), @SVector(zeros(3)), fp1, mp1)

        ### Save the internal point and element states
        for j = 2:ne
            idx = 6 + 18*(j-1)
            # if i==3
            #     @show idx
            # end
            ue = SVector(states[idx+1:idx+3]...)
            te = SVector(states[idx+4:idx+6]...)
            fe = SVector(states[idx+7:idx+9]...)
            me = SVector(states[idx+10:idx+12]...)
            ve = SVector(states[idx+13:idx+15]...)
            oe = SVector(states[idx+16:idx+18]...)
            elementstates[j] = GXBeam.ElementState(ue, te, fe, me, ve, oe)

            L1 = elements[j-1].L/2
            L2 = elements[j].L/2
            L = L1+L2
            
            u = (elementstates[j-1].u).*(L2/L) + (elementstates[j].u).*(L1/L)
            theta = (elementstates[j-1].theta).*(L2/L) + (elementstates[j].theta).*(L1/L)

            # if j==5
            #     @show ue
            #     @show te
            # end

            pointstates[j] = GXBeam.PointState(SVector(u...), SVector(theta...), @SVector(zeros(3)), @SVector(zeros(3)))
        end

        ### Save final point to points vector
        pointstates[end] = GXBeam.PointState(SVector(states[end-5:end-3]...), SVector(states[end-2:end]...), @SVector(zeros(3)), @SVector(zeros(3)))


        ### Save the AssemblyState for this time step
        history[i] = assemblystate(pointstates, elementstates)
    end
    return history
end

function initializegravityloads(gxmodel, env, p; g=9.817)

    ### Create GXBeam Assembly
    assembly = create_gxbeam_assembly(gxmodel, p)

    ### Create prescribed conditions
    prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0))

    omega = SVector(0.0, 0.0, -env.Omega(0.0))
    gvec = SVector(0.0, -g, 0.0)

    system, converged = GXBeam.initial_condition_analysis(assembly, 0.0; prescribed_conditions=prescribed_conditions, gravity = gvec, angular_velocity=omega)

    state = GXBeam.AssemblyState(system, assembly; prescribed_conditions = prescribed_conditions)

    return convert_assemblystate(state)
end