
#=
############## States ###################

    x = [x_p, ..., x_e, ...] 

where x_p is a vector of states for a given point, and x_e is a vector states for a given gxbeam element. There is a x_p for every point in our beam, and a x_e for every gxbeam element. 

x_p has two options, depending if we are constraining displacement or not. Regardless of if we are constraining displacement, there are six elements of x_p; which implies that the total number of elements in x from the point states are 6*n_p (where n_p is the number of points).

    If we're constraining displacement, then the vector of states is:
        x_p = [F_x, F_y, F_z, M_x, M_y, M_z]. 
    These are the forces and moments required to keep the point stationary. 

    If we're not constraining displacement (it has that degree of freedom), then the vector of states is x_p = [delta_x, delta_y, delta_z, theta_x, theta_y, theta_z], where delta and theta are linear and angular displacements respectively. 


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
function gxbeam_residual(residuals, dx, x, assembly, prescribed_conditions) #Not going to use this function, going to use the GXBeam function. 


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


function create_gxbeamfun(gxmodel::gxbeam, env::Environment, distributedload; g = 9.817)
    ### Create prescribed conditions
    prescribed_conditions = Dict(1 => GXBeam.PrescribedConditions(ux=0, uy=0, uz=0, theta_x=0, theta_y=0, theta_z=0)) # root section is fixed

    ### Create point masses
    #TODO: Do I want point masses to be a differentiable thing? If so then this should probably go inside gxbeamfun!(). 
    # create dictionary of point masses
    point_masses = Dict() #Dict(nelem => PointMass(m, p, J))

    force_scaling = 1.0

    ### Create GXBeam pass ins.
    start = 1:gxmodel.n
    stop = 2:gxmodel.n+1
    N, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem =
    GXBeam.system_indices(start, stop, false; prescribed_points=1:1) #The third argument is whether or not the system is static. 
    # @show N, irow_point, irow_elem, icol_point, icol_elem

    ### create the origin
    x0 = @SVector zeros(3)

    ### follower distributed loads
    f1_follower = f2_follower = @SVector zeros(3) 
    m1_follower = m2_follower = @SVector zeros(3)

    # f1 = f2 = @SVector zeros(3) 
    # m1 = m2 = @SVector zeros(3)

    m = @SVector zeros(3)

    function gxbeamfun!(outs, dx, x, p, t)
        ### Create assembly
        assembly = create_gxbeam_assembly(gxmodel, p, start, stop)

        ### Create distributed load 
        elements = view(assembly.elements, :)

        f = distributedload(t)

        # distributed_loads = Dict(ielem => DistributedLoads(assembly, ielem;fx = (s) -> f[1], fy = (s) -> f[2], fz = (s) -> f[3]) for ielem in 1:nelem)
        

        distributed_loads = Dict()
        for ielem in 1:gxmodel.n
            f1 = f2 = elements[ielem].L*f/2
            m1 = m2 = elements[ielem].L*m/2
            distributed_loads[ielem] = GXBeam.DistributedLoads(f1, f2, m1, m2, f1_follower, f2_follower, m1_follower, m2_follower)
        end

        # println(distributed_loads[5])

        # for ielem in 1:gxmodel.n
        #     f1_follower = f2_follower = elements[ielem].L*f/2
        #     m1_follower = m2_follower = elements[ielem].L*m/2
        #     distributed_loads[ielem] = GXBeam.DistributedLoads(f1, f2, m1, m2, f1_follower, f2_follower, m1_follower, m2_follower)
        # end


        ### System velocities and accelerations
        v0 = @SVector zeros(3) #System linear velocity
        omega0 = SVector(env.Omega(t), 0.0, 0.0) #System angular velocity
        a0 = @SVector zeros(3) #System linear acceleration
        alpha0 = SVector(env.Omegadot(t), 0.0, 0.0) #System angular acceleration

        ### Create gravity vector #TODO: Convert this from steady to dependent on actual position? 
        # g = 9.817 #Moving this to the bigger function header so that I can set it to zero if need be. 
        alpha = pi/2 - env.Omega(t)*t #Assuming the blade starts horizontal. 
        gx = -g*cos(alpha)
        gy = -g*sin(alpha)
        gz = 0.0
        gvec = SVector(gx, gy, gz) #[gx gy gz]

        # if any(isnan, dx)
        #     println("dx: ", dx)
        # end

        # if any(isnan, x)
        #     println("x: ", x)
        # end
        # println("")
        # println(outs[end-18:end])
        #Todo: Are the states actually the way that I think that they are? Like... Is GXBeam expecting a different set of states? 
        # println(length(outs)) #This is the length that I expect it to be. 

        GXBeam.dynamic_system_residual!(outs, dx, x, assembly, prescribed_conditions, distributed_loads, point_masses, gvec, force_scaling, irow_point, irow_elem, irow_elem1, irow_elem2, icol_point, icol_elem, x0, v0, omega0, a0, alpha0) #-> There is another function dynamic_element_residual. 
        # I think that this function is actually updating the outs, like I think they are changing, I see a couple instances of the residuals changing when I call it between prints. 

        # println(outs[end-18:end])

        # if any(isnan, outs)
        #     println("t: ", t)
        #     println("outs: ", outs)
        # end
    end
    return gxbeamfun!
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
            idx = 12 + 18*(i-1)
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
        Cab_local = rotate_x(pi/2 - twists[i])*Cab #Todo: Wait... If my beam is extending in the x direction, and I have no rotation... then .... I shouldn't be multiplying by a shift. right? 

        pe[idx+1:idx+3] = Cab_local[1,:] # Store e1
        pe[idx+4:idx+6] = Cab_local[2,:] # Store e2
        pe[idx+7:idx+9] = Cab_local[3,:] # Store e3

        ## Create compliance matrix elements 
        # cross section
        w = chords[i] # inch
        h = thicknesses[i] #0.063 # inch

        A = h*w
        Izz = w*(h^3)/12 #Todo. Who determines what axis these mass moment of inertia should be. For instance, if I'm extending in the z direction, then should I have Iyy and Ixx? I'm going to guess so. -> NO. So, GXBeam assumes that your beam extends in the x direction. So all properties here should be Y and Z properties. 
        Iyy = (w^3)*h/12
        J = Iyy + Izz #Second polar moment of area about axis through centroid

        # apply corrections
        # Ay = A/ky
        # Ax = A/kz
        # Jx = J/kt

        G = E/(2*(1+nu)) ##Todo. Is this to correct relationship? -> For an isotropic material, yes. 

        # pe[idx + 10] = 1/(E*A) #c11 #From one of Taylor's examples. 
        # pe[idx + 16] = 1/(G*Ay) #c22
        # pe[idx + 21] = 1/(G*Az) #c33
        # pe[idx + 25] = 1/(G*Jx) #c44
        # pe[idx + 28] = 1/(E*Iyy) #c55 #TODO. Why does he have the area moment of inertia in the compliance matrix? Should I have that? -> Yes, he uses the integrated unit length compliance matrix, so yes, it has the area moment of inertia in the bottom. 
        # pe[idx + 30] = 1/(E*Izz) #c66

        #Isotropic compliance matrix #Todo. Should I be passing in the stiffness or compliance matrix? -> It seems obvious that it should be the compliance matrix, since that's what is stored in the element and that's what Taylor's commenting says. 
        # pe[idx + 10] = 1/(E*A) #c11 #Todo. It appears that Taylor isn't using the normal compliance matrix. It might be a normalized version... See documentation, getting started. Also, see comments above. 
        # pe[idx + 11] = -nu/(E*A) #c12
        # pe[idx + 12] = -nu/(E*A) #c13
        # pe[idx + 16] = 1/(E*A) #c22
        # pe[idx + 17] = -nu/(E*A) #c23
        # pe[idx + 21] = 1/(E*A) #c33
        # pe[idx + 25] = (1+nu)/(E*J) #c44 *2
        # pe[idx + 28] = (1+nu)/(E*Iyy) #c55 *2
        # pe[idx + 30] = (1+nu)/(E*Izz) #c66 *2 -? acutal strain vs engineering strain

        ### Trying Taylor's isotropic version
        pe[idx + 10] = 1/(E*A) #c11 
        pe[idx + 16] = 1/(G*A) #c22
        pe[idx + 21] = 1/(G*A) #c33
        pe[idx + 25] = (1)/(G*J) #c44 *2
        pe[idx + 28] = (1)/(E*Iyy) #c55 *2
        pe[idx + 30] = (1)/(E*Izz) #c66 *2 -? acutal strain vs engineering strain

        #Todo. What about the other terms. I did just the main diagonal here. That's fairly problematic. -> Added entries for an isotripic material. 
        # The isotropic compliance matrix looks like the orthotropic compliance matrix. But simpler. -> Yes, but this is the integrated unit length compliance matrix, so it won't look like the compliance matrix that I was previously looking at. It will in fact be just the main diagonal. 

        xm2 = 0.0 #For an isotropic beam, the center of mass will be in the center of the beam. 
        xm3 = 0.0

        ## Create mass matrix elements
        pe[idx + 31] = density*A # mu (distributed weight) #Todo. This is going to be in kilograms, not Newtons. Should this be in Newtons? -> No, this should be kg/m, so mass per length. This will be used to calculated accelerations based on forces. 
        pe[idx + 32] = xm2 #xm2 distance to center of mass 
        pe[idx + 33] = xm3 #xm3 distance to center of mass
        pe[idx + 34] = density*Iyy #i22 # Mass moment of inertia #For an isotropic beam, the unit mass moment of inertia is equal to the density times the second moment of area. 
        pe[idx + 35] = density*Izz #i33
        pe[idx + 36] = density*J #i23
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

function initialize_gxbeam2(gxmodel, p, distributedload)

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

function parsesolution(sol, gxmodel::gxbeam)
    t = sol.t
    x = Array(sol)'
    n = length(t)

    # diffvars = differentialvars(gxmodel)

    points = Array{Vector}(undef, n)
    elements = Array{Vector}(undef, n)

    pnts = Array{Vector}(undef, gxmodel.n+1)
    els = Array{Vector}(undef, gxmodel.n)

    
    for i = 1:n

        xl = x[i,:]

        #Extract point states
        for j = 1:gxmodel.n+1
            pidx = 6*(j-1)
            pnts[j] = xl[pidx+1:pidx+6]
            
        end

        #Extract Element states
        for j = 1:gxmodel.n
            elidx = 18*(j-1) + 6*(gxmodel.n+1)
            els[j] = xl[elidx+1:elidx+18]
        end
        points[i] = pnts
        elements[i] = els
    end
    return t, x, points, elements
end