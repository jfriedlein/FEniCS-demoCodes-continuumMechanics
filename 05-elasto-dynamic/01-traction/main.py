#! /usr/bin/env python

from __future__ import print_function
from fenics import *

# Enable C++ optimizations in form compilation
parameters["form_compiler"]["cpp_optimize"] = True 
parameters["form_compiler"]["optimize"] = True


"""-------------------------------------------------2D Elastodynamics of a Rectangular Solid--------------------------------------------------
Problem description: time ddependent prescribed traction

Geometry: Rectangular solid with dimensions 0.4 m x 0.1 m
Boundary conditions:
   -Left face is fully clamped
Loads:
   - A gradually increasing traction is applied on the right face at an angle for 10/3 ms and then removed 
   - No body forces are considered for 2D

Analysis type: Quasi-static model
Material model: Linear, Isothermal, Isotropic elasticity


Main Learnings:

Newmark method
Newmark predictors
Newmark approximations for rates
variational formulation - weak form (newmark)
Performing time integration and function solving        
--------------------------------------------------------------------------------------------------------------------------------------------------

"""
# Time definitions
t = 0.0 # initial time
T = 0.01 # final time in seconds (10 ms)
num_steps = 300 # number of time steps
dt = T / num_steps # time step size

#Parameters for Newmark (Average acceleration)
nm_beta = 0.25 # Newmark beta parameter
nm_gamma = 0.5 # Newmark gamma parameter

#---------------------------------------------------------------------------------------------------------
# Geometry and mesh generation
#---------------------------------------------------------------------------------------------------------

# spatial dimension
dim=2

# Define geometry and mesh
p0 = Point(0.0, 0.0) # bottom-left front corner
p1 = Point(0.4, 0.1) # top-right back corner

#p0 = Point(0.0, 0.0, 0.0) # bottom-left front corner  3D
#p1 = Point(0.4, 0.1, 0.1) # top-right back corner 3D

mesh = RectangleMesh(p0, p1, 20, 5) # if "Domain" is 2D, use this to mesh
#mesh = BoxMesh(p0, p1, 4, 1, 1) # if "Domain" is 3D, use this to mesh


#---------------------------------------------------------------------------------------------------------
# Boundary identification and marking
#---------------------------------------------------------------------------------------------------------


# MeshFunction to store boundary IDs on facets
boundaries = MeshFunction("size_t", mesh, dim-1)  # create a MeshFunction to store boundary IDs on facets (edges in 2D, faces in 3D)
boundaries.set_all(0) # set every face to the value 0
left, right, bottom, top = 1, 2, 3, 4 # define boundary IDs
CompiledSubDomain("near(x[0], side) && on_boundary", side = p0[0]).mark(boundaries, left)  # mark left face with ID 1 
CompiledSubDomain("near(x[0], side) && on_boundary", side = p1[0]).mark(boundaries, right) #  mark right face with ID 2
CompiledSubDomain("near(x[1], side) && on_boundary", side = p0[1]).mark(boundaries, bottom) # mark bottom face with ID 3
CompiledSubDomain("near(x[1], side) && on_boundary", side = p1[1]).mark(boundaries, top) # mark top face with ID 4

# Coordinates and surface integral element
x = SpatialCoordinate(mesh) # spatial coordinates

#  Surface integral element
ds = Measure('ds', domain=mesh, subdomain_data=boundaries) 


#-------------------------------------------------------------------------------------------------------
# Function spaces and variational functions
#-------------------------------------------------------------------------------------------------------

# spatial dimension
dim = 2

# Polynomial degree of Lagrange finite elements
p = 2

# Define vector function space for displacement field.
V = VectorFunctionSpace(mesh, "Lagrange", p)


# unknown displacement field to be solved for
u = TrialFunction(V)

# virtual displacement field for the variational formulation
delta_u = TestFunction(V)

# Interpolate initial displacement, velocity
u_n = Function(V) # previous displacement
u_n.interpolate(Constant((0.0, 0.0))) # assign initial displacement to zero
vel_n = Function(V) # previous velocity
vel_n.interpolate(Constant((0.0, 0.0))) # assign initial velocity to zero

# Define additional rates and predictor functions for Newmark
acc_n = Function(V) # previous acceleration
u_pred = Function(V) # predictor displacement
vel_pred = Function(V) # predictor velocity


#--------------------------------------------------------------------------------------------------------
# Material properties
#--------------------------------------------------------------------------------------------------------

# Material parameters
E = 200.e9 # Young's modulus in Pascals
nu = 0.3 # Poisson's ratio 
#mu    = E/(2.0*(1.0 + nu)) 
#lmbda = E*nu/((1.0 + nu)*(1.0 - 2.0*nu))
mu    = E/(2.0*(1.0 + nu)) #PLANE STRESS # shear modulus # lame second parameter
lmbda = E*nu/((1.0 + nu)*(1.0 - nu)) #PLANE STRESS  # lame first parameter

rho = 8.e3 # density in kg/m^3
g = 9.81 # gravitational acceleration in m/s^2

# Volume force 
b = Constant((0.0, 0.0)) 

#b = as_vector((0.0, -rho*g)) #if volume forces are considered

# Expression for prescribed tractions

# Expression is written such that traction increases linearly to a maximum value at t=t1 and then becomes zero for t>t1
t_p = Expression(('(t<1.0*t1)?1.0*(-m/t1*fabs(t-t1)+m):0.0', '(t<1.0*t1)?-0.1*(-m/t1*fabs(t-t1)+m):0.0'), degree=1, t1=0.3333*T, m=1.e6, t=0) # time-dependent traction load

#t_p = Expression(('(t<1.0*t1)?0.0*(-m/t1*fabs(t-t1)+m):0.0', '(t<1.0*t1)?-0.9*(-m/t1*fabs(t-t1)+m):0.0'), degree=1, t1=0.3333*T, m=1.e6, t=0) # alternative load case


# Prescribed Dirichlet boundary data
#u_p = Expression(('m*t', '0.0', '0.0'), degree=1, m=-0.2/T, t=0)
#theta_p = Expression(('m*t'), degree=1, m=100.0/T, t=0)


# Dirichlet boundary conditions
bcs = [DirichletBC(V, Constant((0.0, 0.0)), boundaries, left)] # fully fixed left face


#---------------------------------------------------------------------------------------------------------
#  Variational formulation (weak form) (newmark)
#---------------------------------------------------------------------------------------------------------

# Strain tensor 
def epsilon(u):
      return sym(grad(u))

# Stress tensor
def sigma(u):
    return lmbda*tr(epsilon(u))*Identity(dim) + 2.0*mu*epsilon(u)

# Newmark approximations for rates

# acceleration approximation
def acc(u, u_pred):
	return (u-u_pred)/(dt*dt*nm_beta)

# velocity approximation
def vel(u, u_pred, v_pred):
	return v_pred+dt*nm_gamma*a(u, u_pred)

# Variational form (F(u, delta_u) = a(u, delta_u) - l(delta_u))
F = rho*dot(delta_u, acc(u, u_pred))*dx + inner(grad(delta_u), sigma(u))*dx - dot(b, delta_u)*dx - dot(t_p, delta_u)*ds(right)

# Project initial stress field
Z = TensorFunctionSpace(mesh, "Lagrange", p) # function space for stress field
stress_n = Function(Z)
stress_n = project(sigma(u_n)/1.e6, Z, solver_type="mumps") # compute and project initial stress field in MPa
	
# Create output files
u_n.rename("u","displacement")
xdmf_u = XDMFFile("displacement_in_meters.xdmf")
xdmf_u.parameters["flush_output"] = True # ensures that data is written to file immediately after each write call
xdmf_u.parameters["functions_share_mesh"] = True # allows multiple functions to share the same mesh in the output file, reducing file size and improving performance when writing multiple functions defined on the same mesh

stress_n.rename("sigma","stress")
xdmf_stress = XDMFFile("stress_in_MPa.xdmf")
xdmf_stress.parameters["flush_output"] = True
xdmf_stress.parameters["functions_share_mesh"] = True


#---------------------------------------------------------------------------------------------------------
# Solve the linear system
#---------------------------------------------------------------------------------------------------------

# assemble bilinear form of variational formulation
A = assemble(lhs(F))
[bc.apply(A) for bc in bcs] # apply boundary conditions to A

# Set solver parameters
print("Setup solver...")
solver = LUSolver(A, "mumps") # using linear solver mumps
solver.parameters["symmetric"] = True #it tells solver that matrix is symmetric
#solver.parameters["reuse_factorization"] = True # reuses the LU factorization not needing to factorize everytime

# displacement solution space
u = Function(V) 

# Time integration and function solving
""" 

1. Increment simulation time t += dt, update traction t_p.t = t for the current timestep
2. Predict displacement and velocity at the current timestep using Newmark formulas:
        Use previous timestep values (u_n, vel_n, acc_n) to compute u_pred and vel_pred
3. Assemble and solve the linear system
        Solve A u = b to obtain the corrected displacement u at this timestep
4. Update acceleration and velocity using the corrected displacement (Newmark method):
5. Update displacement for the next timestep: u_n = u
6. Compute and project stresses onto the tensor function space in MPa:
7. Write results to output files (.xdmf) for displacement and stress.
Repeat for all timesteps (num_steps)

"""
for n in range(num_steps):	
    t += dt # updates t every time step by dt
    print("Step ", n, " (t=", t, ")", sep="")
    
    # Update loads, boundary data, ...
    #u_p.t = t
    #theta_p.t = t
    t_p.t = t  # update the time-dependent traction to the current simulation time t

    # Define predictors
    # computes the current displacemnt and velocity using the previous time step displacemeent,velocity, acceleration using newmark method
    u_pred.vector()[:] = u_n.vector()+dt*vel_n.vector()+ 0.5*dt*dt*(1-2*nm_beta)*acc_n.vector() # Newmark predictor for displacement
    vel_pred.vector()[:] = vel_n.vector()+dt*(1-nm_gamma)*acc_n.vector() # Newmark predictor for velocity
    
    # Assemble and solve
    b = assemble(rhs(F)) # turn the expression to matrix values
    [bc.apply(b) for bc in bcs] # apply boundary conditions to b
    solver.solve(u.vector(), b) # solve u for b  #corrected diplacement u
        
    # Update and correct functions
    #(u, theta) = y.split(deepcopy=True)

    # Compute and update current step acceleration using corrected displacement (Newmark formula)
    acc_n.vector()[:] = (u.vector()-u_pred.vector())/(nm_beta*dt*dt)

    # Compute and update current step velocity using updated acceleration (Newmark formula)
    vel_n.vector()[:] = vel_pred.vector()+nm_gamma*dt*acc_n.vector()

    # update displacement values
    u_n.vector()[:] = u.vector()
    
    # Project stresses at 
    stress_n.assign(project(sigma(u_n)/1.e6, Z, solver_type="mumps")) #compute and project stress
    
    # Write step to files
    xdmf_u.write(u_n, t) # store displacement at current time step
    xdmf_stress.write(stress_n, t) #store stress at current time step
    
