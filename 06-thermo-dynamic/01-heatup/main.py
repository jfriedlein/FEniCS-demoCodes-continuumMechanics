#! /usr/bin/env python

from __future__ import print_function
from fenics import *
import numpy as np

"""------------------------------------------2D Transient Heat Conduction Crank-Nicolson Diffusion-------------------------------------------------------------------------------------------------

Problem Description: Linear Thermal Response of a Square Plate
Geometry: A 2D square domain (0.1 m x 0.1 m) discretized with a structured mesh of 2,500 quadrilateral elements (50x50).

Boundary Conditions:
Left Boundary: Fixed temperature (Dirichlet BC) maintained at 25°C.
Right Boundary: Non-uniform heat flux (Neumann BC).
Top/Bottom: Naturally insulated (Adiabatic) as no conditions are specified.

Initial Condition: The entire domain starts at a uniform temperature of 25°C.
Analysis Type: Transient Thermal Analysis using the Crank-Nicolson time-integration scheme.
Material Model: Linear isotropic thermal

Main Learnings
Transient Thermal Analysis - Captures the temporal evolution of temperature within a 2D domain.
Crank-Nicolson Time Integration -Implements an implicit predictor-corrector scheme for stable and accurate time stepping.

Possibilities for extension:
Can incorporate time-dependent heat sources or fluxes which we will see in the next code.
Adapt to more complex geometries or 3D analysis

-----------------------------------------------------------------------------------------------------------------------------------------------------------------"""

# Time definitions
t = 0.0 #intial time in seconds
T = 3600 #final time in seconds
num_steps = 200 # number of time steps
dt = T / num_steps # time step size

# Parameters for Crank-Nicolson time integration scheme
cn_alpha = 0.5

#---------------------------------------------------------------------------------------------------------
# Geometry and mesh generation
#---------------------------------------------------------------------------------------------------------

# Define geometry and mesh
dim = 2
p0 = Point(0.0, 0.0, 0.0) # lower left corner 
p1 = Point(0.1, 0.1, 0.1) # upper right corner 
mesh = RectangleMesh(p0, p1, 50, 50) # apply a structured mesh to the rectangle domain 

#mesh = BoxMesh(p0, p1, 20, 10, 10) # if 3D geometry is desired, use a box mesh instead of a rectangle mesh

#---------------------------------------------------------------------------------------------------------
# Boundary identification and marking
#---------------------------------------------------------------------------------------------------------

# MeshFunction to store boundary IDs on facets
boundaries = MeshFunction("size_t", mesh, dim-1) #function to store boundary IDs on facets
boundaries.set_all(0) # initialize all boundaries to 0
left, right, bottom, top = 1, 2, 3, 4 # boundary IDs
CompiledSubDomain("near(x[0], side) && on_boundary", side = p0[0]).mark(boundaries, left) # mark left boundary with ID 1
CompiledSubDomain("near(x[0], side) && on_boundary", side = p1[0]).mark(boundaries, right) # mark right boundary with ID 2
CompiledSubDomain("near(x[1], side) && on_boundary", side = p0[1]).mark(boundaries, bottom) # mark bottom boundary with ID 3
CompiledSubDomain("near(x[1], side) && on_boundary", side = p1[1]).mark(boundaries, top) # mark top boundary with ID 4

# Coordinates and surface integral element
x = SpatialCoordinate(mesh) # spatial coordinates for defining expressions

# Surface integral element 
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

#-------------------------------------------------------------------------------------------------------
# Function spaces and variational functions
#-------------------------------------------------------------------------------------------------------

p = 2  # polynomial degree

V = FunctionSpace(mesh, "Lagrange", p) # vector function space for displacements


# Define trial and test functions
theta = TrialFunction(V) # trial
delta_theta = TestFunction(V) # test

# Interpolate initial temperature
theta_n = Function(V) # function to store the solution at the previous time step
theta_n.interpolate(Constant(25.0)) # set initial temperature to 25 degrees Celsius

# if a non-uniform initial temperature distribution is desired, use an Expression instead of a Constant to interpolate the initial temperature.
# For example, the following code defines a Gaussian distribution of initial temperature centered at (0.0, 0.5) with a peak value of 100 degrees Celsius 
# and a standard deviation of 0.01.

# theta_n.interpolate(Expression("1.0e2*exp(-((x[0]-rx)*(x[0]-rx)+(x[1]-ry)*(x[1]-ry))/0.001)", rx=0.0, ry=0.5, degree=4))

# Define additional rates and predictor functions for the Crank-Nicolson time integration scheme.
dtheta_n = Function(V) # function to store the rate of change of temperature at the previous time step
theta_pred = Function(V) # function to store the predicted temperature at the current time step based on previous time step values and Crank-Nicolson scheme

#--------------------------------------------------------------------------------------------------------
# Material properties
#--------------------------------------------------------------------------------------------------------


rho = 8000.0 # density in kg/m^3
kappa = 20.0 # thermal conductivity in W/(m-K)
cv = rho*500.0 # specific heat capacity in J/(kg-K)

# If a spatially varying material property distribution is desired, 
# use conditional expressions to define the material properties as functions of spatial coordinates. 
# For example, the following code defines a two-phase material distribution where the left half of the domain has one set of properties 
# and the right half has another set of properties.

#rho = 1.0*conditional(x[1] < 0.5*p1[1], 7800.0, 7000.0) # density distribution where the left half of the domain has a density of 7800 kg/m^3 and the right half has a density of 7000 kg/m^3
#kappa = 1.0*conditional(x[1] < 0.5*p1[1], 81.0, 121.0) # thermal conductivity distribution where the left half of the domain has a thermal conductivity of 81 W/(m-K) and the right half has a thermal conductivity of 121 W/(m-K)
#cv = 1.0*rho*conditional(x[1] < 0.5*p1[1], 452.0, 387.0) # specific heat capacity distribution where the left half of the domain has a specific heat capacity of 452 J/(kg-K) and the right half has a specific heat capacity of 387 J/(kg-K)

# Heat source and prescribed heat fluxes
r = Constant(0.0) # volumetric heat source in W/(m^3) (no internal heat generation)
q_p = Expression("x[1]*1.e5", degree=4) # prescribed heat flux in W/(m^2) applied on the right boundary (ID 2) with a linear distribution that increases with y-coordinate
# degree=4 is used for the Expression to ensure sufficient accuracy in representing the spatial variation of the heat flux on the boundary. 
# Adjust the degree as needed based on the desired accuracy and computational cost.


# If a time-varying heat source or heat flux is desired, use an Expression that depends on time.

#r = Expression("1.0e6*sin(m*t)", m=2*pi/T, t=0, degree=1) 
#r = Expression("1.0e7*exp(-((x[0]-x0)*(x[0]-x0)+(x[1]-y0)*(x[1]-y0))/0.001)*m*t", x0=0.25, y0=0.25, m=1/T, t=0.0, degree=4)
#q_p = Expression("1.0e3*m*t", m=1/T, t=0, degree=1)
# If a spatially varying heat source or heat flux distribution is desired, use conditional expressions to define the heat source or heat flux as functions of spatial coordinates.
#q_p_bottom = Expression("x[0]<0.5 ? -1.e3*m*t : 0.0", m=1/T, t=0.0, degree=4)

# Dirichlet boundary conditions
bcs = [DirichletBC(V, Constant(25.0), boundaries, left)] # set temperature to 25 degrees Celsius on the left boundary

#---------------------------------------------------------------------------------------------------------
#  Variational formulation of the Heat Equation (weak form)
#---------------------------------------------------------------------------------------------------------

# Crank-N. approximations for rates

# compute the change in temperature (dtheta) based on the previous temperature (theta) and the predicted temperature (theta_pred) using the Crank-Nicolson scheme.
def dtheta(th, th_pred): 
	return (th-th_pred)/(dt*cn_alpha) 

# Form definition (F(u, delta_u) = a(u, delta_u) - l(delta_u))
# F= internal energy rate + Fourier's heat conduction + volumetric heat source - prescribed heat flux on the right boundary 
F = cv*dtheta(theta, theta_pred)*delta_theta*dx \
	+ kappa*dot(grad(delta_theta), grad(theta))*dx \
	- r*delta_theta*dx \
	- q_p*delta_theta*ds(right) # - q_p_bottom*delta_theta*ds(bottom) #https://en.wikiversity.org/wiki/Nonlinear_finite_elements/Weak_form_of_heat_equation

# Create output files
theta_n.rename("theta","temperature")
file_theta = File("temperature.pvd", "compressed")
xdmf_theta = XDMFFile("temperature_in_degrees.xdmf")
xdmf_theta.parameters["flush_output"] = True # ensures that data is written to file immediately after each write call
xdmf_theta.parameters["functions_share_mesh"] = True # allows multiple functions to share the same mesh in the output file, reducing file size and improving performance when writing multiple functions defined on the same mesh
xdmf_theta.write(theta_n, t)


#---------------------------------------------------------------------------------------------------------
# Time integration and function solving
#---------------------------------------------------------------------------------------------------------



# assemble bilinear form of variational formulation
A = assemble(lhs(F))
[bc.apply(A) for bc in bcs] # apply boundary conditions

# Set solver parameters
solver = LUSolver(A, "mumps") # Use MUMPS direct solver for efficiency and robustness in solving the linear systems at each time step
solver.parameters["symmetric"] = True # The stiffness matrix A is symmetric due to the nature of the thermal diffusion operator
#solver.parameters["reuse_factorization"] = True # Reuse the symbolic factorization of the matrix structure across time steps, which can significantly speed up the solution process since the sparsity pattern of A does not change over time.

theta = Function(V) # solution function to store the temperature at the current time step

# Time integration and function solving
""" 
1. Increment simulation time t += dt
2. Predict temperature at the current timestep using crank-nicolson formulas:
        Use previous timestep values (theta_n , dtheta_n ) to compute theta_pred
3. Assemble and solve the linear system
        Solve A theta = b to obtain the corrected temperature theta at this timestep
4.. Write results to output files (.xdmf) for temperature.
Repeat for all timesteps (num_steps)

"""
for n in range(num_steps):	
    t += dt # update simulation time
    print("Step ", n, " (t=", t, ")", sep="")

    # Update loads, boundary data, ...
    #u_p.t = t # update bc 
    #theta_p.t = tc # update bc for the prescribed temperature
    #t_p.t = t # update bc for the prescribed traction

#    Update the heat source and prescribed heat flux for the current time step if they are time-dependent.
#    r.rx = (p1[0]-p0[0])*t/(0.01*T) 
#    r.ry = (p1[1]-p0[1])*t/(0.01*T)

    #r.rx = 0.3*np.cos(2*pi*t/(0.3*T))+0.5
    #r.ry = 0.3*np.sin(2*pi*t/(0.4*T))+0.5
#    if t<0.5*T:
#        q_p.s = t*2/T
#    else:
#        q_p.m = 1.0
    #q_p.t = t
    #q_p_bottom.t = t
    
    # Define predictors for the current time step using crank-nicolson formulas
    theta_pred.vector()[:] = theta_n.vector()+dt*(1-cn_alpha)*dtheta_n.vector()
    
    # Assemble and solve
    b = assemble(rhs(F))
    [bc.apply(b) for bc in bcs] # apply boundary conditions to b
    solver.solve(theta.vector(), b)  # solve the linear system to find the corrected temperature theta at the current time step, using the assembled matrix A and the load vector b. 
    
    # Update and correct functions
    # Update the rate of change of temperature (dtheta_n) based on the difference between the previous temperature (theta) and the predicted temperature (theta_pred) using the Crank-Nicolson scheme.
    dtheta_n.vector()[:] = (theta.vector()-theta_pred.vector())/(cn_alpha*dt)    
    # Update the temperature at the previous time step (theta_n) to the current temperature (theta) for use in the next iteration of the time-stepping loop.
    theta_n.vector()[:] = theta.vector()    
    
    # Write step to files
    xdmf_theta.write(theta_n, t)

#-----------------------------------------------------------------------------------