#! /usr/bin/env python

from __future__ import print_function
from fenics import *
import numpy as np

"""--------------------------------------2D Transient Heat Conduction with Moving Gaussian Heat Source--------------------------------------

Problem Description:
This script performs a transient thermal analysis of a 2D square domain using the Crank-Nicolson time-integration scheme. 
The simulation models heat diffusion within a solid material subjected to a time-dependent moving Gaussian heat source.

Geometry:
- 2D square domain of 0.5 m y 0.5 m
- Structured mesh with 50 x 50 elements

Material Model:
- Linear isotropic elasticity:
    - Density (rho): 4500 kg/m³
- Thermal properties:
    - Thermal conductivity (k): 6.6 W/(m·K)
    - Specific heat capacity (cᵥ): 2.9 W/((m^3)·K)

Boundary Conditions:
- Dirichlet BC: Temperature fixed at 0°C on all boundaries
- Neumann BC: Prescribed heat flux (here set to 0) at the top boundary
- Adiabatic boundaries: naturally insulated (handled implicitly)

Analysis Features:
- Crank-Nicolson scheme for second-order temporal accuracy
- Time-dependent moving Gaussian heat source simulating a localized heat load moving across the domain in a square pattern
- Prediction-correction approach for time integration

Main Learnings:
- Transient thermal response of a material to a moving heat source
- Handling time-dependent heat sources and boundary conditions

Possibilities for extension:
Adapt to more complex geometries or 3D analysis

----------------------------------------------------------------------------------------------------------------------------------------------"""
# Time definitions
t = 0.0 #intial time in seconds
T = 600  #final time in seconds
num_steps = 200 # number of time steps
dt = T / num_steps # time step size

# Parameters Crank-Nicolson
cn_alpha = 0.5
#---------------------------------------------------------------------------------------------------------
# Geometry and mesh generation
#---------------------------------------------------------------------------------------------------------

# Define geometry and mesh
dim = 2
p0 = Point(0.0, 0.0, 0.0) # lower left corner 
p1 = Point(0.5, 0.5, 0.1) # upper right corner
#mesh = BoxMesh(p0, p1, 20, 10, 10) # if 3D geometry is desired, use a box mesh instead of a rectangle mesh
mesh = RectangleMesh(p0, p1, 50, 50) # apply a structured mesh to the rectangle domain

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


p = 2 # polynomial degree
V = FunctionSpace(mesh, "Lagrange", p) # Define function space for temperature field using Lagrange finite elements of degree p on the mesh

# Define trial and test functions
theta = TrialFunction(V) #trial
delta_theta = TestFunction(V) #test

# Interpolate initial temperature
theta_n = Function(V) # function to store the temperature at the previous time step
theta_n.interpolate(Constant(0.0)) # set initial temperature to 0.0

# If desired, an initial temperature distribution can be defined using an expression. For example, a Gaussian distribution centered at (0.0, 0.5) can be used as the initial condition for the temperature field:
#theta_n.interpolate(Expression("1.0e2*exp(-((x[0]-rx)*(x[0]-rx)+(x[1]-ry)*(x[1]-ry))/0.001)", rx=0.0, ry=0.5, degree=4))

# Define additional rates and predictor functions for Newmark
dtheta_n = Function(V) # function to store the rate of change of temperature at the previous time step
theta_pred = Function(V) # function to store the predicted temperature at the current time step based on previous time step values and Newmark scheme

#--------------------------------------------------------------------------------------------------------
# Material properties
#--------------------------------------------------------------------------------------------------------

# define conditionals or constants for material properties.
#rho = 1.0*conditional(x[1] < 0.5*p1[1], 7800.0, 7000.0)
#kappa = 1.0*conditional(x[1] < 0.5*p1[1], 81.0, 121.0)
#cv = 1.0*rho*conditional(x[1] < 0.5*p1[1], 452.0, 387.0)

rho = 4500.0 # density in kg/m^3
kappa = 6.6 # thermal conductivity in W/(m-K)
cv = rho*650.0 # specific heat capacity in W/((m^3)-K)

# Heat source and prescribed heat fluxes
#r = Expression("1.0e6*sin(m*t)", m=2*pi/T, t=0, degree=1)

# Define a moving Gaussian heat source that simulates a heat source moving across the surface of the material.
# The heat source is defined as an exponential function that depends on the distance from the center of the Gaussian (x0, y0) and the time t. 
# The intensity of the heat source is scaled by a factor of 1.0e8 and modulated by a time-dependent term m*t, where m is set to 1/T to ensure
# that the heat source completes one full cycle over the total simulation time T.
r = Expression("1.0e8*exp(-((x[0]-x0)*(x[0]-x0)+(x[1]-y0)*(x[1]-y0))/0.001)*m*t", x0=0.125, y0=0.125, m=1/T, t=0.0, degree=4)
#q_p = Expression("1.0e3*m*t", m=1/T, t=0, degree=1)
#r = Constant(0.0)
q_p = Constant(0.0) # prescribed heat flux 

# Dirichlet boundary conditions for temperature 
bcs = [DirichletBC(V, Constant(0.0), boundaries, left), # apply zero temperature boundary condition on all boundaries
       DirichletBC(V, Constant(0.0), boundaries, right),
       DirichletBC(V, Constant(0.0), boundaries, bottom),
       DirichletBC(V, Constant(0.0), boundaries, top)]

#---------------------------------------------------------------------------------------------------------
#  Variational formulation of the Heat Equation (weak form)
#---------------------------------------------------------------------------------------------------------

# Crank-N. approximations for rates
# compute the change in temperature (dtheta) based on the previous temperature (theta) and the predicted temperature (theta_pred) using the Crank-Nicolson scheme.

def dtheta(th, th_pred):
	return (th-th_pred)/(dt*cn_alpha)

# Define a vector function space to store the heat flux
V_flux = VectorFunctionSpace(mesh, "Lagrange", 1)
q_actual = Function(V_flux)

# Form definition (F(u, delta_u) = a(u, delta_u) - l(delta_u))
# F= internal energy rate + Fourier's heat conduction + volumetric heat source - prescribed heat flux on the right boundary 

F = cv*dtheta(theta, theta_pred)*delta_theta*dx \
	+ kappa*dot(grad(delta_theta), grad(theta))*dx \
	- r*delta_theta*dx - q_p*delta_theta*ds(top)

# Create output files
theta_n.rename("theta","temperature")
xdmf_theta = XDMFFile("temperature_in_degrees.xdmf")
xdmf_theta.parameters["flush_output"] = True # ensures that data is written to file immediately after each write call
xdmf_theta.parameters["functions_share_mesh"] = True # allows multiple functions to share the same mesh in the output file, reducing file size and improving performance when writing multiple functions defined on the same mesh
xdmf_theta.write(theta_n, t)

# Rename for visualization
q_actual.rename("HeatFlux", "q")

# Save to file
xdmf_q = XDMFFile("heat_flux.xdmf")
xdmf_q.parameters["flush_output"] = True
xdmf_q.parameters["functions_share_mesh"] = True
xdmf_q.write(q_actual,t)
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
    t += dt
    print("Step ", n, " (t=", t, ")", sep="")

    # Update loads, boundary data, ...
    #u_p.t = t
    #theta_p.t = t
    #t_p.t = t
#    r.rx = (p1[0]-p0[0])*t/(0.01*T)
#    r.ry = (p1[1]-p0[1])*t/(0.01*T)

    # Define a moving Gaussian heat source that simulates a heat source moving across the surface of the material.
    DT = T/16 # simulation is divided into 16 segments 
    L = p1[0] # length of the domain in the x-direction
    DX = 0.25*L # distance from the center of the Gaussian to the boundaries
    r.t = T # set the time variable in the heat source expression to the current simulation time t, 
    # which allows the heat source to evolve over time according to its defined time dependence. 
    # This is important for simulating the movement of the heat source across the surface of the material, 
    # as the position of the heat source will change with time based on the defined expression for r.
    
    # The heat source moves in a square pattern across the surface of the material, 
    # starting from the center and moving towards the right boundary, 
    # then down to the bottom boundary, then left to the left boundary, 
    # and finally up to the top boundary before repeating the cycle. 
    # The position of the heat source is updated at each time step based on the current simulation time t 
    # and the defined time duration DT for each segment of the movement.
    if t<DT: 
        t0 = t/DT
        r.x0 = (1-t0)*0.25*L+t0*0.75*L # source starts at 0.25*L and moves to 0.75*L in the x-direction while staying at 0.25*L in the y-direction
    elif t<2*DT:
        t0 = (t-DT)/DT
        r.y0 = (1-t0)*0.25*L+t0*0.75*L # source moves from 0.25*L to 0.75*L in the y-direction while staying at 0.75*L in the x-direction
    elif t<3*DT:
        t0 = (t-2*DT)/DT
        r.x0 = t0*0.25*L+(1-t0)*0.75*L # source moves from 0.25*L to 0.75*L in the x-direction while staying at 0.25*L in the y-direction
    elif t<4*DT:
        t0 = (t-3*DT)/DT
        r.y0 = t0*0.25*L+(1-t0)*0.75*L # source moves from 0.25*L to 0.75*L in the y-direction while staying at 0.25*L in the x-direction
    #elif t<6*DT:
    #    pass
    else:
        r.m = 0
    #r.rx = 0.3*np.cos(2*pi*t/(0.3*T))+0.5
    #r.ry = 0.3*np.sin(2*pi*t/(0.4*T))+0.5
    #q_p.t = t
    
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
    
    # Compute actual heat flux after updating theta_n
    q_actual.assign(project(-kappa*grad(theta_n), V_flux))

    # Write step to files
    xdmf_theta.write(theta_n, t)
    xdmf_q.write(q_actual, t)

#-----------------------------------------------------------------------------------