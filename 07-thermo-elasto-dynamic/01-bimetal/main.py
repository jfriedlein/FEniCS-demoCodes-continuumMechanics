#! /usr/bin/env python

from __future__ import print_function
from fenics import *

parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["optimize"] = True

"""--------------------------------------Coupled Thermo-Mechanical Analysis Using Newmark/Crank-Nicolson--------------------------------------

Problem Description:
This script performs a transient coupled thermo-mechanical simulation of a 2D rectangular domain using the Newmark method for structural dynamics and the Crank-Nicolson scheme for heat conduction. 
The simulation computes time-dependent displacements, velocities, accelerations, temperatures, and optionally stresses under prescribed boundary conditions and thermal loading.

Geometry:
- 2D rectangular domain: 0.03 m x 0.002 m x 0.005 m
- Structured mesh: 50 x 10 elements

Material Model:
- Linear isotropic elasticity:
    - Young's modulus (E): 200 GPa
    - Poisson's ratio (nu): 0.3
    - Density (rho): 8000 kg/m³
- Thermal properties:
    - Thermal conductivity (κ): 80 W/(m-K)
    - Specific heat capacity (cᵥ): 3.2 MW/((m^3)-K)
    - Thermal expansion coefficient (alpha): 1e-6 to 15e-6 (spatially varying)
    - Thermo-elastic coupling coefficient (beta)

Boundary Conditions:
- Dirichlet BCs:
    - Left boundary: fixed in x & y directions
    - Prescribed temperature θ_p on all boundaries (time-dependent)
- Neumann BCs:
    - Prescribed traction/heat flux (initialized to zero)
- Adiabatic/insulated conditions handled implicitly

Analysis :
- Newmark predictor-corrector scheme for displacement and velocity
- Crank-Nicolson scheme for temperature evolution
- Coupling between thermal expansion and elastic response

Main Learnings:
- Coupled thermo-mechanical transient response under time-dependent thermal loading
- Handling time-dependent Dirichlet boundary conditions
- Predictor-corrector Newmark method with Crank-Nicolson temperature update

---------------------------------------------------------------------------------------------------------------------------------------------"""

# Time definitions
t = 0.0 # initial time in seconds
T = 10 # final time in seconds
num_steps = 200 # number of time steps
dt = T / num_steps # time step size

# Parameters for Newmark (Average acceleration)/Crank-Nicolson
cn_alpha = 0.5 # alpha parameter for Crank-Nicolson method
nm_beta = 0.25 # beta parameter for Newmark method
nm_gamma = 0.5 # gamma parameter for Newmark method

#---------------------------------------------------------------------------------------------------------
# Geometry and mesh generation
#---------------------------------------------------------------------------------------------------------

# Define geometry and mesh
p0 = Point(0.0, 0.0, 0.0) # lower left corner 
p1 = Point(0.03, 0.002, 0.005) # upper right corner
mesh = RectangleMesh(p0, p1, 50, 10) # apply a structured mesh to the rectangle domain

#---------------------------------------------------------------------------------------------------------
# Boundary identification and marking
#---------------------------------------------------------------------------------------------------------

# MeshFunction to store boundary IDs on facets
boundaries = MeshFunction("size_t", mesh, 1) #function to store boundary IDs on facets
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

dim = 2 
p = 2 # polynomial degree for finite element space
P = VectorElement("Lagrange", mesh.ufl_cell(), p, dim) # vector-valued finite element space for displacement
Q = FiniteElement("Lagrange", mesh.ufl_cell(), p) # scalar-valued finite element space for temperature
V = FunctionSpace(mesh, P*Q) # mixed function space for displacement and temperature
# mesh.ufl_cell() returns the type of cell used in the mesh (e.g., triangle, quadrilateral, tetrahedron, etc.) and is used to define the finite element spaces appropriately.

# Define trial and test functions
(u, theta) = TrialFunctions(V) # trial functions for displacement and temperature
(delta_u, delta_theta) = TestFunctions(V) # test functions for displacement and temperature

# Collapse mixed function space
V_u = V.sub(0).collapse() # extract the subspace for displacement and collapse it to a standard function space
V_theta = V.sub(1).collapse() # extract the subspace for temperature and collapse it to a standard function space

# Interpolate initial displacement, velocity and temperature
u_n = Function(V_u) # function to store the displacement at the previous time step, initialized to zero
u_n.interpolate(Constant((0.0, 0.0))) 
#u_n.interpolate(Expression(("0.0", "x[0]*x[0]/l/l*0.1"), l=p1[0], degree=1)) # initial displacement with a parabolic profile in x-direction, scaled by 0.1 for visualization purposes
v_n = Function(V_u) # function to store the velocity at the previous time step, initialized to zero
v_n.interpolate(Constant((0.0, 0.0)))

theta0 = 273.0 # reference temperature in Kelvin used to define the zero-thermal-stress state and to scale the heat generated by mechanical strain rates.
theta_n = Function(V_theta) # function to store the temperature at the previous time step
theta_n.interpolate(Constant(0.0)) # initial temperature is zero everywhere 

# Define additional rates and predictor functions for Newmark
a_n = Function(V_u) # function to store the acceleration at the previous time step
dtheta_n = Function(V_theta) # function to store the rate of change of temperature at the previous time step
u_pred = Function(V_u) # function to store the predicted displacement for the current time step
v_pred = Function(V_u) # function to store the predicted velocity for the current time step
theta_pred = Function(V_theta) # function to store the predicted temperature for the current time step

#--------------------------------------------------------------------------------------------------------
# Material properties
#--------------------------------------------------------------------------------------------------------

E = 200.e9 # Young's modulus in Pascals
nu = 0.3 # Poisson's ratio
mu    = E/(2.0*(1.0 + nu)) # shear modulus
lmbda = E*nu/((1.0 + nu)*(1.0 - 2.0*nu)) # first Lamé parameter 
#mu    = E/(2.0*(1.0 + nu)) #PLANE STRESS
#lmbda = E*nu/((1.0 + nu)*(1.0 - nu)) #PLANE STRESS
rho = 8.e3 # mass density in kg/m^3
g = 9.81 # gravitational acceleration in m/s^2

kappa = 80.0 # thermal conductivity in W/(m*K)

# Thermal expansion coefficient (alpha) and thermo-elastic coupling coefficient (beta)
alpha = 1.0*conditional(lt(x[1], p1[1]/2.0), 1.e-6, 15.0e-6) 
beta = alpha*E/(1-2.0*nu)
cv = rho*400 # specific heat capacity in MW/((m^3)*K)

# Volume force/ heat source and prescribed tractions/ prescribed heat fluxes
#b = as_vector((0.0, -rho*g, 0.0))
b = Constant((0.0, 0.0)) # no body force for simplicity
t_p = Constant((0.0, 0.0)) # prescribed traction on the top boundary, initialized to zero
#t_p = Expression(('(t<1.0*t1)?1.0*(-m/t1*fabs(t-t1)+m):0.0', '(t<1.0*t1)?-0.1*(-m/t1*fabs(t-t1)+m):0.0'), degree=2, t1=0.001, m=1.e9, t=0)
r = Constant(0.0) # volumetric heat source, initialized to zero
q_p = Constant(0.0) # prescribed heat flux, initialized to zero

# Prescribed Dirichlet boundary data
#u_p = Expression("u0", degree=1, u0=0.0)
theta_p = Expression(("t0+m*t"), degree=1, m=25.0/T, t0=0, t=0) # temperature in Kelvin prescribed on all boundaries, 
# increasing linearly with time from t0 to t0+m at the final time T. The parameter m controls the rate of increase of temperature, 
# and t is the current simulation time that will be updated in the time-stepping loop.

# Dirichlet boundary conditions
bcs = [DirichletBC(V.sub(0), Constant((0.0, 0.0)), boundaries, left), # fix left boundary in both x and y directions
       DirichletBC(V.sub(1), theta_p, boundaries, left), # prescribe temperature on left boundary
       DirichletBC(V.sub(1), theta_p, boundaries, right), # prescribe temperature on right boundary
       DirichletBC(V.sub(1), theta_p, boundaries, bottom), # prescribe temperature on bottom boundary
       DirichletBC(V.sub(1), theta_p, boundaries, top) # prescribe temperature on top boundary
      ]

#---------------------------------------------------------------------------------------------------------
#  Variational formulation of the Heat Equation (weak form)
#---------------------------------------------------------------------------------------------------------
# Strain tensor
def epsilon(u):
    return sym(grad(u))

# Stress tensor (linear isotropic elasticity)
def sigma(u, theta):
    return lmbda*tr(epsilon(u))*Identity(dim) + 2.0*mu*epsilon(u) - alpha*(3.0*lmbda + 2.0*mu)*Identity(dim)*theta

# Newmark approximations for rates
# compute the acceleration, velocity and rate of change of temperature based on the Newmark predictor values and the current solution.
def a(u, u_pred):
	return (u-u_pred)/(dt*dt*nm_beta)
def v(u, u_pred, v_pred):
	return v_pred+dt*nm_gamma*a(u, u_pred)
def dtheta(th, th_pred):
	return (th-th_pred)/(dt*cn_alpha)

# Form definition (F(u, delta_u) = a(u, delta_u) - l(delta_u))
# F= internal energy rate + kinetic energy rate + Fourier's heat conduction + volumetric heat source - prescribed heat flux on the top boundary - work done by body forces and tractions
F = rho*dot(delta_u, a(u, u_pred))*dx \
	+ cv*dtheta(theta, theta_pred)*delta_theta*dx \
	+ inner(grad(delta_u), sigma(u, theta))*dx \
	+ kappa*dot(grad(delta_theta), grad(theta))*dx \
	+ theta0*beta*tr(sym(grad(v(u, u_pred, v_pred))))*delta_theta*dx \
	- dot(b, delta_u)*dx - dot(t_p, delta_u)*ds(top) \
	- r*delta_theta*dx - q_p*delta_theta*ds(top)

# Project initial stress field
#V_stress = TensorFunctionSpace(mesh, "Lagrange", p)
#stress_n = Function(V_stress)
#stress_n = project(sigma(u_n, theta_n)/1.e6, V_stress, solver_type="mumps")

# Create output files
u_n.rename("u", "displacement") # rename the function u_n to "u" 
xdmf_u = XDMFFile("displacement_in_meters.xdmf")
xdmf_u.parameters["flush_output"] = True # ensures that data is written to file immediately after each write call
xdmf_u.parameters["functions_share_mesh"] = True # allows multiple functions to share the same mesh in the output file, reducing file size and improving performance when writing multiple functions defined on the same mesh
xdmf_u.write(u_n, t)

xdmf_theta = XDMFFile("temperature_in_Kelvin.xdmf")
xdmf_theta.parameters["flush_output"] = True # ensures that data is written to file immediately after each write call
xdmf_theta.parameters["functions_share_mesh"] = True # allows multiple functions to share the same mesh in the output file, reducing file size and improving performance when writing multiple functions defined on the same mesh
xdmf_theta.write(theta_n, t)
#stress_n.rename("stress","vonMises")
#xdmf_stress = XDMFFile("stress.xdmf")
#xdmf_stress.write(stress_n, t)

#---------------------------------------------------------------------------------------------------------
# Time integration and function solving
#---------------------------------------------------------------------------------------------------------

A = assemble(lhs(F)) # assemble the left-hand side of the variational form F
[bc.apply(A) for bc in bcs] # apply the boundary conditions to the matrix A, modifying it to account for the Dirichlet conditions specified in bcs. This typically involves setting rows corresponding to Dirichlet boundary conditions to enforce the prescribed values and ensuring that the system of equations is consistent with 

print("Setup solver...")
solver = LUSolver(A, "mumps")
solver.parameters["symmetric"] = True # The stiffness matrix A is symmetric due to the nature of the thermal diffusion operator
#solver.parameters["reuse_factorization"] = True # reuse the symbolic factorization of the matrix A across time steps, 
# which can significantly reduce computational cost when the structure of A does not change over time 
# (as is the case here since A is derived from a linear variational form that does not depend on the solution at previous time steps).

y = Function(V) # function to store the solution at the current time step, defined on the mixed function space V, which includes both displacement and temperature components. This function will be updated at each time step with the new solution obtained from solving the linear system defined by the variational form F.

# Time integration and function solving

""" 

1. Increment simulation time t += dt, update temperature theta_p.t = t for the current timestep
2. Predict displacement, temperature and velocity at the current timestep using Newmark formulas:
        Use previous timestep values (u_n, theta_n, v_n, a_n) to compute u_pred and vel_pred
3. Assemble and solve the linear system
        Solve A y = b to obtain the corrected displacement u and temperature theta at this timestep and split the solution into its components (u, theta)
4. Update acceleration and velocity using the corrected displacement and update the rate of change of temperature using the corrected temperature:
5. Update displacement and temperature for the next timestep: u_n = u, theta_n = theta
6. If desired, compute and project stresses onto the tensor function space in MPa:
7. Write results to output files (.xdmf) for displacement and stress.
Repeat for all timesteps (num_steps)

"""
for n in range(num_steps):	
    t += dt # increment simulation time by the time step size dt
    print("Step ", n, " (t=", t, ")", sep="")
    
    # Update loads, boundary data, ...    
    theta_p.t = t # update the time variable t in the expression theta_p to reflect the current simulation time, ensuring that the prescribed temperature boundary condition changes appropriately at each time step according to the defined expression for theta_p.
    #if t > T/2:
    #    theta_p.m = 0

    # Define predictors for the current time step using Newmark formulas
    u_pred.vector()[:] = u_n.vector()+dt*v_n.vector()\
    					 +0.5*dt*dt*(1-2*nm_beta)*a_n.vector()
    v_pred.vector()[:] = v_n.vector()+dt*(1-nm_gamma)*a_n.vector()
    theta_pred.vector()[:] = theta_n.vector()+dt*(1-cn_alpha)*dtheta_n.vector()
    
    # Assemble and solve
    b = assemble(rhs(F))
    [bc.apply(b) for bc in bcs] # apply the boundary conditions to the right-hand side vector b, modifying it to account for the Dirichlet conditions specified in bcs. This typically involves setting entries corresponding to Dirichlet boundary conditions to enforce the prescribed values and ensuring that the system of equations is consistent with the modified matrix A.
    solver.solve(y.vector(), b) # solve the linear system A y = b to obtain the corrected solution y at the current time step, which includes both displacement and temperature components. The solution is stored in the vector of the function y.
        
    # Update and correct functions
    (u, theta) = y.split(deepcopy=True) # split the solution y into its components u (displacement) and theta (temperature). 
    #The deepcopy=True argument ensures that the resulting functions u and theta are independent copies of the corresponding components of y, 
    # allowing for modifications without affecting the original solution y.
    a_n.vector()[:] = (u.vector()-u_pred.vector())/(nm_beta*dt*dt) # update acceleration
    v_n.vector()[:] = v_pred.vector()+nm_gamma*dt*a_n.vector() # update velocity
    u_n.vector()[:] = u.vector() # update displacement for the next time step
    dtheta_n.vector()[:] = (theta.vector()-theta_pred.vector())/(cn_alpha*dt)    # update rate of change of temperature
    theta_n.vector()[:] = theta.vector() # update temperature for the next time step
    
    # Project stresses
#    stress_n.assign(project(sigma(u_n, theta_n)/1.e6, V_stress, solver_type="mumps"))
    
    # Write step to files
    xdmf_u.write(u_n, t) # write the updated displacement u_n to the XDMF file for displacement.
    xdmf_theta.write(theta_n, t) # write the updated temperature theta_n to the XDMF file for temperature.
#   xdmf_stress.write(stress_n, t)

#------------------------------------------------------------------------------------------------------------------------------------------