#! /usr/bin/env python

from __future__ import print_function

from sympy import Function
from fenics import *

parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["optimize"] = True

"""----------------------------------Transient Thermo-Mechanical Analysis Using Newmark/Crank-Nicolson----------------------------------

Problem Overview:
This script performs a 2D coupled thermo-mechanical transient analysis of a rectangular domain using:
- Newmark method for structural dynamics (displacement, velocity, acceleration)
- Crank-Nicolson scheme for transient heat conduction (temperature evolution)

Geometry & Mesh:
- 2D rectangular domain: 0.4 m x 0.1 m
- Structured mesh: 100 x 25 quadrilateral elements

Material Model:
- Linear isotropic elasticity
  - Young's modulus (E): 200 GPa
  - Poisson's ratio (nu): 0.3
  - Density (rho): 8000 kg/m³
- Thermal properties
  - Thermal conductivity (κ): 80 W/(m-K)
  - Specific heat (cᵥ): 3.2 W/((m^3)-K)
  - Thermal expansion coefficient (alpha): 12 x 10⁻⁶ 1/K
  - Thermo-mechanical coupling coefficient (beta)

Boundary Conditions:
- Dirichlet BCs: fixed on left boundary
- Neumann BCs: time-dependent sinusoidal traction on right boundary
- Heat flux: zero prescribed flux

Loading:
Time-dependent traction: oscillates sinusoidally in the x-direction with magnitude over period .


Main Learnings:
- Handling time-dependent boundary conditions
- Coupling between mechanical deformation and thermal response

"""
# Time definitions
t = 0.0 # initial time in seconds
T = 1.0 # final time in seconds
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
p1 = Point(0.4, 0.1, 0.1) # upper right corner
mesh = RectangleMesh(p0, p1, 100, 25) # apply a structured mesh to the rectangle domain

#---------------------------------------------------------------------------------------------------------
# Boundary identification and marking
#---------------------------------------------------------------------------------------------------------


boundaries = MeshFunction("size_t", mesh, 1) # MeshFunction to store boundary IDs on facets
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

# Define function space
dim = 2
p = 2  # polynomial degree
P = VectorElement("Lagrange", mesh.ufl_cell(), p, dim) # vector element for displacement
Q = FiniteElement("Lagrange", mesh.ufl_cell(), p) # scalar element for temperature
V = FunctionSpace(mesh, P*Q) # mixed function space for displacement and temperature

# Define trial and test functions
(u, theta) = TrialFunctions(V) # trial functions for displacement and temperature
(delta_u, delta_theta) = TestFunctions(V) # test functions for virtual work and virtual heat transfer

# Collapse mixed function space
V_u = V.sub(0).collapse() # function space for displacement
V_theta = V.sub(1).collapse() # function space for temperature

# Interpolate initial displacement, velocity and temperature
u_n = Function(V_u) # function to store displacement at previous time step
u_n.interpolate(Constant((0.0, 0.0))) # initial displacement is zero everywhere
v_n = Function(V_u) # function to store velocity at previous time step
v_n.interpolate(Constant((0.0, 0.0)))  # initial velocity is zero everywhere
#v_n.interpolate(Expression("x[0]*x[0]/l/l*0.1"), l=p1[0], degree=1)) # if desired, prescribe an initial velocity field (e.g., quadratic distribution in x-direction)

theta0 = 298 # reference temperature in Kelvin used to define the zero-thermal-stress state and to scale the heat generated by mechanical strain rates.
theta_n = Function(V_theta) # function to store temperature at previous time step
theta_n.interpolate(Constant(0.0)) # initial temperature is zero everywhere

# Define additional rates and predictor functions for Newmark
a_n = Function(V_u) # function to store acceleration at previous time step
dtheta_n = Function(V_theta) # function to store rate of change of temperature at previous time step
u_pred = Function(V_u) # function to store predicted displacement at current time step based on previous time step values and Newmark formulas
v_pred = Function(V_u) # function to store predicted velocity at current time step based on previous time step values and Newmark formulas
theta_pred = Function(V_theta) # function to store predicted temperature at current time step based on previous time step values and Crank-Nicolson formula

#--------------------------------------------------------------------------------------------------------
# Material properties
#--------------------------------------------------------------------------------------------------------

# Material parameters
E = 200.e9 # Young's modulus in Pascals
nu = 0.  # Poisson's ratio
mu    = E/(2.0*(1.0 + nu)) # shear modulus # second Lamé parameter
lmbda = E*nu/((1.0 + nu)*(1.0 - 2.0*nu)) # first Lamé parameter 
#mu    = E/(2.0*(1.0 + nu)) #PLANE STRESS
#lmbda = E*nu/((1.0 + nu)*(1.0 - nu)) #PLANE STRESS
rho = 8.e3 # mass density in kg/m^3
g = 9.81 # gravitational acceleration in m/s^2

kappa = 80.0 # thermal conductivity in W/(m*K)
alpha = 12.e-6 # thermal expansion coefficient in 1/K
#alpha = 0.0
beta = alpha*E/(1-2.0*nu) # thermo-mechanical coupling coefficient in W/(m^3*K) 
# scales the heat generated by mechanical strain rates, derived from the thermal expansion coefficient and elastic constants. It represents the conversion of mechanical work into heat due to thermoelastic effects.
cv = rho*400 # specific heat capacity in J/(kg*K)

# Volume force/ heat source and prescribed tractions/ prescribed heat fluxes
#b = as_vector((0.0, -rho*g, 0.0))
b = Constant((0.0, 0.0)) # body force per unit volume 
#t_p = Constant((0.0, 0.0)) 
t_p = Expression(("m*sin(2*3.141*t/t0)", "0.0"), degree=4, t0=T, m=1.e8, t=0) # Time-dependent traction: oscillates sinusoidally in the x-direction with magnitude 'm' over period 't0'.
#t_p = Expression(("(t<3.0*t1)?1.0*(-m/t1*fabs(t-t1)+m):0.0", "0.0"), degree=2, t1=0.25*T, m=1.e9, t=0)
#t_p = Expression(("(t<t1)?m:0.0", "0.0"), degree=2, t1=0.33*T, m=1.e9, t=0)
r = Constant(0.0) # volumetric heat source per unit volume in W/m^3
q_p = Constant(0.0) # prescribed heat flux per unit area in W/m^2

# Prescribed Dirichlet boundary data
#u_p = Expression(("0.03*sin(2*3.141*t/t0)", "0.0"), t=t, t0=T, degree=3)
#theta_p = Expression(('m*t'), degree=1, m=100.0/T, t=0)

# Dirichlet boundary conditions
bcs = [DirichletBC(V.sub(0), Constant((0.0, 0.0)), boundaries, left), # prescribe zero displacement on the left boundary
       #DirichletBC(V.sub(0), u_p, boundaries, right),
       #DirichletBC(V.sub(0).sub(0), u0, boundaries, right),
       #DirichletBC(V.sub(1), theta_p, boundaries, left)]#,
       #DirichletBC(V.sub(1), Constant(theta0), boundaries, right)
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


#dirname = "ev"+str(j)+"/"
#if not path.exists(dirname):
#	mkdir(dirname)


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
	- dot(b, delta_u)*dx - dot(t_p, delta_u)*ds(right) \
	- r*delta_theta*dx - q_p*delta_theta*ds(right)

# Project initial stress field
V_stress = TensorFunctionSpace(mesh, "Lagrange", p)
stress_n = Function(V_stress)
stress_n = project(sigma(u_n, theta_n)/1.e6, V_stress, solver_type="mumps")

# Define a vector function space to store the heat flux
V_flux = VectorFunctionSpace(mesh, "Lagrange", 1)
q_actual = Function(V_flux)

# Create output files
u_n.rename("u", "displacement") # rename the function u_n to "u" 
xdmf_u = XDMFFile("displacement_in_meters.xdmf")
xdmf_u.parameters["flush_output"] = True # ensures that data is written to file immediately after each write call
xdmf_u.parameters["functions_share_mesh"] = True # allows multiple functions to share the same mesh in the output file, reducing file size and improving performance when writing multiple functions defined on the same mesh
xdmf_u.write(u_n, t)

theta_n.rename("theta", "temperature")
xdmf_theta = XDMFFile("temperature_in_Kelvin.xdmf")
xdmf_theta.parameters["flush_output"] = True # ensures that data is written to file immediately after each write call
xdmf_theta.parameters["functions_share_mesh"] = True # allows multiple functions to share the same mesh in the output file, reducing file size and improving performance when writing multiple functions defined on the same mesh
xdmf_theta.write(theta_n, t)

stress_n.rename("stress","vonMises")
xdmf_stress = XDMFFile("stress.xdmf")
xdmf_stress.write(stress_n, t)

# Rename for visualization
q_actual.rename("HeatFlux", "q")

# Save to file (XDMF)
xdmf_q = XDMFFile("heat_flux.xdmf")
xdmf_q.parameters["flush_output"] = True
xdmf_q.parameters["functions_share_mesh"] = True
xdmf_q.write(q_actual, t)


#---------------------------------------------------------------------------------------------------------
# Time integration and function solving
#---------------------------------------------------------------------------------------------------------

A = assemble(lhs(F)) # assemble the left-hand side of the variational form F
[bc.apply(A) for bc in bcs] # apply the boundary conditions to the matrix A, modifying it to account for the Dirichlet conditions specified in bcs. This typically involves setting rows corresponding to Dirichlet boundary conditions to enforce the prescribed values and ensuring that the system of equations is consistent with 

print("Setup solver...")
solver = LUSolver(A, "mumps")
#solver.parameters["symmetric"] = True
#solver.parameters["reuse_factorization"] = True
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
    t += dt
    print("Step ", n, " (t=", t, ")", sep="")
    
    # Update loads, boundary data, ...
    #u_p.t = t
    #theta_p.t = t
    t_p.t = t

    # Define predictors 
    u_pred.vector()[:] = u_n.vector()+dt*v_n.vector()\
    					 +0.5*dt*dt*(1-2*nm_beta)*a_n.vector()
    v_pred.vector()[:] = v_n.vector()+dt*(1-nm_gamma)*a_n.vector()
    theta_pred.vector()[:] = theta_n.vector()+dt*(1-cn_alpha)*dtheta_n.vector()
    
    # Assemble and solve
    b = assemble(rhs(F))
    [bc.apply(b) for bc in bcs] # apply the boundary conditions to the right-hand side vector b, modifying it to account for the Dirichlet conditions specified in bcs. This typically involves setting entries corresponding to Dirichlet boundary conditions to enforce the prescribed values and ensuring that the system of equations is consistent with the modified matrix A.
    solver.solve(y.vector(), b) # solve the linear system A y = b to obtain the corrected solution y at the current time step, which includes both displacement and temperature components. The solution is stored in the vector of the function y.
    
        
    # Update and correct functions # split the solution y into its components u (displacement) and theta (temperature). 
    #The deepcopy=True argument ensures that the resulting functions u and theta are independent copies of the corresponding components of y, 
    # allowing for modifications without affecting the original solution y.
    (u, theta) = y.split(deepcopy=True)
    a_n.vector()[:] = (u.vector()-u_pred.vector())/(nm_beta*dt*dt) # update acceleration 
    v_n.vector()[:] = v_pred.vector()+nm_gamma*dt*a_n.vector() # update velocity
    u_n.vector()[:] = u.vector() # update displacement for the next time step
    dtheta_n.vector()[:] = (theta.vector()-theta_pred.vector())/(cn_alpha*dt)  # update rate of change of temperature  
    theta_n.vector()[:] = theta.vector() # update temperature for the next time step
    
    # Project stresses
    stress_n.assign(project(sigma(u_n, theta_n)/1.e6, V_stress, solver_type="mumps"))
    
     # Compute actual heat flux after updating theta_n
    q_actual.assign(project(-kappa*grad(theta_n), V_flux))

    # Write step to files
    xdmf_u.write(u_n, t) # write the updated displacement u_n to the XDMF file for displacement.
    xdmf_theta.write(theta_n, t) # write the updated temperature theta_n to the XDMF file for temperature.
    xdmf_stress.write(stress_n, t)
    xdmf_q.write(q_actual, t)

#-------------------------------------------------------------------------------------------------------------------------------