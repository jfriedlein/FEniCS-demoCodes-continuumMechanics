#! /usr/bin/env python

from __future__ import print_function
from fenics import *
from os import path, mkdir

parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["optimize"] = True

"""-------------------------------------Time-dependent free vibration analysis of a 2D cantilever beam------------------------------------

Description: This analysis accounts for inertia, allowing observation of stress wave propagation 
and structural vibrations due to initial displacement along a selected mode (no external loads applied).

Geometry: 
A 2D rectangular solid  (Plane Stress formulation).

Material Model: Linear Isotropic elasticity.

Boundary Conditions:
Left Face: Fully clamped (Fixed: u = 0, v = 0).
Right Face: Free (no traction applied).

Loading:
No external loads or body forces. The system is excited by the initial displacement set along a selected eigenmode.

Main Learnings:
- Eigenvalue Analysis (Modal Analysis)
- Finds natural frequencies and mode shapes and filters them
- Sets up initial conditions for transient free vibration simulation based on a selected mode
- Performs time integration using the Newmark method

Possible Extensions:
- Introduce time-dependent or impulse loads for forced response
- Extend to 3D dynamic analysis
- Combine multiple modes for complex initial excitations
------------------------------------------------------------------------------------------------------------------------------ """


# Time definitions
t = 0.0 #intial time in seconds
T = 0.005 #final time in seconds
num_steps = 200 # number of time steps
dt = T / num_steps # time step size

#Parameters for Newmark (Average acceleration)
nm_beta = 0.25 # Newmark beta parameter
nm_gamma = 0.5 # Newmark gamma parameter

#---------------------------------------------------------------------------------------------------------
# Geometry and mesh generation
#---------------------------------------------------------------------------------------------------------


# Define geometry and mesh
p0 = Point(0.0, 0.0, 0.0) # corner point 
p1 = Point(0.4, 0.1, 0.1) # opposite corner point 
mesh = RectangleMesh(p0, p1, 32, 8) # apply structured mesh 32x8
#mesh = BoxMesh(p0, p1, 4, 1, 1)

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
x = SpatialCoordinate(mesh) # spatial coordinates

# Surface integral element
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

#-------------------------------------------------------------------------------------------------------
# Function spaces and variational functions
#-------------------------------------------------------------------------------------------------------

 
dim = 2 # spatial dimension
p = 2 # polynomial degree

V = VectorFunctionSpace(mesh, "Lagrange", p) # vector function space for displacements

# Define trial and test functions for variational formulation
u = TrialFunction(V) # trial
delta_u = TestFunction(V) # test 

# Interpolate displacement, velocity
u_n = Function(V) # previous displacement
u_n.interpolate(Constant((0.0, 0.0))) # set previous displacement is zero
v_n = Function(V) # previous velocity
v_n.interpolate(Constant((0.0, 0.0)))  # set previous velocity to zero

# Define additional rates and predictor functions for Newmark
a_n = Function(V) # previous acceleration
u_pred = Function(V) # predicted displacement
v_pred = Function(V) # predicted velocity


#--------------------------------------------------------------------------------------------------------
# Material properties
#--------------------------------------------------------------------------------------------------------

# Material parameters
E = 200.e9 # Young's modulus in Pa
nu = 0.3 # Poisson's ratio
mu    = E/(2.0*(1.0 + nu)) # shear modulus for plane strain # Pascals
lmbda = E*nu/((1.0 + nu)*(1.0 - 2.0*nu)) # first Lamé parameter for plane strain # Pascals
#mu    = E/(2.0*(1.0 + nu)) #PLANE STRESS
#lmbda = E*nu/((1.0 + nu)*(1.0 - nu)) #PLANE STRESS

rho = 8.e3 # mass density in kg/m^3
g = 9.81 # gravitational acceleration in m/s^2

# Volume force/ heat source and prescribed tractions
#b = as_vector((0.0, -rho*g)) 
#b = Constant((0.0, 0.0))
#t_p = Constant((0.0, 0.0))

# Prescribed Dirichlet boundary data
#u_p = Expression(('m*t', '0.0', '0.0'), degree=1, m=-0.2/T, t=0)
#theta_p = Expression(('m*t'), degree=1, m=100.0/T, t=0)

# Dirichlet boundary conditions
bcs = [DirichletBC(V, Constant((0.0, 0.0)), boundaries, left)] # fix left boundary

#--------------------------------------------------------------------------------------------------------
# Modal Analysis (Eigenvalue Problem)
#-------------------------------------------------------------------------------------------------------

# Strain tensor 
def epsilon(u):
      return sym(grad(u))

# Stress tensor (linear isotropic elasticity)
def sigma(u):
    return lmbda*tr(epsilon(u))*Identity(dim) + 2.0*mu*epsilon(u)

# Purpose: Find the natural frequencies and mode shapes of the structure 
# by solving the generalized eigenvalue problem K*rx = r*M*rx, 
# where K is the stiffness matrix, M is the mass matrix, 
# r are the eigenvalues (squared angular frequencies), and rx are the eigenvectors (mode shapes).

print("Assemble mass and stiffness matrix...")
m = rho*dot(delta_u, u)*dx # Bilinear form for mass matrix 
a = inner(grad(delta_u), sigma(u))*dx # Bilinear form for stiffness matrix

# Initialize PETSc matrices to store the assembled numerical data
M = PETScMatrix()
K = PETScMatrix()

#Transform the symbolic math (UFL) into numerical matrices
assemble(m, tensor=M)
assemble(a, tensor=K)

# Apply Dirichlet Boundary Conditions 
[bc.apply(M) for bc in bcs]
[bc.apply(K) for bc in bcs]

print("Solve for eigenvalues...")
# Setup the SLEPc solver for the generalized eigenvalue problem: K*rx = r*M*rx
eigensolver = SLEPcEigenSolver(K, M)  # K*rx = r*M*rx

# Create a directory and a PVD file to save all discovered vibration modes for ParaView
dirname2 = "all_evs"+"/"
if not path.exists(dirname2):
	mkdir(dirname2)
file_rx = File(dirname2+"modes.pvd", "compressed")

# Execute the solver to find all eigenvalues/eigenvectors
#eigensolver.parameters["problem_type"] = "gen_hermitian"
eigensolver.solve()

# Loop through all calculated pairs to filter and save valid physical modes
for i in range(0, M.size(0)):
	r, c, rx, cx = eigensolver.get_eigenpair(i) ## r = real eigenvalue (square of angular frequency - omega^2), rx = real eigenvector (displacement shape)
	if r < 1.1: # Filter out non-physical modes (e.g., rigid body modes with zero frequency)
		break
	print(" Eigenvalue ", i, " = +-i", sqrt(r))

    # Convert the raw eigenvector (rx) into a FEniCS Function for exporting
	rx_n = Function(V)
	rx_n.vector()[:]=rx
	rx_n.rename("u","mode")	
	file_rx << (rx_n, float(i)) # Save mode shape with its index as the timestamp

# Prompt user to select a specific mode for dynamic simulation
j = 0
print("Eigen value index taken: 0. You can change the value in main.py file line 184")

# Extract the selected eigenpair (frequency and mode shape) for the dynamic simulation
r, c, rx, cx = eigensolver.get_eigenpair(j)

# Set the initial displacement to the selected mode shape, scaled to a small amplitude for stability
u_n.vector()[:] = rx/norm(rx, "linf")*0.00001*0.5 # Scale the initial mode shape to a small amplitude (0.00001) for numerical stability in the dynamic simulation
#linf norm is used to find the maximum absolute value in the eigenvector rx, ensuring the mode shape is normalized before scaling.

T = 4*2*pi/sqrt(r) # Set total simulation time to cover multiple periods of the selected mode (4 full periods) to observe the dynamic response clearly. 
#The period is calculated as 2*pi/sqrt(r), where r is the eigenvalue corresponding to the selected mode.

dt = T / num_steps # set time step size

## Create a directory to save the dynamic simulation results for the selected mode
dirname = "ev"+str(j)+"/" 
if not path.exists(dirname):
	mkdir(dirname)

#---------------------------------------------------------------------------------------------------------
#  Variational formulation (weak form) (newmark)
#---------------------------------------------------------------------------------------------------------

# Newmark approximations for rates

# acceleration and velocity update functions based on Newmark-beta method
def a(u, u_pred):
	return (u-u_pred)/(dt*dt*nm_beta)
def v(u, u_pred, v_pred):
	return v_pred+dt*nm_gamma*a(u, u_pred)

# Form definition (F(u, delta_u) = a(u, delta_u) - l(delta_u))
F = rho*dot(delta_u, a(u, u_pred))*dx \
	+ inner(grad(delta_u), sigma(u))*dx #\
#	- dot(b, delta_u)*dx - dot(t_p, delta_u)*ds(top)

# Project initial stress field
Z = TensorFunctionSpace(mesh, "Lagrange", p) # Tensor function space for stress projection
stress_n = Function(Z) 
stress_n = project(sigma(u_n)/1.e6, Z, solver_type="mumps") # Project the initial stress field
	
# Create output files
u_n.rename("u","displacement")
xdmf_u = XDMFFile(dirname+"displacement_in_meters.xdmf")
xdmf_u.parameters["flush_output"] = True # ensures that data is written to file immediately after each write call
xdmf_u.parameters["functions_share_mesh"] = True # allows multiple functions to share the same mesh in the output file, reducing file size and improving performance when writing multiple functions defined on the same mesh

stress_n.rename("sigma","stress")
xdmf_stress = XDMFFile(dirname+"stress_in_MPa.xdmf")
xdmf_stress.parameters["flush_output"] = True
xdmf_stress.parameters["functions_share_mesh"] = True

xdmf_u.write(u_n, t)
xdmf_stress.write(stress_n, t)

#---------------------------------------------------------------------------------------------------------
# Time integration and function solving
#---------------------------------------------------------------------------------------------------------

# assemble bilinear form of variational formulation
A = assemble(lhs(F))
[bc.apply(A) for bc in bcs] # apply boundary conditions to A

# Set solver parameters
print("Setup solver...")
solver = LUSolver(A, "mumps") # Use MUMPS direct solver for efficiency and robustness in solving the linear systems at each time step
solver.parameters["symmetric"] = True # The stiffness matrix A is symmetric due to the nature of the elasticity problem
#solver.parameters["reuse_factorization"] = True # Reuse the symbolic factorization of the matrix structure across time steps, which can significantly speed up the solution process since the sparsity pattern of A does not change over time.

u = Function(V) # solution function to store the displacement at the current time step

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
    t += dt  # updates t every time step by dt
    print("Step ", n, " (t=", t, ")", sep="")
    
    # Update loads, boundary data, ...
    #u_p.t = t
    #theta_p.t = t
    #print("here")
	
    # Define predictors
    u_pred.vector()[:] = u_n.vector()+dt*v_n.vector()\
    					 +0.5*dt*dt*(1-2*nm_beta)*a_n.vector()
    v_pred.vector()[:] = v_n.vector()+dt*(1-nm_gamma)*a_n.vector()
    

    # Assemble and solve
    b = assemble(rhs(F)) # assemble the right-hand side vector b based on the current predicted displacement and velocity, as well as any time-dependent loads or boundary conditions. This step computes the load vector that corresponds to the current state of the system, which will be used in the linear solve to find the corrected displacement u at this time step.
    [bc.apply(b) for bc in bcs] # apply boundary conditions to b
    solver.solve(u.vector(), b) # solve the linear system A u = b to find the corrected displacement u at the current time step, using the assembled matrix A and the load vector b. The solution is stored in the function u, which represents the displacement field at this time step.
        
    # Update and correct functions
    #(u, theta) = y.split(deepcopy=True)
    
     # Compute and update current step acceleration using corrected displacement (Newmark formula)
    a_n.vector()[:] = (u.vector()-u_pred.vector())/(nm_beta*dt*dt)
    # Compute and update current step velocity using updated acceleration (Newmark formula)
    v_n.vector()[:] = v_pred.vector()+nm_gamma*dt*a_n.vector()
    #update displacement for the next time step
    u_n.vector()[:] = u.vector()
    
    # Project stresses at
    stress_n.assign(project(sigma(u_n)/1.e6, Z, solver_type="mumps")) # in MPa
    
    # Write step to files
    xdmf_u.write(u_n, t) # store displacement and stress at current time step in .xdmf file for visualization in ParaView
    xdmf_stress.write(stress_n, t)


#-----------------------------------------------------------------------------------------------------------------------------------------