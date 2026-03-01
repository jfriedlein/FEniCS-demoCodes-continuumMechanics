from __future__ import print_function
from fenics import *

#parameters['form_compiler']['optimize'] = True
#parameters['form_compiler']['cpp_optimize'] = True

"""-------------------------------------------------3D Thermo Elastostatics of a Rectangular Solid--------------------------------------------------
Problem description:
Geometry: Rectangular solid with dimensions 0.03 x 0.002 x 0.005
Boundary conditions:
	-Left face is fully clamped
	-All faces are held at 25 degree Celsius
      
Loads:
	- No surface traction is applied
    - No heat flux is applied 

	- No body forces applied
    - No Volumetric heat source is applied 


Analysis type: Quasi-static coupled thermoelastic model
Material model: Linear Isotropic thermal elastic material

Main learnings:
mesh.ufl_cell() - Usage of mesh.ufl_cell() for defining finite elements on the same cell type as the mesh
Mixed Function Space - Definition and usage of mixed function spaces for coupled problems
Defining and Solving variational formulation - Definition of combined variational formulation for coupled thermo-elastic problems
--------------------------------------------------------------------------------------------------------------------------------------------------

"""

#list_linear_solver_methods()

print("Preprocessing... ")


#---------------------------------------------------------------------------------------------------------
# Geometry and mesh generation
#---------------------------------------------------------------------------------------------------------

# spatial dimension
dim = 3 

# Define opposite corners of the rectangular domain
p0 = Point(0.0, 0.0, 0.0) #	bottom-left front corner
p1 = Point(0.03, 0.002, 0.005)  # top-right back corner
mesh = BoxMesh(p0, p1, 50, 10, 5) # Create a structured mesh of the rectangular solid

#---------------------------------------------------------------------------------------------------------
# Boundary identification and measures
#---------------------------------------------------------------------------------------------------------

# MeshFunction to store boundary IDs on facets
boundaries = MeshFunction("size_t", mesh, dim-1)
boundaries.set_all(0) # set every face to the value 0
left, right, bottom, top, back, front = 1, 2, 3, 4, 5, 6 # define boundary IDs

# filters every face based on its position and assigns the corresponding boundary ID
CompiledSubDomain("near(x[0], side) && on_boundary", side = p0[0]).mark(boundaries, left)
CompiledSubDomain("near(x[0], side) && on_boundary", side = p1[0]).mark(boundaries, right)
CompiledSubDomain("near(x[1], side) && on_boundary", side = p0[1]).mark(boundaries, bottom)
CompiledSubDomain("near(x[1], side) && on_boundary", side = p1[1]).mark(boundaries, top)
CompiledSubDomain("near(x[2], side) && on_boundary", side = p0[2]).mark(boundaries, back)
CompiledSubDomain("near(x[2], side) && on_boundary", side = p1[2]).mark(boundaries, front)

# Coordinate and surface integral element
x = SpatialCoordinate(mesh) # spatial coordinates for defining expressions over the domain - x[1], x[2], x[3]

#  Surface integral element
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

#-------------------------------------------------------------------------------------------------------
# Function spaces and variational functions
#-------------------------------------------------------------------------------------------------------

# Polynomial degree of Lagrange finite elements
p = 2

# Define mixed function space for displacement and temperature field.
P = VectorElement("Lagrange", mesh.ufl_cell(), p, 3) # vector value element for displacement
Q = FiniteElement("Lagrange", mesh.ufl_cell(), p) # scalar value element for temperature
V = FunctionSpace(mesh, P*Q) # mixed function space

# Define trial and test functions
(u, theta) = TrialFunctions(V)
(delta_u, delta_theta) = TestFunctions(V)


#--------------------------------------------------------------------------------------------------------
# Material properties
#--------------------------------------------------------------------------------------------------------

# Material parameters
E = 200.e9 # Young's modulus in Pascals
rho = 8.e3 # density in kg/m^3
g = 9.81 # acceleration due to gravity in m/s^2
nu = 0.3 # Poisson's ratio
mu    = E/(2.0*(1.0 + nu)) # Shear modulus # lame's second parameter
lmbda = E*nu/((1.0 + nu)*(1.0 - 2.0*nu)) # lame's first parameter
kappa = 80.0 # thermal conductivity in W/mK
alpha = 1.0*conditional(lt(x[1], p1[1]/2.0), 1.e-6, 15.0e-6) # thermal expansion coefficient in 1/K
# alpha = 1.e-6 for lower half, 15.e-6 for upper half


#--------------------------------------------------------------------------------------------------------
# Loads and boundary conditions
#--------------------------------------------------------------------------------------------------------

# Volume force  
b = Constant((0.0, 0.0, 0.0))

# heat source 
r = Constant(0.0)

# prescribed tractions
t_p = Constant((0.0, 0.0, 0.0))

# prescribed heat fluxes
q_p = Constant(0.0)

# Dirichlet boundary conditions
bcs = [DirichletBC(V.sub(0), Constant((0.0, 0.0, 0.0)), boundaries, left), # fully fixed left face
       DirichletBC(V.sub(1), Constant(25.0), boundaries, left), # All faces are held at 25 degree Celsius
       DirichletBC(V.sub(1), Constant(25.0), boundaries, right),
       DirichletBC(V.sub(1), Constant(25.0), boundaries, bottom),
       DirichletBC(V.sub(1), Constant(25.0), boundaries, top),
       DirichletBC(V.sub(1), Constant(25.0), boundaries, back),
       DirichletBC(V.sub(1), Constant(25.0), boundaries, front)
       ]


#---------------------------------------------------------------------------------------------------------
#  Variational formulation (weak form of heat conduction)
#---------------------------------------------------------------------------------------------------------

# Strain tensor
def epsilon(u):
    return sym(grad(u))

# Stress tensor (linear isotropic elasticity) with thermal strain
def sigma(u, theta):
    return lmbda*tr(epsilon(u))*Identity(3) + 2.0*mu*epsilon(u) - alpha*(3.0*lmbda + 2.0*mu)*Identity(3)*theta

# bilinear form: internal work + internal heat conduction
# left-hand side of the variational formulation
a = inner(grad(delta_u), sigma(u, theta))*dx + kappa*dot(grad(delta_theta), grad(theta))*dx

# linear form: external work + volumetric heat source + boundary heat flux
# right-hand side of the variational formulation
l = dot(b, delta_u)*dx + dot(t_p, delta_u)*ds(top) + r*delta_theta*dx + q_p*delta_theta*ds(top)

#---------------------------------------------------------------------------------------------------------
# Solve the linear system
#---------------------------------------------------------------------------------------------------------


y = Function(V)

#print("Assembling... ")
#A = assemble(a)
#B = assemble(l)
#for bc in bcs:
#	bc.apply(A, B)
#Y = y.vector()

print("Solving... ")
#solver = LUSolver(A, "mumps")
#solver.solve(Y, B)

# Solve Ku = f using the direct MUMPS solver
solve(a == l, y, bcs=bcs, 
	      solver_parameters={"linear_solver": "mumps"},
	      form_compiler_parameters={"optimize": True}) # solve the variational problem

(u, theta) = y.split()

#---------------------------------------------------------------------------------------------------------
#  Post processing:
#---------------------------------------------------------------------------------------------------------

print("Postprocessing... ")

# Create displacement and temperature file
u.rename("u","displacement") # rename displacement for output
File("displacement_in_meters.pvd", "compressed") << u # save displacement to file

# Compute and save stress
S = FunctionSpace(mesh, "Lagrange", p)  # scalar function space for stress component
stress = project(sigma(u, theta)[0,0]/1.e6, S) # stress component in MPa
stress.rename("sigma", "xx") # rename stress component for output

File("stress_in_MPa.pvd", "compressed") << stress # save stress to file 

# Create temperature file
theta.rename("theta", "temperature") # rename temperature for output
File("temperature_in_degrees.pvd", "compressed") << theta # save temperature to file

#------------------------------------------------------------------------------------------------------------------------------------------
