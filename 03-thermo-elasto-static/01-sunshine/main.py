from __future__ import print_function
from fenics import *

#parameters['form_compiler']['optimize'] = True
#parameters['form_compiler']['cpp_optimize'] = True

"""-------------------------------------------------3D Thermo Elastostatics of a Rectangular Solid--------------------------------------------------
Problem description:
Geometry: Rectangular solid with dimensions 0.3 x 0.1 x 0.1
Boundary conditions:
	-Left  and Right face is fully clamped
	-Left and Right face are held at 0 K
      
Loads:
	- No surface traction is applied
    - A prescribed heat flux is applied on the top face and bottom face

	- No body forces applied
    - No Volumetric heat source is applied over the domain


      
Analysis type: Quasi-static coupled thermoelastic model
Material model: Linear Isotropic thermal elastic material

Main learnings:
mesh.ufl_cell() - Usage of mesh.ufl_cell() for defining finite elements on the same cell type as the mesh
Mixed Function Space - Definition and usage of mixed function spaces for coupled problems
Defining and Solving variational formulation - Definition of combined variational formulation for coupled thermo-elastic problems

Possibilities for extensions:
heat convection and radiation to environment
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
p1 = Point(0.30, 0.10, 0.10)  # top-right back corner
mesh = BoxMesh(p0, p1, 30, 10, 10) # Create a structured mesh of the rectangular solid

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
# Create a vector-valued Lagrange finite element of degree p on the same cell type as the mesh
# mesh.ufl_cell() tells UFL what kind of reference cell (triangle, tetrahedron, etc.) your mesh uses so the finite element is defined correctly.

Q = FiniteElement("Lagrange", mesh.ufl_cell(), p) # scalar value element for temperature
# Create a scalar-valued Lagrange finite element of degree p on the same cell type as the mesh
# mesh.ufl_cell() tells UFL what kind of reference cell (triangle, tetrahedron, etc.) your mesh uses so the finite element is defined correctly.

V = FunctionSpace(mesh, P*Q) # mixed function space
# Thermo-elasticity is coupled:
#       Temperature affects displacement through thermal strain
#       Displacement and temperature are solved simultaneously

# Define trial and test functions
(u, theta) = TrialFunctions(V)
(delta_u, delta_theta) = TestFunctions(V)


#--------------------------------------------------------------------------------------------------------
# Material properties
#--------------------------------------------------------------------------------------------------------

E = 200.e9 # Young's modulus in Pascals
rho = 8.e3 # density in kg/m^3
g = 9.81 # acceleration due to gravity in m/s^2
nu = 0.3 # Poisson's ratio
mu    = E/(2.0*(1.0 + nu)) # shear modulus # second Lamé parameter # Pascals
lmbda = E*nu/((1.0 + nu)*(1.0 - 2.0*nu)) # first Lamé parameter # Pascals
kappa = 80.0 # thermal conductivity in W/(m·K)
alpha = 15.0e-6 # coefficient of linear thermal expansion in 1/K


#--------------------------------------------------------------------------------------------------------
# Loads and boundary conditions
#--------------------------------------------------------------------------------------------------------

# Volume force in N/m^3
b = Constant((0.0, 0.0, 0.0))

# heat source in W/m^3
r = Constant(0.0)

# prescribed tractions in Pascals
t_p = Constant((0.0, 0.0, 0.0))

# prescribed heat fluxes in W/m^2
#if heat flux points inwards (q_p_top), it points opposite to surface normal -> positive heat flux (energy gain for system, system egoistic viewpoint)
#if heat flux points outwards (q_p_bottom), it points in the direction of the surface normal -> negative heat flux (energy loss from the system's point of view)
q_p_top = Constant(10000.0) 
q_p_bottom = Constant(-10000.0)

# reference temperature in Kelvin
theta_0=0

# Dirichlet boundary conditions
bcs = [DirichletBC(V.sub(0), Constant((0.0, 0.0, 0.0)), boundaries, left), # fully fixed left face
	   DirichletBC(V.sub(0), Constant((0.0, 0.0, 0.0)), boundaries, right), # fully fixed right face
       DirichletBC(V.sub(1), Constant(theta_0), boundaries, left), # left face at 0 Kelvin
       DirichletBC(V.sub(1), Constant(theta_0), boundaries, right), # right face at 0 Kelvin
       #DirichletBC(V.sub(1), Constant(20.0), boundaries, bottom),
       #DirichletBC(V.sub(1), Constant(20.0), boundaries, top),
       #DirichletBC(V.sub(1), Constant(20.0), boundaries, back),
       #DirichletBC(V.sub(1), Constant(20.0), boundaries, front)
       ]


#---------------------------------------------------------------------------------------------------------
#  Variational formulation (weak form of heat conduction)
#---------------------------------------------------------------------------------------------------------

# Strain tensor
def epsilon(u):
    return sym(grad(u))

# Stress tensor for isotropic linear elasticity with thermal strain
def sigma(u, theta):
    return lmbda*tr(epsilon(u))*Identity(3) + 2.0*mu*epsilon(u) - alpha*(3.0*lmbda + 2.0*mu)*Identity(3)*theta

# bilinear form: internal work + internal heat conduction
# left-hand side of the variational formulation
a = inner(grad(delta_u), sigma(u, theta))*dx + kappa*dot(grad(delta_theta), grad(theta))*dx


# linear form: external virtual work + volumetric heat source + boundary heat flux
# right-hand side of the variational formulation
l = dot(b, delta_u)*dx + dot(t_p, delta_u)*ds(top) + r*delta_theta*dx + q_p_top*delta_theta*ds(top) + q_p_bottom*delta_theta*ds(bottom)


#---------------------------------------------------------------------------------------------------------
# Solve the linear system
#---------------------------------------------------------------------------------------------------------

# Solution function to store displacement and temperature field
y = Function(V)

# Solve Ku = f using the direct MUMPS solver
solve(a == l, y, bcs=bcs, 
	      solver_parameters={"linear_solver": "mumps"},
	      form_compiler_parameters={"optimize": True}) # solve the variational problem

#print("Assembling... ")
#A = assemble(a)
#B = assemble(l)
#for bc in bcs:
#	bc.apply(A, B)
#Y = y.vector()

#print("Solving... ")
#solver = LUSolver(A, "mumps")
#solver.solve(Y, B)

(u, theta) = y.split()
# y.split() splits the mixed solution function y into its components

#---------------------------------------------------------------------------------------------------------
#  Post processing:
#---------------------------------------------------------------------------------------------------------

print("Postprocessing... ")

# Create displacement and temperature file
u.rename("u","displacement") # rename displacement for output
File("displacement_in_meters.pvd", "compressed") << u # save displacement to file

# Define function space for stress component
S = FunctionSpace(mesh, "Lagrange", p)
stress = project(sigma(u, theta)[0,0]/1.e6, S)  # stress component in MPa
stress.rename("sigma", "xx") # rename stress component for output

File("stress_in_MPa.pvd", "compressed") << stress # save stress to file

# Create temperature file
theta.rename("theta", "temperature") # rename temperature for output
File("temperature_in_Kelvin.pvd", "compressed") << theta # save temperature to file

#--------------------------------------------------------------------------------------------------------------------------------------------------

