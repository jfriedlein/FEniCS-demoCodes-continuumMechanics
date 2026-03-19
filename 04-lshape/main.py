from fenics import *

"""-------------------------------------------------2D Elastostatics of a L shaped Solid--------------------------------------------------
Problem description: Stress singularity at the re-entrant corner

Geometry: L shaped solid

Boundary conditions:
	-Top face is fully clamped
      
Loads:
	- surface traction is applied in 45 degree clockwise direction on the right face
	- No body forces applied
      
Analysis type: Quasi-static elastic model, plane stress condition
Material model: Linear Isotropic elastic material

Main learnings:

sharp corners cause singularities
singularities cause infinite stress (cannot be evaluated or interpreted)

--------------------------------------------------------------------------------------------------------------------------------------------------

"""
#---------------------------------------------------------------------------------------------------------
# Geometry and mesh generation
#---------------------------------------------------------------------------------------------------------

# spatial dimension
dim = 2

# load mesh
mesh = Mesh("lshape.xml.gz")
# load pre-defined mesh of L-shape domain

mesh = refine(mesh)
# Refines all cells uniformly:
# Each triangle (2D) - 4 smaller triangles
# Each tetrahedron (3D) - 8 smaller tetrahedra
# Improves resolution everywhere.

margin = 0.5 # This defines a horizontal band around the middle of the domain

# refine mesh locally around the re-entrant corner
for i in range (0, 4):   # refine mesh 4 times around the re-entrant corner
	margin = margin*2.0/3.0 # reduce margin size
	markers = MeshFunction("bool", mesh, dim) # cell markers for refinement - true/false
	markers.set_all(False) # initialize all cells to false
	CompiledSubDomain("x[1]>(0.5-m) && x[1]<(0.5+m)", m=margin).mark(markers, True) # mark cells within the margin band as true
	mesh = refine(mesh, markers) # refine only the marked cells

#for i in range (0, 2):
#	margin = 0.1
#	markers = CellFunction("bool", mesh)
#	markers.set_all(False)
#	CompiledSubDomain("x[0]>(0.5-m) && x[0]<(0.5+m) && x[1]>(0.5-m) && x[1]<(0.5+m)", m=margin).mark(markers, True)
#	mesh = refine(mesh, markers)

#---------------------------------------------------------------------------------------------------------
# Boundary identification and measures
#---------------------------------------------------------------------------------------------------------

# MeshFunction to store boundary IDs on facets
boundaries = MeshFunction("size_t", mesh, dim-1)
boundaries.set_all(0) # set every face to the value 0
left, right, bottom, top, = 1, 2, 3, 4 # define boundary IDs
CompiledSubDomain("near(x[0], side) && on_boundary", side = 1.0).mark(boundaries, left) # left face
CompiledSubDomain("near(x[0], side) && on_boundary", side = 1.0).mark(boundaries, right) # right face
CompiledSubDomain("near(x[1], side) && on_boundary", side = 1.0).mark(boundaries, bottom) #	 bottom face
CompiledSubDomain("near(x[1], side) && on_boundary", side = 1.0).mark(boundaries, top) # top face

# Surface integral element
ds = Measure('ds', domain=mesh, subdomain_data=boundaries) 


#-------------------------------------------------------------------------------------------------------
# Function spaces and variational functions
#-------------------------------------------------------------------------------------------------------

# Polynomial degree of Lagrange finite elements
p = 2

# Define function space for displacement field
V = VectorFunctionSpace(mesh, "Lagrange", p) 

# Define trial and test functions

# unknown displacement field to be solved for
u = TrialFunction(V)

# virtual displacement field for the variational formulation
delta_u = TestFunction(V)

#--------------------------------------------------------------------------------------------------------
# Material properties
#--------------------------------------------------------------------------------------------------------

E = 200.e9 # Young's modulus in Pascals
nu = 0.3 # Poisson's ratio
#mu    = E/(2.0*(1.0 + nu))
#lmbda = E*nu/((1.0 + nu)*(1.0 - 2.0*nu))
mu    = E/(2.0*(1.0 + nu)) #PLANE STRESS # shear modulus # lame second parameter # Pascals
lmbda = E*nu/((1.0 + nu)*(1.0 - nu)) #PLANE STRESS # lame first parameter #Pascals

#rho = 8.e3
#g = 9.81

#--------------------------------------------------------------------------------------------------------
# Loads and boundary conditions
#--------------------------------------------------------------------------------------------------------

# Volume force in N/m^3
b = Constant((0.0, 0.0))

# prescribed tractions in Pascals
t_p = Constant((1.e6, -1.e6)) 

# Dirichlet boundary conditions 
bcs = [DirichletBC(V, Constant((0.0, 0.0)), boundaries, top) # fixed top edge
       ]

#---------------------------------------------------------------------------------------------------------
#  Variational formulation (weak form of linear elasticity)
#---------------------------------------------------------------------------------------------------------

# Strain tensor
def epsilon(u):
	return sym(grad(u))

# Stress tensor (linear isotropic elasticity)
def sigma(u):
    return lmbda*tr(epsilon(u))*Identity(dim) + 2.0*mu*epsilon(u)

# Bilinear form:  internal virtual work
a = inner(grad(delta_u), sigma(u))*dx

# linear form: external virtual work
l = dot(b, delta_u)*dx + dot(t_p, delta_u)*ds(right)

# Potential energy (M) = Strain energy - Work done by external forces
M = 0.5*inner(grad(u), sigma(u))*dx - dot(b, u)*dx - dot(t_p, u)*ds(right)


#---------------------------------------------------------------------------------------------------------
# Solve the linear system
#---------------------------------------------------------------------------------------------------------

# Solution function to store displacement field
u = Function(V) 


# Solve Ku = f using the direct MUMPS solver
solve(a == l, u, bcs=bcs, 
	      solver_parameters={"linear_solver": "mumps"},
	      form_compiler_parameters={"optimize": True}) # solve the variational problem

#---------------------------------------------------------------------------------------------------------
#  Post processing:
#---------------------------------------------------------------------------------------------------------

# Save displacement
u.rename("u", "displacement") # rename displacement for output
File("displacement_in_meters.pvd", "compressed") << u # save displacement to file

# Compute stress
T = TensorFunctionSpace(mesh, "Lagrange", p) 

# Calculates the stress based on your solved displacements.
stress = project(sigma(u)/1.e6, T, solver_type="mumps") # stress in MPa

# saves the stress to file
stress.rename("sigma", "stress") # rename stress for output
File("stress_in_MPa.pvd", "compressed") << stress 

#------------------------------------------------------------------------------------------------------------------------------------------


