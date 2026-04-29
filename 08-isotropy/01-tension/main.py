from fenics import *

"""----------------------------------------------------3D Static Linear Elastic Analysis of a Unit Cube-----------------------------------------------------------------------------------------

Problem Overview:
This script performs a 3D steady-state structural analysis of a cube under prescribed displacement using FEniCS. It calculates the displacement field 
and the resulting internal stress distribution.

Geometry:
- 3D Cube: 1.0 m x 1.0 m x 1.0 m (defined by Point(0,0,0) and Point(1,1,1))
- Structured Mesh: 3 x 3 x 3 hexahedral elements (BoxMesh)

Material Model:
- Linear Isotropic Elasticity (Hooke's Law)
- Young's Modulus (E): 200 GPa
- Poisson's Ratio (nu): 0.3
- Lame Parameters: Derived for 3D/Plane Strain conditions

Boundary Conditions:
- Dirichlet BCs: 
    * Left face (x=0): Fixed x-displacement
    * Bottom face (y=0): Fixed y-displacement
    * Back face (z=0): Fixed z-displacement
    * Right face (x=1): Prescribed tensile displacement of 0.2 m
- Neumann BCs: Zero traction (t_p) on all other surfaces
-----------------------------------------------------------------------------------------------------------------------------------------"""


#---------------------------------------------------------------------------------------------------------
# Geometry and mesh generation
#---------------------------------------------------------------------------------------------------------

# Define geometry and mesh
dim = 3 # spatial dimension
p0 = Point(0.0, 0.0, 0.0)  # front lower left corner 
p1 = Point(1.0, 1.0, 1.0)  # back upper right corner 
mesh = BoxMesh(p0, p1, 3, 3, 3) # apply a structured mesh with 3 elements in each direction

#---------------------------------------------------------------------------------------------------------
# Boundary identification and marking
#---------------------------------------------------------------------------------------------------------

# Define boundaries
boundaries = MeshFunction("size_t", mesh, dim-1) # MeshFunction to store boundary IDs on facets
boundaries.set_all(0) # initialize all boundaries to 0
left, right, bottom, top, back, front = 1, 2, 3, 4, 5, 6 # boundary IDs
CompiledSubDomain("near(x[0], side) && on_boundary", side = p0[0]).mark(boundaries, left) # mark left boundary with ID 1
CompiledSubDomain("near(x[0], side) && on_boundary", side = p1[0]).mark(boundaries, right) # mark right boundary with ID 2
CompiledSubDomain("near(x[1], side) && on_boundary", side = p0[1]).mark(boundaries, bottom) # mark bottom boundary with ID 3
CompiledSubDomain("near(x[1], side) && on_boundary", side = p1[1]).mark(boundaries, top) # mark top boundary with ID 4
CompiledSubDomain("near(x[2], side) && on_boundary", side = p0[2]).mark(boundaries, back) # mark back boundary with ID 5
CompiledSubDomain("near(x[2], side) && on_boundary", side = p1[2]).mark(boundaries, front) # mark front boundary with ID 6

# Surface integral element
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

#-------------------------------------------------------------------------------------------------------
# Function spaces and variational functions
#-------------------------------------------------------------------------------------------------------

# Define function space
p = 2 # polynomial degree for finite element 
V = VectorFunctionSpace(mesh, "Lagrange", p) # vector function space for displacement field

# Define trial and test functions
u = TrialFunction(V) # trial function for unknown displacement field
delta_u = TestFunction(V) # test function for virtual displacement field

#--------------------------------------------------------------------------------------------------------
# Material properties
#--------------------------------------------------------------------------------------------------------

E = 200.e9 # Young's modulus in Pa
nu = 0.3 # Poisson's ratio
mu    = E/(2.0*(1.0 + nu)) # shear modulus # Lame's second parameter
lmbda = E*nu/((1.0 + nu)*(1.0 - 2.0*nu)) # Lame's first parameter
#mu    = E/(2.0*(1.0 + nu)) #PLANE STRESS
#lmbda = E*nu/((1.0 + nu)*(1.0 - nu)) #PLANE STRESS

# Volume force/ heat source and prescribed tractions/ prescribed heat fluxes
b = Constant((0.0, 0.0, 0.0)) # body force in N/m^3
t_p = Constant((0.0, 0.0, 0.0)) # surface traction vector  in N/m^2 (Pa) - positive for tension case, negative for compression case

# Dirichlet boundary conditions
#class Near(SubDomain):
#	def __init__(self, point):
#		self.p = point
#		SubDomain.__init__(self)
#	def inside(self, x, on_boundary):
#		return (abs(x[0]-self.p[0]) < DOLFIN_EPS and \
#						       abs(x[1]-self.p[1]) < DOLFIN_EPS and \
#						       abs(x[2]-self.p[2]) < DOLFIN_EPS)

# Define a subdomain class representing the origin point (0,0,0)
#class Origin(SubDomain):
#	def inside(self, x, on_boundary):
		# This function checks whether a given point x lies at the origin
        # DOLFIN_EPS is a small tolerance used to avoid floating-point errors
#		return (abs(x[0]) < DOLFIN_EPS) and (abs(x[1]) < DOLFIN_EPS) and (abs(x[2]) < DOLFIN_EPS)

# Define a subdomain class representing the Z-axis (x = 0 and y = 0)
# class AxisZ(SubDomain):
#	def inside(self, x, on_boundary):
		# Checks if a point lies on the Z-axis
        # i.e., x-coordinate = 0 and y-coordinate = 0 (z can be anything)
#		return (abs(x[0]) < DOLFIN_EPS) and (abs(x[1]) < DOLFIN_EPS)

# Define Dirichlet boundary conditions (prescribed displacement)	       
bcs = [DirichletBC(V.sub(0), Constant(0.0), boundaries, left), # fix x-displacement on left boundary
	   DirichletBC(V.sub(1), Constant(0.0), boundaries, bottom), # fix y-displacement on bottom boundary
	   DirichletBC(V.sub(2), Constant(0.0), boundaries, back), # fix z-displacement on back boundary
	   DirichletBC(V.sub(0), Constant(0.2), boundaries, right) # prescribe x-displacement of 0.2 m on right boundary
       ]

#---------------------------------------------------------------------------------------------------------
#  Variational formulation of the balance of linear momentum (weak form)
#---------------------------------------------------------------------------------------------------------
# Strain tensor
def epsilon(u):
	return sym(grad(u))

# Stress tensor (linear isotropic elasticity)
def sigma(u):
    return lmbda*tr(epsilon(u))*Identity(dim) + 2.0*mu*epsilon(u)
		   
# Weak form a==l
a = inner(grad(delta_u), sigma(u))*dx # bilinear form
l = dot(b, delta_u)*dx + dot(t_p, delta_u)*ds(right) # linear form

#---------------------------------------------------------------------------------------------------------
# Solve the linear system
#---------------------------------------------------------------------------------------------------------

# Solution function to store displacement field
u = Function(V)

# Solve Ku = f using the direct MUMPS solver
solve(a == l, u, bcs=bcs,
      solver_parameters={"linear_solver": "mumps"},
      form_compiler_parameters={"optimize": True})

#---------------------------------------------------------------------------------------------------------
# Post processing
#---------------------------------------------------------------------------------------------------------

# Create displacement file
u.rename("u", "displacement") # rename displacement for output
File("displacement_in_meters.pvd", "compressed") << u # save displacement to file

# Project stress field and create stress file
T = TensorFunctionSpace(mesh, "Lagrange", p) # tensor function space for stress field
stress = project(sigma(u)/1.e6, T, solver_type="mumps") # project stress in MPa for output
stress.rename("sigma", "stress") #rename stress for output

File("stress_in_MPa.pvd", "compressed") << stress # save stress to file


#---------------------------------------------------------------------------------------------------------
