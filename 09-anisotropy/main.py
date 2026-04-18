from sympy import Identity

from fenics import *

"""
----------------------------------Static Analysis of a Transversely Isotropic 3D Cube----------------------------------

Problem Overview:
This script performs a 3D steady-state structural analysis of a unit cube featuring 
transversely isotropic material properties. It specifically examines the effect of 
fiber orientation on the mechanical response under prescribed tensile displacement.

Geometry & Mesh:
- 3D Unit Cube: 1.0 m x 1.0 m x 1.0 m
- Structured Mesh: 5 x 5 x 5  elements (BoxMesh)

Material Model:
- Transversely Isotropic Linear Elasticity:
    * Longitudinal Young's Modulus (E_L): 44 GPa (Fiber direction)
    * Transverse Young's Modulus (E_T): 13 GPa
    * Fiber Direction (a1): Oriented along the Y-axis (0.0, 1.0, 0.0)

Boundary Conditions:
- Dirichlet BCs (Mixed Support):
    * Left face (x=0): Fixed x-displacement
    * Z-Axis (x=0, y=0): Fixed y-displacement using pointwise method to prevent translation
    * Origin (0,0,0): Fixed z-displacement using pointwise method to eliminate rigid body motion
    * Right face (x=1): Prescribed tensile displacement of 0.02 m (2 cm)
- Neumann BCs: Zero traction (t_p) on all other surfaces

Main Learnings:
- Implementing invariant-based constitutive laws for anisotropic materials in FEniCS
- Using 'SubDomain' classes (Origin, AxisZ) for precise pointwise rigid body constraints
- Defining structural tensors (A1 = a1 ⊗ a1) to handle directional material stiffness

---------------------------------------------------------------------------------------------------------------------------"""

#---------------------------------------------------------------------------------------------------------
# Geometry and mesh generation
#---------------------------------------------------------------------------------------------------------


# Define geometry and mesh
dim = 3 # spatial dimension
p0 = Point(0.0, 0.0, 0.0)  # front lower left corner 
p1 = Point(1.0, 1.0, 1.0)  # back upper right corner 
mesh = BoxMesh(p0, p1, 5, 5, 5) # apply a structured mesh with 3 elements in each direction


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

# Material parameters
E_L = 44000*1.e6  #=E_11  # Young's modulus in longitudinal direction in Pa
E_T = 13000*1.e6  #=E_22=E_33  # Young's modulus in transverse directions in Pa
nu_LT = 0.25 #=nu=nu_12=nu_13 # Poisson's ratio for strain in longitudinal direction when stretched in transverse direction
nu_TT = 0.3  #=nu_23=nu_32 # Poisson's ratio in transverse plane 
mu_LT = 5600*1.e6 #=G_12=G_13 # Longitudinal shear modulus in Pa
mu_TT = E_T/(2*(1+nu_TT)) #=G_23 # Transverse shear modulus in Pa
nu_TL = E_T/E_L*nu_LT # Poisson's ratio for strain in transverse direction when stretched in longitudinal direction

lmbda = (nu_LT*nu_TL+nu_TT)/(1-nu_TT-2*nu_LT*nu_TL)/(1+nu_TT)*E_T # Modified Lamé parameter for anisotropic material # Used in the stress-strain relationship
alpha = (nu_LT*(1+nu_TT-nu_TL)-nu_TT)/(1-nu_TT-2*nu_LT*nu_TL)/(1+nu_TT)*E_T # Coupling parameter controlling interaction between longitudinal and transverse strains
beta = E_L*(1-nu_TT*nu_TT)/(1-nu_TT-2*nu_LT*nu_TL)/(1+nu_TT)-lmbda-4*mu_LT # Additional stiffness parameter required for transversely isotropic elasticity

#a1 = Constant((1.0, 0.0, 0.0))
a1 = Constant((0.0, 1.0, 0.0)) # Fiber direction along the y-axis (longitudinal direction)

# alternate fiber direction along 45 degree in x-y plane
#a1 = Constant((1.0/sqrt(2), 1.0/sqrt(2), 0.0)) # direction of anisotropy along 45 degree in x-y plane # This means the material is stiffest along the 45 degree direction in the x-y plane

A1 = outer(a1, a1)
# Creates a second-order tensor representing the fiber orientation
# A1 = a1 ⊗ a1
# This tensor is used in anisotropic stress calculations

#rho = 8.e3
#g = 9.81


# Volume force/ heat source and prescribed tractions/ prescribed heat fluxes
b = Constant((0.0, 0.0, 0.0)) # body force in N/m^3
t_p = Constant((0.0, 0.0, 0.0)) # surface traction vector  in N/m^2 (Pa) 

# Dirichlet boundary conditions
#class Near(SubDomain):
#	def __init__(self, point):
#		self.p = point
#		SubDomain.__init__(self)
#	def inside(self, x, on_boundary):
#		return (abs(x[0]-self.p[0]) < DOLFIN_EPS and \
#						       abs(x[1]-self.p[1]) < DOLFIN_EPS and \
#						       abs(x[2]-self.p[2]) < DOLFIN_EPS)

class Origin(SubDomain):
	def inside(self, x, on_boundary):
		#checks whether a given point x lies at the origin
		return (abs(x[0]) < DOLFIN_EPS) and (abs(x[1]) < DOLFIN_EPS) and (abs(x[2]) < DOLFIN_EPS) 
	# DOLFIN_EPS is a small, predefined constant used as a numerical tolerance. Its value is typically 1.11 x 10^(-16)

class AxisZ(SubDomain):
	def inside(self, x, on_boundary): 
		# Checks if a point lies on the Z-axis
		return (abs(x[0]) < DOLFIN_EPS) and (abs(x[1]) < DOLFIN_EPS)

# Define Dirichlet boundary conditions (prescribed displacement)			       
bcs = [DirichletBC(V.sub(0), Constant(0.0), boundaries, left), # fix x-displacement on left boundary
	   DirichletBC(V.sub(1), Constant(0.0), AxisZ(), method="pointwise"), # fix y-displacement along Z-axis (x=0 and y=0)
	   DirichletBC(V.sub(2), Constant(0.0), Origin(), method="pointwise"), # fix z-displacement at the origin point (0,0,0)
	   DirichletBC(V.sub(0), Constant(0.02), boundaries, right) # prescribe x-displacement of 0.02 m on right boundary
       ]
#---------------------------------------------------------------------------------------------------------
#  Variational formulationn (weak form)
#---------------------------------------------------------------------------------------------------------

# Strain tensor
def epsilon(u):
	return sym(grad(u))

# Stress tensor # transversely isotropic elasticity
def sigma(u):
    return lmbda*tr(epsilon(u))*Identity(dim) + 2.0*mu_TT*epsilon(u) \
    	   + alpha*(dot(a1, dot(epsilon(u), a1))*Identity(dim) + tr(epsilon(u))*A1) \
		   + 2.0*(mu_LT-mu_TT)*(dot(A1, epsilon(u))+dot(epsilon(u), A1)) \
		   + beta*dot(a1, dot(epsilon(u), a1))*A1
		   
# Material parameters
#E = 200.e9
#nu = 0.3
#mu    = E/(2.0*(1.0 + nu))
#lmbda = E*nu/((1.0 + nu)*(1.0 - 2.0*nu))

# Stress tensor (linear isotropic elasticity)
#def sigma(u):
#    return lmbda*tr(sym(grad(u)))*Identity(d) + 2.0*mu*sym(grad(u))
		   
# Weak form a==l
a = inner(grad(delta_u), sigma(u))*dx # Bilinear form (internal virtual work)
l = dot(b, delta_u)*dx + dot(t_p, delta_u)*ds(right) # Linear form (external virtual work) - body forces and surface tractions on the right boundary

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

# Create tensor function space for stress field
T = TensorFunctionSpace(mesh, "Lagrange", p) 

# project stress in MPa for output using an iterative solver with AMG preconditioner for better performance
stress = project(sigma(u)/1.e6, T, solver_type="cg", preconditioner_type="petsc_amg")
stress.rename("sigma", "stress") # rename stress for output

File("stress_in_MPa.pvd", "compressed") << stress # save stress to file

#---------------------------------------------------------------------------------------------------------