from fenics import *

"""-------------------------------------------------3D Elastostatics of a Rectangular Solid--------------------------------------------------
Problem description:
Geometry:Rectangular solid with dimensions 3.0 x 1.0 x 1.0
Boundary conditions:
    -Left face is fully clamped
Loads:
    -No surface traction is applied
    -body forces are applied (gravity load in negative y direction)
Analysis type: Quasi-static model
Material model: Linear isotropic elasticity

Main learnings:
GRAVITY LOADS - Implementation of body forces in the variational formulation
--------------------------------------------------------------------------------------------------------------------------------------------------
"""

#---------------------------------------------------------------------------------------------------------
# Geometry and mesh generation
#---------------------------------------------------------------------------------------------------------

# spatial dimension
dim = 3 

# Define opposite corners of the rectangular domain
p0 = Point(0.0, 0.0, 0.0)  # bottom-left front corner
p1 = Point(3.0, 1.0, 1.0)  # top-right back corner
mesh = BoxMesh(p0, p1, 30, 10, 10) # Create a structured mesh of the rectangular solid

#---------------------------------------------------------------------------------------------------------
# Boundary identification and marking
#---------------------------------------------------------------------------------------------------------

# MeshFunction to store boundary IDs on facets
boundaries = MeshFunction("size_t", mesh, dim-1)
boundaries.set_all(0) # set every face to the value 0
left, right, bottom, top = 1, 2, 3, 4 # define boundary IDs

# filters every face based on its position and assigns the corresponding boundary ID
CompiledSubDomain("near(x[0], side) && on_boundary", side = p0[0]).mark(boundaries, left)   
CompiledSubDomain("near(x[0], side) && on_boundary", side = p1[0]).mark(boundaries, right)  
CompiledSubDomain("near(x[1], side) && on_boundary", side = p0[1]).mark(boundaries, bottom) 
CompiledSubDomain("near(x[1], side) && on_boundary", side = p1[1]).mark(boundaries, top)    

#  Surface integral element
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)


#-------------------------------------------------------------------------------------------------------
# Function spaces and variational functions
#-------------------------------------------------------------------------------------------------------

# Polynomial degree of Lagrange finite elements
p = 2 


# Define vector function space for displacement field.  
V = VectorFunctionSpace(mesh, "Lagrange", p)

# Define trial and test functions

# unknown displacement field to be solved for
u = TrialFunction(V)

# virtual displacement field for the variational formulation
delta_u = TestFunction(V)


#--------------------------------------------------------------------------------------------------------
# Material properties
#--------------------------------------------------------------------------------------------------------

E = 200.e9  # Young's modulus in Pascals
rho = 8.e3  # density in kg/m^3
g = 9.81    # acceleration due to gravity in m/s^2
nu = 0.3    # Poisson's ratio
mu    = E/(2.0*(1.0 + nu)) # Lamé's second parameter (μ), commonly called the shear modulus; # Pascals
lmbda = E*nu/((1.0 + nu)*(1.0 - 2.0*nu)) # Lamé's first parameter (λ); # Pascals (see Lamé parameters — Wikipedia https://en.wikipedia.org/wiki/Lam%C3%A9_parameters)


#--------------------------------------------------------------------------------------------------------
# Loads and boundary conditions
#--------------------------------------------------------------------------------------------------------

# Volume force (N/m^3)
b = Constant((0.0, -rho*g, 0.0))

# Prescribed tractions # Pascals
t_p = Constant((0.0, 0.0, 0.0))  

# Define Dirichlet boundary conditions
bcs = [DirichletBC(V, Constant((0.0, 0.0, 0.0)), boundaries, left), # fixed left edge 
       ]

#---------------------------------------------------------------------------------------------------------
#  Variational formulation (weak form of linear elasticity)
#---------------------------------------------------------------------------------------------------------
# Strain tensor
def epsilon(u):
    return sym(grad(u))

# Stress tensor for isotropic linear elasticity
def sigma(u):
    return lmbda*tr(epsilon(u))*Identity(dim) + 2.0*mu*epsilon(u) # stress tensor (see Stress tensor — Wikipedia https://en.wikipedia.org/wiki/Lam%C3%A9_parameters)


# Bilinear form: internal virtual work
a = inner(grad(delta_u), sigma(u))*dx

# Linear form: external virtual work
l = dot(b, delta_u)*dx + dot(t_p, delta_u)*ds(top)


#---------------------------------------------------------------------------------------------------------
# Solve the linear system
#---------------------------------------------------------------------------------------------------------

# Solution function to store displacement field
u = Function(V) 


# Solve K u = f using the direct MUMPS solver
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
