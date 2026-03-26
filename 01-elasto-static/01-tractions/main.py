from fenics import *
#import  FEniCS Package 

"""-------------------------------------------------3D Elastostatics of a Rectangular Solid--------------------------------------------------
Problem description:

Geometry: Rectangular solid with dimensions 3.0 x 1.0 x 1.0
Boundary conditions:
   -Left face is fully clamped
   -Symmetry condition at x=1.5 (u_x = 0)
Loads:
   -Surface traction is applied on the top face
   -No body forces
Analysis type: Quasi-static model
Material model: Linear isotropic elasticity

Main Learnings:

Application of symmetry boundary conditions
Linear isotropic elasticity (Hooke's law)
Geometry and mesh setup, 
Boundary identification and marking
Function space definition (vector-valued FEM)
Variational (weak) formulation
Application of Dirichlet boundary conditions
Linear solver usage (MUMPS)
Stress recovery and output

--------------------------------------------------------------------------------------------------------------------------------------------------

"""

#---------------------------------------------------------------------------------------------------------
# Geometry and mesh generation
#---------------------------------------------------------------------------------------------------------


# spatial dimension
dim = 3 

# Define opposite corners of the rectangular domain
p0 = Point(0.0, 0.0, 0.0) #	bottom-left front corner
p1 = Point(1.5, 1.0, 1.0)  # top-right back corner

mesh = BoxMesh(p0, p1, 15, 10, 10)
# generates a structured 3D tetrahedral mesh within a box. 
# 30, 10 and 10 are the number of elements in x-, y- and z-directions


#---------------------------------------------------------------------------------------------------------
# Boundary identification and surface measures
#---------------------------------------------------------------------------------------------------------

# MeshFunction to store boundary IDs on facets
boundaries = MeshFunction("size_t", mesh, dim-1) 
# create mesh function for boundary domains, function stores unsigned integers "size_t" which act as boundary IDs
# "mesh"- links it to the mesh, "dim-1" indicates that the function is defined on the facets (2D entities) of the 3D mesh

boundaries.set_all(0) # Default all boundary IDs to 0 
left, right, bottom, top = 1, 2, 3, 4 # define boundary ID labels for easy reference

# filters every face based on its position and assigns the corresponding boundary ID
CompiledSubDomain("near(x[0], side) && on_boundary", side = p0[0]).mark(boundaries, left)   # near(x[0], side),&& on_boundary - checks if the x-coordinate of a point 
CompiledSubDomain("near(x[0], side) && on_boundary", side = p1[0]).mark(boundaries, right)  # on the boundary is near the specified side (p0[0] = 0.0 for left face) 
CompiledSubDomain("near(x[1], side) && on_boundary", side = p0[1]).mark(boundaries, bottom) # and on_boundary ensures that only boundary points are considered.
CompiledSubDomain("near(x[1], side) && on_boundary", side = p1[1]).mark(boundaries, top)    # If both conditions are met, that face is marked with the ID (eg'left' (1))

#  Surface integral element
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)
# 'ds' refers to the infinitesimal area elements on the exterior boundaries of the mesh
# 'domain=mesh' links this measure to our specific 3D box geometry.
# 'subdomain_data=boundaries' tells the measure to use our tagged IDs (1, 2, 3, 4).
# This allows us to integrate over specific faces, e.g., ds(top), ds(bottom).


#-------------------------------------------------------------------------------------------------------
# Function spaces and variational functions
#-------------------------------------------------------------------------------------------------------

# Polynomial degree of Lagrange finite elements
p = 2 
# adds "mid-side nodes" to every edge of every tetrahedron. 
# These mid-side nodes allow the polynomial basis functions to be quadratic (degree 2). 
# This enhances the solution's accuracy and smoothness.

# Define vector function space for displacement field.  
V = VectorFunctionSpace(mesh, "Lagrange", p)
# It creates three interconnected components (u_x, u_y, u_z) for every single point in the mesh.
# Lagrange elements ensure that the solution is tied to nodes and remains continuous across the entire 3D body, preventing physical gaps. 


# Define trial and test functions
# unknown displacement field
u = TrialFunction(V)
# The Trial Function represents the actual physical solution you are looking for 
# (in this case, the displacement of the beam)(Unknown)


# virtual displacement field
delta_u = TestFunction(V)
# delta_u represents an arbitrary virtual displacement.
# ensures the equilibrium of the system. It converts the differential equations 
# into the Weak Form that the computer can solve.



#--------------------------------------------------------------------------------------------------------
# Material properties
#--------------------------------------------------------------------------------------------------------

E = 200.e9  # Young's modulus in Pascals
rho = 8.e3  # density in kg/m^3
g = 9.81    # acceleration due to gravity in m/s^2
nu = 0.3    # Poisson's ratio
mu    = E/(2.0*(1.0 + nu)) # Lamé's second parameter (μ), commonly called the shear modulus; #Pascals
lmbda = E*nu/((1.0 + nu)*(1.0 - 2.0*nu)) # Lamé's first parameter (λ); #Pascals (see Lamé parameters — Wikipedia https://en.wikipedia.org/wiki/Lam%C3%A9_parameters)


#--------------------------------------------------------------------------------------------------------
# Loads and boundary conditions
#--------------------------------------------------------------------------------------------------------

# Volume force (N/m^3)
b = Constant((0.0, 0.0, 0.0)) 
# force that acts on every particle within the material. (eg gravity).

#prescribed tractions in Pascals
t_p = Constant((0.0, -1.e6, 0.0))  
# force applied on the surface of the material (no traction force applied here)


# Define Dirichlet boundary conditions
bcs = [DirichletBC(V, Constant((0.0, 0.0, 0.0)), boundaries, left), # fixed left edge 
       DirichletBC(V.sub(0), Constant(0.0), boundaries, right),  # symmetry condition on right edge (u_x = 0)
       ]
# DirichletBC - used to impose fixed displacements on specific boundaries of the mesh.


#---------------------------------------------------------------------------------------------------------
#  Variational formulation (weak form of linear elasticity)
#---------------------------------------------------------------------------------------------------------

# Strain tensor
def epsilon(u):
    return sym(grad(u))

# Stress tensor for isotropic linear elasticity
def sigma(u):
    return lmbda*tr(epsilon(u))*Identity(dim) + 2.0*mu*epsilon(u) # stress tensor (see Stress tensor — Wikipedia https://en.wikipedia.org/wiki/Lam%C3%A9_parameters)


# Bilinear form:  internal virtual work
a = inner(grad(delta_u), sigma(u))*dx
# left-hand side of the variational formulation (see Pg 144/148 LKM Slides) 

# linear form: external virtual work 
l = dot(b, delta_u)*dx + dot(t_p, delta_u)*ds(top)
# right-hand side of the variational formualation 

# Explanation of terms:
# grad - gradient operator, measures how a field changes in space
# tr - trace of a tensor, sum of diagonal elements
# Identity(d) - identity tensor of dimension d (3 here)
# inner - double contraction of two tensors
# dx - volume integral over the entire 3D domain
# ds(top) - surface integral over the top face only
# dot - standard dot product between vectors
# sym - symmetric part of a tensor


#---------------------------------------------------------------------------------------------------------
# Solve the linear system
#---------------------------------------------------------------------------------------------------------

# Solution function to store displacement field
u = Function(V) 
# It creates an array of zeros with the exact size needed to store the 
# results (u_x, u_y, u_z) for every node in mesh. 
# After the solve command runs, this u will hold the final answer.


# Solve Ku = f using the direct MUMPS solver
solve(a == l, u, bcs=bcs, 
	      solver_parameters={"linear_solver": "mumps"},
	      form_compiler_parameters={"optimize": True}) # solve the variational problem

# Explanation of terms:
# a == l: This tells the solver to find a u such that the internal work (a) perfectly balances the external work (l).
# bcs=bcs: This injects your Dirichlet conditions into the system of equations (KU = F), removing the "unknowns" at those specific nodes.
# MUMPS is particularly efficient for large, sparse systems like those arising in finite element analyses.
# form_compiler_parameters={"optimize": True}: This enables optimizations during the compilation of the variational forms, which can speed up the assembly and solution process.

#---------------------------------------------------------------------------------------------------------
# Post processing
#---------------------------------------------------------------------------------------------------------

# Save displacement
u.rename("u", "displacement") # rename displacement for output
File("displacement_in_meters.pvd", "compressed") << u # save displacement to file
# displacement.pvd is the output file that can be visualized in ParaView or other compatible visualization software.

# Compute stress
T = TensorFunctionSpace(mesh, "Lagrange", p) 
# To save the stress, new Function Space that can hold tensors instead of just vectors is used.

# Calculates the stress based on your solved displacements.
stress = project(sigma(u)/1.e6, T, solver_type="mumps") # stress in MPa
# Converts the units from Pascals to MegaPascals (MPa)
# Stress is physically calculated at internal quadrature (integration) points. 
# We 'project' it to the nodes of space T to create a smooth, continuous field for visualization in Paraview.

stress.rename("sigma", "stress") # rename stress for output
File("stress_in_MPa.pvd", "compressed") << stress # save stress to file
# stress_in_MPa.pvd is the output file that can be visualized in ParaView or other compatible visualization software.


#------------------------------------------------------------------------------------------------------------------------------------------
