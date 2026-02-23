from fenics import *


"""-------------------------------------------------3D Thermostatics of a Rectangular Solid--------------------------------------------------
Problem description: prescribed heat flux and heat source

Geometry: Rectangular solid with dimensions 0.30 x 0.10 x 0.10
Boundary conditions:
	-Left face is held at 0 degree Celsius
	
Loads:
	- Prescribed heat flux is applied on the right face
	- A volumetric heat source is applied over the domain

Analysis type: Quasi-static model
Material model: Linear Isotropic heat conduction


Main learnings:

Weak Formulation of Heat Conduction 
Spatial Coordinate - Usage of spatial coordinates for defining expressions over the domain
--------------------------------------------------------------------------------------------------------------------------------------------------

"""

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
left, right, bottom, top = 1, 2, 3, 4 # define boundary IDs

# filters every face based on its position and assigns the corresponding boundary ID
CompiledSubDomain("near(x[0], side) && on_boundary", side = p0[0]).mark(boundaries, left)   
CompiledSubDomain("near(x[0], side) && on_boundary", side = p1[0]).mark(boundaries, right)  
CompiledSubDomain("near(x[1], side) && on_boundary", side = p0[1]).mark(boundaries, bottom) 
CompiledSubDomain("near(x[1], side) && on_boundary", side = p1[1]).mark(boundaries, top)    

# Coordinate and surface integral element
x = SpatialCoordinate(mesh) # spatial coordinates for defining expressions over the domain - x[1], x[2], x[3]

#  Surface integral element
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

#-------------------------------------------------------------------------------------------------------
# Function spaces and variational functions
#-------------------------------------------------------------------------------------------------------

# Polynomial degree of Lagrange finite elements
p = 2

# Define function space for temperature field.
V = FunctionSpace(mesh, "Lagrange", p)

# Define trial and test functions
theta = TrialFunction(V) # unknown temperature field to be solved for
delta_theta = TestFunction(V) # virtual temperature field for the variational formulation



#--------------------------------------------------------------------------------------------------------
# Material properties
#--------------------------------------------------------------------------------------------------------

# thermal conductivity in W/mK
kappa = 1.0 



#--------------------------------------------------------------------------------------------------------
# Loads and boundary conditions
#--------------------------------------------------------------------------------------------------------

# heat source term
r = Expression("100000.0*exp(-((x[0]-0.1)*(x[0]-0.1)+(x[1]-0.033)*(x[1]-0.033))/0.0001)", degree=4)
# r is a localized heat source
# Has a maximum value of 100,000 at x0i = (0.1, 0.033)
# rapidly decays away from that point
# degree=4 controls the interpolation accuracy of the Expression.

# prescribed heat flux in W/m2
q_p = Constant(1.0)


# Dirichlet boundary conditions
bcs = [DirichletBC(V, Constant((0.0)), boundaries, left),    # left face at 0.0 temperature
	   ]
#---------------------------------------------------------------------------------------------------------
#  Variational formulation (weak form of heat conduction)
#---------------------------------------------------------------------------------------------------------

# bilinear form: internal heat conduction
# left-hand side of the variational formulation
a = kappa*dot(grad(delta_theta), grad(theta))*dx  #see Pg 148 LKM Slides

# linear form: volumetric heat source+ boundary heat flux 
# right-hand side of the variational formulation
l = r*delta_theta*dx + q_p*delta_theta*ds(right)


#---------------------------------------------------------------------------------------------------------
# Solve the linear system
#---------------------------------------------------------------------------------------------------------

# Solution function to store temperature field
theta = Function(V)

# Solve Ku = f using the direct MUMPS solver
solve(a == l, theta, bcs=bcs, 
	      solver_parameters={"linear_solver": "mumps"},
	      form_compiler_parameters={"optimize": True}) # solve the variational problem

#---------------------------------------------------------------------------------------------------------
#  Post processing:
#---------------------------------------------------------------------------------------------------------

# Save temperature
theta.rename("theta", "temperature") # rename temperature for output
File("temperature_in_degrees.pvd", "compressed") << theta # save temperature to file


#---------------------------------------------------------------------------------------------------------