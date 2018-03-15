from fenics import *

################################
#### PROBLEM DEFINITION ########
################################

# Define geometry and mesh
d = 3
p0 = Point(0.0, 0.0, 0.0)
p1 = Point(0.30, 0.10, 0.10)
mesh = BoxMesh(p0, p1, 30, 10, 10)

# Define boundaries
boundaries = MeshFunction("size_t", mesh, d-1)
boundaries.set_all(0)
left, right, bottom, top = 1, 2, 3, 4
CompiledSubDomain("near(x[0], side) && on_boundary", side = p0[0]).mark(boundaries, left)
CompiledSubDomain("near(x[0], side) && on_boundary", side = p1[0]).mark(boundaries, right)
CompiledSubDomain("near(x[1], side) && on_boundary", side = p0[1]).mark(boundaries, bottom)
CompiledSubDomain("near(x[1], side) && on_boundary", side = p1[1]).mark(boundaries, top)

# Coordinate and surface integral element
x = SpatialCoordinate(mesh)
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

# Define function space
p = 2
V = FunctionSpace(mesh, "Lagrange", p)

# Define trial and test functions
theta = TrialFunction(V)
delta_theta = TestFunction(V)

# Material parameters
kappa = 1.0

# Volume force/ heat source and prescribed tractions/ prescribed heat fluxes
r = Expression("100000.0*exp(-((x[0]-0.1)*(x[0]-0.1)+(x[1]-0.033)*(x[1]-0.033))/0.0001)", degree=4)
q_p = Constant(1.0)

# Dirichlet boundary conditions
bcs = [DirichletBC(V, Constant((0.0)), boundaries, left),
	   #DirichletBC(V, Constant((-10.0)), boundaries, right)
	   ]

# Weak form a==l
a = kappa*dot(grad(delta_theta), grad(theta))*dx
l = r*delta_theta*dx + q_p*delta_theta*ds(right)

################################
#### ASSEMBLE AND SOLVE ########
################################


theta = Function(V)

A = assemble(a)
B = assemble(l)
for bc in bcs:
	bc.apply(A, B)
X = theta.vector()
solve(A, X, B)


################################
#### POST-PROCESSING ###########
################################

# Create temperature file
theta.rename("theta", "temperature")
file_theta = File("temperature.pvd", "compressed")
file_theta << theta

