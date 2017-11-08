from __future__ import print_function
from dolfin import *

parameters['form_compiler']['optimize'] = True
parameters['form_compiler']['cpp_optimize'] = True

################################
#### PROBLEM DEFINITION ########
################################

#list_linear_solver_methods()

print("Preprocessing... ")

# Define geometry and mesh
d = 3
p0 = Point(0.0, 0.0, 0.0)
p1 = Point(0.03, 0.002, 0.005)
mesh = BoxMesh(p0, p1, 50, 10, 5)

# Define boundaries
boundaries = FacetFunction("size_t", mesh)
boundaries.set_all(0)
left, right, bottom, top, back, front = 1, 2, 3, 4, 5, 6
CompiledSubDomain("near(x[0], side) && on_boundary", side = p0[0]).mark(boundaries, left)
CompiledSubDomain("near(x[0], side) && on_boundary", side = p1[0]).mark(boundaries, right)
CompiledSubDomain("near(x[1], side) && on_boundary", side = p0[1]).mark(boundaries, bottom)
CompiledSubDomain("near(x[1], side) && on_boundary", side = p1[1]).mark(boundaries, top)
CompiledSubDomain("near(x[2], side) && on_boundary", side = p0[2]).mark(boundaries, back)
CompiledSubDomain("near(x[2], side) && on_boundary", side = p1[2]).mark(boundaries, front)

# Coordinate and surface integral element
x = SpatialCoordinate(mesh)
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

# Define function space
p = 2
P = VectorElement("Lagrange", mesh.ufl_cell(), p, 3)
Q = FiniteElement("Lagrange", mesh.ufl_cell(), p)
V = FunctionSpace(mesh, P*Q)

# Define trial and test functions
(u, theta) = TrialFunctions(V)
(delta_u, delta_theta) = TestFunctions(V)

# Material parameters
E = 200.e9
rho = 8.e3
g = 9.81
nu = 0.3
mu    = E/(2.0*(1.0 + nu))
lmbda = E*nu/((1.0 + nu)*(1.0 - 2.0*nu))
kappa = 80.0
alpha = 1.0*conditional(lt(x[1], p1[1]/2.0), 1.e-6, 15.0e-6)

# Volume force/ heat source and prescribed tractions/ prescribed heat fluxes
b = Constant((0.0, 0.0, 0.0))
t_p = Constant((0.0, 0.0, 0.0))
r = Constant(0.0)
q_p = Constant(0.0)

# Dirichlet boundary conditions
bcs = [DirichletBC(V.sub(0), Constant((0.0, 0.0, 0.0)), boundaries, left),
       DirichletBC(V.sub(1), Constant(25.0), boundaries, left),
       DirichletBC(V.sub(1), Constant(25.0), boundaries, right),
       DirichletBC(V.sub(1), Constant(25.0), boundaries, bottom),
       DirichletBC(V.sub(1), Constant(25.0), boundaries, top),
       DirichletBC(V.sub(1), Constant(25.0), boundaries, back),
       DirichletBC(V.sub(1), Constant(25.0), boundaries, front)
       ]

# Stress tensor (linear isotropic elasticity)
def sigma(u, theta):
    return lmbda*tr(sym(grad(u)))*Identity(3) + 2.0*mu*sym(grad(u)) - alpha*(3.0*lmbda + 2.0*mu)*Identity(3)*theta

# Weak form a==l
a = inner(grad(delta_u), sigma(u, theta))*dx + kappa*dot(grad(delta_theta), grad(theta))*dx
l = dot(b, delta_u)*dx + dot(t_p, delta_u)*ds(top) + r*delta_theta*dx + q_p*delta_theta*ds(top)

################################
#### ASSEMBLE AND SOLVE ########
################################

y = Function(V)

print("Assembling... ")
A = assemble(a)
B = assemble(l)
for bc in bcs:
	bc.apply(A, B)
Y = y.vector()

print("Solving... ")
solver = LUSolver(A, "mumps")
solver.solve(Y, B)

(u, theta) = y.split()

################################
#### POST-PROCESSING ###########
################################

print("Postprocessing... ")

# Create displacement and temperature file
u.rename("u","displacement")
File("displacement.pvd", "compressed") << u

# Project stress field and create stress file
def dev(s):
	return s-tr(s)*Identity(d)/3.0
def von_mises(s):
	return sqrt(3.0/2.0*inner(dev(s), dev(s)))
S = FunctionSpace(mesh, "Lagrange", p)
stress = project(sigma(u, theta)[0,0]/1.e6, S)
stress.rename("stress", "xx")
#stress = project(von_mises(sigma(u))/1.e6, S)
#stress.rename("stress", "vonMises")
#T = TensorFunctionSpace(mesh, "Lagrange", p)
#stress = project(sigma(u, theta)/1.e6, T, solver_type="cg")
#stress.rename("sigma", "stress")

File("stress.pvd", "compressed") << stress

# Create temperature file
theta.rename("theta", "temperature")
File("temperature.pvd", "compressed") << theta

