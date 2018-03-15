from fenics import *

################################
#### PROBLEM DEFINITION ########
################################

# Define geometry and mesh
d = 3
p0 = Point(0.0, 0.0, 0.0)
p1 = Point(1.5, 1.0, 1.0)
mesh = BoxMesh(p0, p1, 15, 10, 10)

# Define boundaries
boundaries = MeshFunction("size_t", mesh, d-1)
boundaries.set_all(0)
left, right, bottom, top = 1, 2, 3, 4
CompiledSubDomain("near(x[0], side) && on_boundary", side = p0[0]).mark(boundaries, left)
CompiledSubDomain("near(x[0], side) && on_boundary", side = p1[0]).mark(boundaries, right)
CompiledSubDomain("near(x[1], side) && on_boundary", side = p0[1]).mark(boundaries, bottom)
CompiledSubDomain("near(x[1], side) && on_boundary", side = p1[1]).mark(boundaries, top)

# Surface integral element
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

# Define function space
p = 2
V = VectorFunctionSpace(mesh, "Lagrange", p)

# Define trial and test functions
u = TrialFunction(V)
delta_u = TestFunction(V)

# Material parameters
E = 200.e9
rho = 8.e3
g = 9.81
nu = 0.3
mu    = E/(2.0*(1.0 + nu))
lmbda = E*nu/((1.0 + nu)*(1.0 - 2.0*nu))

# Volume force/ heat source and prescribed tractions/ prescribed heat fluxes
b = Constant((0.0, 0.0, 0.0))
t_p = Constant((0.0, -1.e6, 0.0))

# Dirichlet boundary conditions
bcs = [DirichletBC(V, Constant((0.0, 0.0, 0.0)), boundaries, left),
	   DirichletBC(V.sub(0), Constant((0.0)), boundaries, right)
       ]

# Stress tensor (linear isotropic elasticity)
def sigma(u):
    return lmbda*tr(sym(grad(u)))*Identity(d) + 2.0*mu*sym(grad(u))

# Weak form a==l
a = inner(grad(delta_u), sigma(u))*dx
l = dot(b, delta_u)*dx + dot(t_p, delta_u)*ds(top)

################################
#### ASSEMBLE AND SOLVE ########
################################


u = Function(V)

#K = assemble(a)
#F = assemble(l)
#for bc in bcs:
#	bc.apply(K, F)
#U = u.vector()
#solve(K, U, F)

solve(a == l, u, bcs=bcs, 
	      solver_parameters={"linear_solver": "mumps"},
	      form_compiler_parameters={"optimize": True})

################################
#### POST-PROCESSING ###########
################################

# Create displacement and temperature file
u.rename("u", "displacement")
File("displacement.pvd", "compressed") << u

# Project stress field and create stress file
def dev(s):
	return s-tr(s)*Identity(d)/3.0
def von_mises(s):
	return sqrt(3.0/2.0*inner(dev(s), dev(s)))
S = FunctionSpace(mesh, "Lagrange", p)
T = TensorFunctionSpace(mesh, "Lagrange", p)

#stress = project(sigma(u)[0,0]/1.e6, S)
#stress.rename("stress", "xx")
#stress = project(von_mises(sigma(u))/1.e6, S)
#stress.rename("stress", "vonMises")
stress = project(sigma(u)/1.e6, T, solver_type="mumps")
stress.rename("sigma", "stress")

File("stress.pvd", "compressed") << stress
