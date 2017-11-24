from dolfin import *

################################
#### PROBLEM DEFINITION ########
################################

# Define geometry and mesh
d = 2
mesh = Mesh("lshape.xml.gz")
mesh = refine(mesh)

margin = 0.5
for i in range (0, 4):
	margin = margin*2.0/3.0
	markers = CellFunction("bool", mesh)
	markers.set_all(False)
	CompiledSubDomain("x[1]>(0.5-m) && x[1]<(0.5+m)", m=margin).mark(markers, True)
	mesh = refine(mesh, markers)
#for i in range (0, 2):
#	margin = 0.1
#	markers = CellFunction("bool", mesh)
#	markers.set_all(False)
#	CompiledSubDomain("x[0]>(0.5-m) && x[0]<(0.5+m) && x[1]>(0.5-m) && x[1]<(0.5+m)", m=margin).mark(markers, True)
#	mesh = refine(mesh, markers)


# Define boundaries
boundaries = FacetFunction("size_t", mesh)
boundaries.set_all(0)
left, right, bottom, top, = 1, 2, 3, 4
CompiledSubDomain("near(x[0], side) && on_boundary", side = 1.0).mark(boundaries, left)
CompiledSubDomain("near(x[0], side) && on_boundary", side = 1.0).mark(boundaries, right)
CompiledSubDomain("near(x[1], side) && on_boundary", side = 1.0).mark(boundaries, bottom)
CompiledSubDomain("near(x[1], side) && on_boundary", side = 1.0).mark(boundaries, top)

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
nu = 0.3
#mu    = E/(2.0*(1.0 + nu))
#lmbda = E*nu/((1.0 + nu)*(1.0 - 2.0*nu))
mu    = E/(2.0*(1.0 + nu)) #PLANE STRESS
lmbda = E*nu/((1.0 + nu)*(1.0 - nu)) #PLANE STRESS

#rho = 8.e3
#g = 9.81

# Volume force/ heat source and prescribed tractions/ prescribed heat fluxes
b = Constant((0.0, 0.0))
t_p = Constant((1.e6, -1.e6))

# Dirichlet boundary conditions
bcs = [DirichletBC(V, Constant((0.0, 0.0)), boundaries, top),
#       DirichletBC(V.sub(1), Constant((0.0)), origin_edge, method='pointwise'),
#	   DirichletBC(V.sub(2), Constant((0.0)), origin_point, method='pointwise')
       ]

# Stress tensor (linear isotropic elasticity)
def sigma(u):
    return lmbda*tr(sym(grad(u)))*Identity(d) + 2.0*mu*sym(grad(u))

# Weak form a==l
a = inner(grad(delta_u), sigma(u))*dx
l = dot(b, delta_u)*dx + dot(t_p, delta_u)*ds(right)

M = 0.5*inner(grad(u), sigma(u))*dx - dot(b, u)*dx - dot(t_p, u)*ds(right)

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
#solver = LUSolver(K, "mumps")
#solver.solve(U, F)

solve(a == l, u, bcs=bcs, 
	  solver_parameters={"linear_solver": "mumps"},
	  form_compiler_parameters={"optimize": True})
#	  tol=1.e-1, M=M)
	  

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

