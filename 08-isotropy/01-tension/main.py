from dolfin import *

################################
#### PROBLEM DEFINITION ########
################################

# Define geometry and mesh
d = 3
p0 = Point(0.0, 0.0, 0.0)
p1 = Point(1.0, 1.0, 1.0)
mesh = BoxMesh(p0, p1, 3, 3, 3)

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
mu    = E/(2.0*(1.0 + nu))
lmbda = E*nu/((1.0 + nu)*(1.0 - 2.0*nu))
#mu    = E/(2.0*(1.0 + nu)) #PLANE STRESS
#lmbda = E*nu/((1.0 + nu)*(1.0 - nu)) #PLANE STRESS

#rho = 8.e3
#g = 9.81


# Volume force/ heat source and prescribed tractions/ prescribed heat fluxes
b = Constant((0.0, 0.0, 0.0))
t_p = Constant((0.0, 0.0, 0.0))

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
		return (abs(x[0]) < DOLFIN_EPS) and (abs(x[1]) < DOLFIN_EPS) and (abs(x[2]) < DOLFIN_EPS)

class AxisZ(SubDomain):
	def inside(self, x, on_boundary):
		return (abs(x[0]) < DOLFIN_EPS) and (abs(x[1]) < DOLFIN_EPS)
						       
bcs = [DirichletBC(V.sub(0), Constant(0.0), boundaries, left),
	   DirichletBC(V.sub(1), Constant(0.0), boundaries, bottom),
	   DirichletBC(V.sub(2), Constant(0.0), boundaries, back),
	   DirichletBC(V.sub(0), Constant(0.2), boundaries, right)
       ]

# Stress tensor (linear isotropic elasticity)
def sigma(u):
    return lmbda*tr(sym(grad(u)))*Identity(d) + 2.0*mu*sym(grad(u))
		   
# Material parameters
#E = 200.e9
#nu = 0.3
#mu    = E/(2.0*(1.0 + nu))
#lmbda = E*nu/((1.0 + nu)*(1.0 - 2.0*nu))

# Stress tensor (linear isotropic elasticity)
#def sigma(u):
#    return lmbda*tr(sym(grad(u)))*Identity(d) + 2.0*mu*sym(grad(u))
		   
# Weak form a==l
a = inner(grad(delta_u), sigma(u))*dx
l = dot(b, delta_u)*dx + dot(t_p, delta_u)*ds(right)

################################
#### ASSEMBLE AND SOLVE ########
################################

u = Function(V)
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

stress = project(sigma(u)/1.e6, T, solver_type="cg", preconditioner_type="petsc_amg")
stress.rename("sigma", "stress")

File("stress.pvd", "compressed") << stress

