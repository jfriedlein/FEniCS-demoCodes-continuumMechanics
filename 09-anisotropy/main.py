from fenics import *

################################
#### PROBLEM DEFINITION ########
################################

# Define geometry and mesh
d = 3
p0 = Point(0.0, 0.0, 0.0)
p1 = Point(1.0, 1.0, 1.0)
mesh = BoxMesh(p0, p1, 5, 5, 5)

# Define boundaries
boundaries = MeshFunction("size_t", mesh, d-1)
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
E_L = 44000*1.e6  #=E_11
E_T = 13000*1.e6  #=E_22=E_33
nu_LT = 0.25 #=nu=nu_12=nu_13
nu_TT = 0.3  #=nu_23=nu_32
mu_LT = 5600*1.e6 #=G_12=G_13
mu_TT = E_T/(2*(1+nu_TT)) #=G_23
nu_TL = E_T/E_L*nu_LT

lmbda = (nu_LT*nu_TL+nu_TT)/(1-nu_TT-2*nu_LT*nu_TL)/(1+nu_TT)*E_T
alpha = (nu_LT*(1+nu_TT-nu_TL)-nu_TT)/(1-nu_TT-2*nu_LT*nu_TL)/(1+nu_TT)*E_T
beta = E_L*(1-nu_TT*nu_TT)/(1-nu_TT-2*nu_LT*nu_TL)/(1+nu_TT)-lmbda-4*mu_LT

#a1 = Constant((1.0, 0.0, 0.0))
a1 = Constant((0.0, 1.0, 0.0))
#a1 = Constant((1.0/sqrt(2), 1.0/sqrt(2), 0.0))
A1 = outer(a1, a1)

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
	   DirichletBC(V.sub(1), Constant(0.0), AxisZ(), method="pointwise"),
	   DirichletBC(V.sub(2), Constant(0.0), Origin(), method="pointwise"),
	   DirichletBC(V.sub(0), Constant(0.02), boundaries, right)
       ]

# Stress tensor (linear isotropic elasticity)
def sigma(u):
    return lmbda*tr(sym(grad(u)))*Identity(d) + 2.0*mu_TT*sym(grad(u)) \
    	   + alpha*(dot(a1, dot(sym(grad(u)), a1))*Identity(d) + tr(sym(grad(u)))*A1) \
		   + 2.0*(mu_LT-mu_TT)*(dot(A1, sym(grad(u)))+dot(sym(grad(u)), A1)) \
		   + beta*dot(a1, dot(sym(grad(u)), a1))*A1
		   
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

#stress = project(sigma(u)[0,0]/1.e6, S)
#stress.rename("stress", "xx")
#stress = project(von_mises(sigma(u))/1.e6, S)
#stress.rename("stress", "vonMises")
stress = project(sigma(u)/1.e6, T, solver_type="cg", preconditioner_type="petsc_amg")
stress.rename("sigma", "stress")

File("stress.pvd", "compressed") << stress
