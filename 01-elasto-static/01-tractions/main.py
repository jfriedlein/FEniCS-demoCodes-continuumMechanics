from fenics import *


#-------------------------------------------------
# Geometry and mesh generation
#-------------------------------------------------


d = 2 # spatial dimension
p0 = Point(0.0, 0.0) #	bottom-left corner
p1 = Point(3.0, 1.0) # top-right corner
mesh = RectangleMesh(p0, p1, 30, 10) # 30 and 10 are the number of elements in x- and y-directions

#-------------------------------------------------
# Boundary identification and measures
#-------------------------------------------------


boundaries = MeshFunction("size_t", mesh, d-1) # create mesh function for boundary domains
boundaries.set_all(0) # initialize all boundaries to 0
left, right, bottom, top = 1, 2, 3, 4 # define boundary IDs
CompiledSubDomain("near(x[0], side) && on_boundary", side = p0[0]).mark(boundaries, left)
CompiledSubDomain("near(x[0], side) && on_boundary", side = p1[0]).mark(boundaries, right)
CompiledSubDomain("near(x[1], side) && on_boundary", side = p0[1]).mark(boundaries, bottom)
CompiledSubDomain("near(x[1], side) && on_boundary", side = p1[1]).mark(boundaries, top)

ds = Measure('ds', domain=mesh, subdomain_data=boundaries) # define measure for boundary integration


#-------------------------------------------------
# Function spaces and variational functions
#-------------------------------------------------


p = 2 # polynomial degree
V = VectorFunctionSpace(mesh, "Lagrange", p) # define function space for displacement field

u = TrialFunction(V) # define trial and test functions
delta_u = TestFunction(V)


#-------------------------------------------------
# Material properties
#-------------------------------------------------


E = 200.e9 # Young's modulus in Pa
rho = 8.e3 # density in kg/m^3
g = 9.81 # acceleration due to gravity in m/s^2
nu = 0.3 # Poisson's ratio
mu    = E/(2.0*(1.0 + nu)) # Lame's first parameter
lmbda = E*nu/((1.0 + nu)*(1.0 - nu)) # Lame's second parameter


#-------------------------------------------------
# Loads and boundary conditions
#-------------------------------------------------


b = Constant((0.0, 0.0)) # body force
t_p = Constant((0.0, -1.e6)) # traction force
u_zero = Constant((0.0, 0.0)) # zero displacement

bcs = [DirichletBC(V, u_zero, boundaries, left), # fixed left edge
	   DirichletBC(V,u_zero, boundaries, right) # fixed right edge
       ]


#-----------------------------------------------
#  Variational formulation of linear elasticity
#-----------------------------------------------


def sigma(u):
    return lmbda*tr(sym(grad(u)))*Identity(d) + 2.0*mu*sym(grad(u)) # stress tensor

a = inner(sym(grad(delta_u)),sigma(u))*dx # bilinear form
l = dot(b, delta_u)*dx + dot(t_p, delta_u)*ds(top) # linear form


#---------------------------------------------
# Solve the linear system
#---------------------------------------------


u = Function(V) # solution function
solve(a == l, u, bcs=bcs, 
	      solver_parameters={"linear_solver": "mumps"},
		  form_compiler_parameters={"optimize": True}) # solve the variational problem


#-------------------------------------------------
# Save displacement
#-------------------------------------------------


u.rename("u", "displacement") # rename solution for output
File("displacement.pvd", "compressed") << u # save displacement to file


#-------------------------------------------------
# Compute and save stress
#-------------------------------------------------


T = TensorFunctionSpace(mesh, "Lagrange", p) # function space for stress output
stress = project(sigma(u)/1.e6, T, solver_type="mumps")

stress.rename("sigma", "stress") # rename stress for output
File("stress.pvd", "compressed") << stress # save stress to file