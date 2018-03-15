#! /usr/bin/env python

from __future__ import print_function
from fenics import *

parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["optimize"] = True

# Time definitions
t = 0.0
T = 1.0
num_steps = 200
dt = T / num_steps

# Parameters for Newmark (Average acceleration)/Crank-Nicolson
cn_alpha = 0.5
nm_beta = 0.25
nm_gamma = 0.5

################################
#### PROBLEM DEFINITION ########
################################

# Define geometry and mesh
p0 = Point(0.0, 0.0, 0.0)
p1 = Point(0.4, 0.1, 0.1)
mesh = RectangleMesh(p0, p1, 100, 25)

# Define boundaries
boundaries = MeshFunction("size_t", mesh, d-1)
boundaries.set_all(0)
left, right, bottom, top = 1, 2, 3, 4
CompiledSubDomain("near(x[0], side) && on_boundary", side = p0[0]).mark(boundaries, left)
CompiledSubDomain("near(x[0], side) && on_boundary", side = p1[0]).mark(boundaries, right)
CompiledSubDomain("near(x[1], side) && on_boundary", side = p0[1]).mark(boundaries, bottom)
CompiledSubDomain("near(x[1], side) && on_boundary", side = p1[1]).mark(boundaries, top)

# Coordinates and surface integral element
x = SpatialCoordinate(mesh)
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

# Define function space
d = 2
p = 2
P = VectorElement("Lagrange", mesh.ufl_cell(), p, d)
Q = FiniteElement("Lagrange", mesh.ufl_cell(), p)
V = FunctionSpace(mesh, P*Q)

# Define trial and test functions
(u, theta) = TrialFunctions(V)
(delta_u, delta_theta) = TestFunctions(V)

# Collapse mixed function space
V_u = V.sub(0).collapse()
V_theta = V.sub(1).collapse()

# Interpolate initial displacement, velocity and temperature
u_n = Function(V_u)
u_n.interpolate(Constant((0.0, 0.0)))
v_n = Function(V_u)
v_n.interpolate(Constant((0.0, 0.0)))
#v_n.interpolate(Expression("x[0]*x[0]/l/l*0.1"), l=p1[0], degree=1))

theta0 = 25.0
theta_n = Function(V_theta)
theta_n.interpolate(Constant(0.0))

# Define additional rates and predictor functions for Newmark
a_n = Function(V_u)
dtheta_n = Function(V_theta)
u_pred = Function(V_u)
v_pred = Function(V_u)
theta_pred = Function(V_theta)

# Material parameters
E = 200.e9
nu = 0.3
mu    = E/(2.0*(1.0 + nu))
lmbda = E*nu/((1.0 + nu)*(1.0 - 2.0*nu))
#mu    = E/(2.0*(1.0 + nu)) #PLANE STRESS
#lmbda = E*nu/((1.0 + nu)*(1.0 - nu)) #PLANE STRESS
rho = 8.e3
g = 9.81

kappa = 80.0
alpha = 12.e-6
#alpha = 0.0
beta = alpha*E/(1-2.0*nu)
cv = rho*400

# Volume force/ heat source and prescribed tractions/ prescribed heat fluxes
#b = as_vector((0.0, -rho*g, 0.0))
b = Constant((0.0, 0.0))
#t_p = Constant((0.0, 0.0))
t_p = Expression(("m*sin(2*3.141*t/t0)", "0.0"), degree=4, t0=T, m=1.e8, t=0)
#t_p = Expression(("(t<3.0*t1)?1.0*(-m/t1*fabs(t-t1)+m):0.0", "0.0"), degree=2, t1=0.25*T, m=1.e9, t=0)
#t_p = Expression(("(t<t1)?m:0.0", "0.0"), degree=2, t1=0.33*T, m=1.e9, t=0)
r = Constant(0.0)
q_p = Constant(0.0)

# Prescribed Dirichlet boundary data
#u_p = Expression(("0.03*sin(2*3.141*t/t0)", "0.0"), t=t, t0=T, degree=3)
#theta_p = Expression(('m*t'), degree=1, m=100.0/T, t=0)

# Dirichlet boundary conditions
bcs = [DirichletBC(V.sub(0), Constant((0.0, 0.0)), boundaries, left),
       #DirichletBC(V.sub(0), u_p, boundaries, right),
       #DirichletBC(V.sub(0).sub(0), u0, boundaries, right),
       #DirichletBC(V.sub(1), theta_p, boundaries, left)]#,
       #DirichletBC(V.sub(1), Constant(theta0), boundaries, right)
       ]
    
# Stress tensor (linear isotropic elasticity)
def sigma(u, theta):
    return lmbda*tr(sym(grad(u)))*Identity(d) + 2.0*mu*sym(grad(u)) - alpha*(3.0*lmbda + 2.0*mu)*Identity(d)*theta


#dirname = "ev"+str(j)+"/"
#if not path.exists(dirname):
#	mkdir(dirname)


# Newmark approximations for rates
def a(u, u_pred):
	return (u-u_pred)/(dt*dt*nm_beta)
def v(u, u_pred, v_pred):
	return v_pred+dt*nm_gamma*a(u, u_pred)
def dtheta(th, th_pred):
	return (th-th_pred)/(dt*cn_alpha)

# Form definition (F(u, delta_u) = a(u, delta_u) - l(delta_u))
F = rho*dot(delta_u, a(u, u_pred))*dx \
	+ cv*dtheta(theta, theta_pred)*delta_theta*dx \
	+ inner(grad(delta_u), sigma(u, theta))*dx \
	+ kappa*dot(grad(delta_theta), grad(theta))*dx \
    + theta0*beta*tr(sym(grad(v(u, u_pred, v_pred))))*delta_theta*dx \
	- dot(b, delta_u)*dx - dot(t_p, delta_u)*ds(right) \
	- r*delta_theta*dx - q_p*delta_theta*ds(right)

# Project initial stress field
def dev(s):
	return s-tr(s)*Identity(d)/3.0
def von_mises(s):
	return sqrt(3.0/2.0*inner(dev(s), dev(s)))
#V_stress = TensorFunctionSpace(mesh, "Lagrange", p)
#stress_n = Function(V_stress)
#stress_n = project(sigma(u_n, theta_n)/1.e6, V_stress, solver_type="mumps")

# Create output files
u_n.rename("u", "displacement")
file_u = File("displacement.pvd", "compressed")
file_u << (u_n, t)
theta_n.rename("theta","temperature")
file_theta = File("temperature.pvd", "compressed")
file_theta << (theta_n, t)
#stress_n.rename("stress","vonMises")
#file_stress = File("stress.pvd", "compressed")
#file_stress << (stress_n, t)


################################
#### ASSEMBLE AND SOLVE ########
################################

A = assemble(lhs(F))
[bc.apply(A) for bc in bcs]

print("Setup solver...")
solver = LUSolver(A, "mumps")
#solver.parameters["symmetric"] = True
solver.parameters["reuse_factorization"] = True
y = Function(V)

# Time integration and function solving
for n in range(num_steps):	
    t += dt
    print("Step ", n, " (t=", t, ")", sep="")
    
    # Update loads, boundary data, ...
    #u_p.t = t
    #theta_p.t = t
    t_p.t = t

    # Define predictors
    u_pred.vector()[:] = u_n.vector()+dt*v_n.vector()\
    					 +0.5*dt*dt*(1-2*nm_beta)*a_n.vector()
    v_pred.vector()[:] = v_n.vector()+dt*(1-nm_gamma)*a_n.vector()
    theta_pred.vector()[:] = theta_n.vector()+dt*(1-cn_alpha)*dtheta_n.vector()
    
    # Assemble and solve
    b = assemble(rhs(F))
    [bc.apply(b) for bc in bcs]
    solver.solve(y.vector(), b)
        
    # Update and correct functions
    (u, theta) = y.split(deepcopy=True)
    a_n.vector()[:] = (u.vector()-u_pred.vector())/(nm_beta*dt*dt)
    v_n.vector()[:] = v_pred.vector()+nm_gamma*dt*a_n.vector()
    u_n.vector()[:] = u.vector()
    dtheta_n.vector()[:] = (theta.vector()-theta_pred.vector())/(cn_alpha*dt)    
    theta_n.vector()[:] = theta.vector()
    
    # Project stresses
    #stress_n.assign(project(sigma(u_n, theta_n)/1.e6, V_stress, solver_type="mumps"))
    
    # Write step to files
    file_u << (u_n, t)
    file_theta << (theta_n, t)
    #file_stress << (stress_n, t)

