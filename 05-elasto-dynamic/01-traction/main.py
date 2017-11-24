#! /usr/bin/env python

from __future__ import print_function
from dolfin import *

parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["optimize"] = True

# Time definitions
t = 0.0
T = 0.001
num_steps = 200
dt = T / num_steps

#Parameters for Newmark (Average acceleration)
nm_beta = 0.25
nm_gamma = 0.5


################################
#### PROBLEM DEFINITION ########
################################

# Define geometry and mesh
p_llb = Point(0.0, 0.0, 0.0)
p_rtf = Point(0.04, 0.01, 0.01)
mesh = RectangleMesh(p_llb, p_rtf, 20, 5)
#mesh = BoxMesh(p_llb, p_rtf, 4, 1, 1)

# Define boundaries
boundaries = FacetFunction("size_t", mesh)
boundaries.set_all(0)
left, right, bottom, top = 1, 2, 3, 4
CompiledSubDomain("near(x[0], side) && on_boundary", side = p_llb[0]).mark(boundaries, left)
CompiledSubDomain("near(x[0], side) && on_boundary", side = p_rtf[0]).mark(boundaries, right)
CompiledSubDomain("near(x[1], side) && on_boundary", side = p_llb[1]).mark(boundaries, bottom)
CompiledSubDomain("near(x[1], side) && on_boundary", side = p_rtf[1]).mark(boundaries, top)

# Coordinates and surface integral element
x = SpatialCoordinate(mesh)
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

# Define function space
d = 2
p = 2
V = VectorFunctionSpace(mesh, "Lagrange", p)

# Define trial and test functions
u = TrialFunction(V)
delta_u = TestFunction(V)

# Interpolate initial displacement, velocity
u_n = Function(V)
u_n.interpolate(Constant((0.0, 0.0)))
v_n = Function(V)
v_n.interpolate(Constant((0.0, 0.0)))

# Define additional rates and predictor functions for Newmark
a_n = Function(V)
u_pred = Function(V)
v_pred = Function(V)

# Material parameters
E = 200.e9
nu = 0.45
mu    = E/(2.0*(1.0 + nu))
lmbda = E*nu/((1.0 + nu)*(1.0 - 2.0*nu))
#mu    = E/(2.0*(1.0 + nu)) #PLANE STRESS
#lmbda = E*nu/((1.0 + nu)*(1.0 - nu)) #PLANE STRESS

rho = 8.e3
g = 9.81

# Volume force/ heat source and prescribed tractions
#b = as_vector((0.0, -rho*g))
b = Constant((0.0, 0.0))
t_p = Expression(('(t<1.0*t1)?10.0*(-m/t1*fabs(t-t1)+m):0.0', '(t<1.0*t1)?-(-m/t1*fabs(t-t1)+m):0.0'), degree=2, t1=0.1*T, m=1.e5, t=0)

# Prescribed Dirichlet boundary data
#u_p = Expression(('m*t', '0.0', '0.0'), degree=1, m=-0.2/T, t=0)
#theta_p = Expression(('m*t'), degree=1, m=100.0/T, t=0)

# Dirichlet boundary conditions
bcs = [DirichletBC(V, Constant((0.0, 0.0)), boundaries, left)]
    
# Stress tensor (linear isotropic elasticity)
def sigma(u):
    return lmbda*tr(sym(grad(u)))*Identity(d) + 2.0*mu*sym(grad(u))

################################
#### WEAK FORM (+NEWMARK) ######
################################

# Newmark approximations for rates
def a(u, u_pred):
	return (u-u_pred)/(dt*dt*nm_beta)
def v(u, u_pred, v_pred):
	return v_pred+dt*nm_gamma*a(u, u_pred)

# Form definition (F(u, delta_u) = a(u, delta_u) - l(delta_u))
F = rho*dot(delta_u, a(u, u_pred))*dx \
	+ inner(grad(delta_u), sigma(u))*dx \
	- dot(b, delta_u)*dx - dot(t_p, delta_u)*ds(right)

# Project initial stress field
def dev(s):
	return s-tr(s)*Identity(3)/3.0
def von_mises(s):
	return sqrt(3.0/2.0*inner(dev(s), dev(s)))
Z = TensorFunctionSpace(mesh, "Lagrange", p)
stress_n = Function(Z)
stress_n = project(sigma(u_n)/1.e6, Z, solver_type="mumps")
	
# Create output files
u_n.rename("u","displacement")
file_u = File("displacement.pvd", "compressed")
file_u << (u_n, t)
stress_n.rename("sigma","stress")
file_stress = File("stress.pvd", "compressed")
file_stress << (stress_n, t)

################################
#### ASSEMBLE AND SOLVE ########
################################

A = assemble(lhs(F))
[bc.apply(A) for bc in bcs]

print("Setup solver...")
solver = LUSolver(A, "mumps")
solver.parameters["symmetric"] = True
solver.parameters["reuse_factorization"] = True
u = Function(V)

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
    
    # Assemble and solve
    b = assemble(rhs(F))
    [bc.apply(b) for bc in bcs]
    solver.solve(u.vector(), b)
        
    # Update and correct functions
    #(u, theta) = y.split(deepcopy=True)
    a_n.vector()[:] = (u.vector()-u_pred.vector())/(nm_beta*dt*dt)
    v_n.vector()[:] = v_pred.vector()+nm_gamma*dt*a_n.vector()
    u_n.vector()[:] = u.vector()
    
    # Project stresses
    stress_n.assign(project(sigma(u_n)/1.e6, Z, solver_type="mumps"))
    
    # Write step to files
    file_u << (u_n, t)
    file_stress << (stress_n, t)
