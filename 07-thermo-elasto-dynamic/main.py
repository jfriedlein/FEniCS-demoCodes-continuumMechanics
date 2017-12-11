#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from dolfin import *

# Time definitions
t = 0.0
T = 10.0
num_steps = 100
dt = T / num_steps

# Parameters for Newmark (Average acceleration)/Crank-Nicolson
cn_alpha = 0.5
nm_beta = 0.25
nm_gamma = 0.5



################################
#### PROBLEM DEFINITION ########
################################

# Define geometry and mesh
p = Point(0.0, 0.0, 0.0)
q = Point(0.2, 0.05, 0.1)
mesh = BoxMesh(p, q, 20, 10, 10)

# Define boundaries
boundaries = FacetFunction("size_t", mesh)
boundaries.set_all(0)
left, right, bottom, top = 1, 2, 3, 4
CompiledSubDomain("near(x[0], side) && on_boundary", side = p[0]).mark(boundaries, left)
CompiledSubDomain("near(x[0], side) && on_boundary", side = q[0]).mark(boundaries, right)
CompiledSubDomain("near(x[1], side) && on_boundary", side = p[1]).mark(boundaries, bottom)
CompiledSubDomain("near(x[1], side) && on_boundary", side = q[1]).mark(boundaries, top)

# Coordinates and surface integral element
x = SpatialCoordinate(mesh)
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

# Define function space
P = VectorElement("Lagrange", mesh.ufl_cell(), 2, 3)
Q = FiniteElement("Lagrange", mesh.ufl_cell(), 2)
V = FunctionSpace(mesh, P*Q)

# Define trial and test functions
(u, theta) = TrialFunctions(V)
(delta_u, delta_theta) = TestFunctions(V)

# Collapse mixed function space
V_u = V.sub(0).collapse()
V_theta = V.sub(1).collapse()

# Interpolate initial displacement, velocity and temperature
u_n = Function(V_u)
u_n.interpolate(Constant((0.0, 0.0, 0.0)))
v_n = Function(V_u)
v_n.interpolate(Constant((0.0, 0.0, 0.0)))
theta_n = Function(V_theta)
theta_n.interpolate(Constant(0.0))

# Define additional rates and predictor functions for Newmark
a_n = Function(V_u)
dtheta_n = Function(V_theta)
u_pred = Function(V_u)
v_pred = Function(V_u)
theta_pred = Function(V_theta)

# Material parameters
rho = 1.0*conditional(x[1] < 0.5*q[1], 7874.0, 7140.0)
g = 9.81
E  = 1.0*conditional(x[1] < 0.5*q[1], 120.0e9, 92.0e9)
nu = 0.3
mu    = E/(2.0*(1.0 + nu))
lmbda = E*nu/((1.0 + nu)*(1.0 - 2.0*nu))
alpha = 1.0*conditional(x[1] < 0.5*q[1], 11.8e-6, 30.2e-6)
beta = alpha*E/(1-2.0*nu)
theta0 = 273.0
kappa = 1.0*conditional(x[1] < 0.5*q[1], 80.0, 112.0)
#cv = 461.0*rho#*conditional(x[1] < 0.5*q[1], 10.0, 1.0)
cv = 1.0*rho*conditional(x[1] < 0.5*q[1], 449.0, 388.0)

# Volume force/ heat source and prescribed tractions/ prescribed heat fluxes
#b = as_vector((0.0, -rho*g, 0.0))
b = Constant((0.0, 0.0, 0.0))
t_p = Constant((0.0, 0.0, 0.0))
r = Constant(1.0e6)
q_p = Constant(0.0)

# Prescribed Dirichlet boundary data
#u_p = Expression(('m*t', '0.0', '0.0'), degree=1, m=-0.2/T, t=0)
#theta_p = Expression(('m*t'), degree=1, m=100.0/T, t=0)

# Dirichlet boundary conditions
bcs = [DirichletBC(V.sub(0), Constant((0.0, 0.0, 0.0)), boundaries, left),
       #DirichletBC(V.sub(0), u_p, boundaries, right),
       #DirichletBC(V.sub(1), theta_p, boundaries, left)]#,
       DirichletBC(V.sub(1), Constant(0.0), boundaries, left)]
    
# Stress tensor (linear isotropic elasticity)
def sigma(u, theta):
    return lmbda*tr(sym(grad(u)))*Identity(3) + 2.0*mu*sym(grad(u)) - alpha*(3.0*lmbda + 2.0*mu)*Identity(3)*theta

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
	- dot(b, delta_u)*dx - dot(t_p, delta_u)*ds(top) \
	- r*delta_theta*dx - q_p*delta_theta*ds(top)

# Project initial stress field
def dev(s):
	return s-tr(s)*Identity(3)/3.0
def von_mises(s):
	return sqrt(3.0/2.0*inner(dev(s), dev(s)))
V_stress = FunctionSpace(mesh, "Lagrange", 1)
stress_n = Function(V_stress)
stress_n.assign(project(von_mises(sigma(u_n, theta_n)), V_stress))

# Create output files
u_n.rename("u","displacement")
file_u = File("thermoelasticitydynamic_newmark_displacement.pvd", "compressed")
file_u << (u_n, t)
theta_n.rename("theta","temperature")
file_theta = File("thermoelasticitydynamic_newmark_temperature.pvd", "compressed")
file_theta << (theta_n, t)
stress_n.rename("stress","vonMises")
file_stress = File("thermoelasticitydynamic_newmark_stress.pvd", "compressed")
file_stress << (stress_n, t)



################################
#### ASSEMBLE AND SOLVE ########
################################

A = assemble(lhs(F))
solver = LUSolver(A, 'mumps')
solver.parameters['reuse_factorization'] = True
y = Function(V)

# Time integration and function solving
for n in range(num_steps):	
    t += dt
    print("Step ", n, " (t=", t, ")", sep="")
    
    # Update loads, boundary data, ...
    #u_p.t = t
    #theta_p.t = t

    # Define predictors
    u_pred.vector()[:] = u_n.vector()+dt*v_n.vector()\
    					 +0.5*dt*dt*(1-2*nm_beta)*a_n.vector()
    v_pred.vector()[:] = v_n.vector()+dt*(1-nm_gamma)*a_n.vector()
    theta_pred.vector()[:] = theta_n.vector()+dt*(1-cn_alpha)*dtheta_n.vector()
    
    # Assemble and solve
    b = assemble(rhs(F))
    [bc.apply(A, b) for bc in bcs]
    solver.solve(y.vector(), b)
        
    # Update and correct functions
    (u, theta) = y.split(deepcopy=True)
    a_n.vector()[:] = (u.vector()-u_pred.vector())/(nm_beta*dt*dt)
    v_n.vector()[:] = v_pred.vector()+nm_gamma*dt*a_n.vector()
    u_n.vector()[:] = u.vector()
    dtheta_n.vector()[:] = (theta.vector()-theta_pred.vector())/(cn_alpha*dt)    
    theta_n.vector()[:] = theta.vector()
    
    # Project stresses
    stress_n.assign(project(von_mises(sigma(u_n, theta_n)), V_stress))
    
    # Write step to files
    file_u << (u_n, t)
    file_theta << (theta_n, t)
    file_stress << (stress_n, t)
