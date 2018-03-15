#! /usr/bin/env python

from __future__ import print_function
from fenics import *
import numpy as np

# Time definitions
t = 0.0
T = 3600
num_steps = 200
dt = T / num_steps

# Parameters Crank-Nicolson
cn_alpha = 0.5

################################
#### PROBLEM DEFINITION ########
################################

# Define geometry and mesh
d = 2
p0 = Point(0.0, 0.0, 0.0)
p1 = Point(0.1, 0.1, 0.1)
#mesh = BoxMesh(p0, p1, 20, 10, 10)
mesh = RectangleMesh(p0, p1, 50, 50)

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
p = 2
V = FunctionSpace(mesh, "Lagrange", p)

# Define trial and test functions
theta = TrialFunction(V)
delta_theta = TestFunction(V)

# Interpolate initial temperature
theta_n = Function(V)
theta_n.interpolate(Constant(25.0))
#theta_n.interpolate(Expression("1.0e2*exp(-((x[0]-rx)*(x[0]-rx)+(x[1]-ry)*(x[1]-ry))/0.001)", rx=0.0, ry=0.5, degree=4))

# Define additional rates and predictor functions for Newmark
dtheta_n = Function(V)
theta_pred = Function(V)

# Material parameters
#rho = 1.0*conditional(x[1] < 0.5*p1[1], 7800.0, 7000.0)
#kappa = 1.0*conditional(x[1] < 0.5*p1[1], 81.0, 121.0)
#cv = 1.0*rho*conditional(x[1] < 0.5*p1[1], 452.0, 387.0)
rho = 8000.0
kappa = 20.0
cv = rho*500.0

# Heat source and prescribed heat fluxes
#r = Expression("1.0e6*sin(m*t)", m=2*pi/T, t=0, degree=1)
#r = Expression("1.0e7*exp(-((x[0]-x0)*(x[0]-x0)+(x[1]-y0)*(x[1]-y0))/0.001)*m*t", x0=0.25, y0=0.25, m=1/T, t=0.0, degree=4)
#q_p = Expression("1.0e3*m*t", m=1/T, t=0, degree=1)
r = Constant(0.0)
q_p = Expression("x[1]*1.e5", degree=4)
#q_p_bottom = Expression("x[0]<0.5 ? -1.e3*m*t : 0.0", m=1/T, t=0.0, degree=4)

# Dirichlet boundary conditions
bcs = [DirichletBC(V, Constant(25.0), boundaries, left)]

# Crank-N. approximations for rates
def dtheta(th, th_pred):
	return (th-th_pred)/(dt*cn_alpha)

# Form definition (F(u, delta_u) = a(u, delta_u) - l(delta_u))
F = cv*dtheta(theta, theta_pred)*delta_theta*dx \
	+ kappa*dot(grad(delta_theta), grad(theta))*dx \
	- r*delta_theta*dx \
	- q_p*delta_theta*ds(right) # - q_p_bottom*delta_theta*ds(bottom)

# Create output files
theta_n.rename("theta","temperature")
file_theta = File("temperature.pvd", "compressed")
file_theta << (theta_n, t)


################################
#### ASSEMBLE AND SOLVE ########
################################

A = assemble(lhs(F))
[bc.apply(A) for bc in bcs]

solver = LUSolver(A, 'mumps')
solver.parameters["symmetric"] = True
solver.parameters['reuse_factorization'] = True
theta = Function(V)

# Time integration and function solving
for n in range(num_steps):	
    t += dt
    print("Step ", n, " (t=", t, ")", sep="")

    # Update loads, boundary data, ...
    #u_p.t = t
    #theta_p.t = t
    #t_p.t = t
#    r.rx = (p1[0]-p0[0])*t/(0.01*T)
#    r.ry = (p1[1]-p0[1])*t/(0.01*T)

    #r.rx = 0.3*np.cos(2*pi*t/(0.3*T))+0.5
    #r.ry = 0.3*np.sin(2*pi*t/(0.4*T))+0.5
#    if t<0.5*T:
#        q_p.s = t*2/T
#    else:
#        q_p.m = 1.0
    #q_p.t = t
    #q_p_bottom.t = t
    
    # Define predictors
    theta_pred.vector()[:] = theta_n.vector()+dt*(1-cn_alpha)*dtheta_n.vector()
    
    # Assemble and solve
    b = assemble(rhs(F))
    [bc.apply(b) for bc in bcs]
    solver.solve(theta.vector(), b)
        
    # Update and correct functions
    dtheta_n.vector()[:] = (theta.vector()-theta_pred.vector())/(cn_alpha*dt)    
    theta_n.vector()[:] = theta.vector()    
    
    # Write step to files
    file_theta << (theta_n, t)
