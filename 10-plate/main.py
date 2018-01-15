from fenics import *
from mshr import *
import numpy as np


class Stresses(Expression):
    def __init__(self, **kwargs):
        #Expression.__init__(self)
        self.a = 1
        self.b = 1
        self.m = (self.a-self.b)/(self.a+self.b)
        self.beta = 0
        self.p = 1
        
    def value_shape(self):
        return (3,)
        
    def eval(self, value, x):
        ure = (x[0]*x[0]-x[1]*x[1])-(self.a*self.a-self.b*self.b)
        urho = np.sqrt(ure*ure+4.0*x[0]*x[0]*x[1]*x[1])
        zetare = (x[0]+np.sign(x[0])*np.sqrt(0.5*np.maximum(urho+ure, 0.0)))/(self.a+self.b)
        zetaim = (x[1]+np.sign(x[1])*np.sqrt(0.5*np.maximum(urho-ure, 0.0)))/(self.a+self.b)
        rho = np.sqrt(zetare*zetare+zetaim*zetaim)
        theta = np.arctan2(zetaim, zetare)
        
        rho2 = rho*rho
        rho4 = rho2*rho2
        rho6 = rho4*rho2
        
        S1 = (rho4-2.0*rho2*np.cos(2.0*theta-2.0*self.beta)+2.0*self.m*np.cos(2.0*self.beta)-self.m*self.m)/(rho4-2.0*self.m*rho2*np.cos(2.0*theta)+self.m*self.m)
	    
        S2zr = (rho6*np.cos(-6.0*theta+2.0*self.beta) - 2.0*self.m*rho2*(np.cos(-4.0*theta-2.0*self.beta) - self.m*np.cos(-4.0*theta)) - 3.0*self.m*rho4*np.cos(-4.0*theta+2.0*self.beta) - (self.m*self.m-2.0*np.cos(2.0*self.beta)*self.m+1.0)*rho4*np.cos(4.0*theta)	- 2.0*rho4*(np.cos(-2.0*theta-2.0*self.beta) - self.m*np.cos(-2.0*theta)) + 3.0*rho2*np.cos(2.0*theta+2.0*self.beta) - self.m*(self.m*self.m-2.0*np.cos(2.0*self.beta)*self.m+1.0)*rho2*np.cos(2.0*theta) - self.m*np.cos(2.0*self.beta))
        S2zi = (rho6*np.sin(-6.0*theta+2.0*self.beta) - 2.0*self.m*rho2*(np.sin(-4.0*theta-2.0*self.beta) - self.m*np.sin(-4.0*theta)) - 3.0*self.m*rho4*np.sin(-4.0*theta+2.0*self.beta) - (self.m*self.m-2.0*np.cos(2.0*self.beta)*self.m+1.0)*rho4*np.sin(-4.0*theta) - 2.0*rho4*(np.sin(-2.0*theta-2.0*self.beta) - self.m*np.sin(-2.0*theta)) + 3.0*rho2*np.sin(-2.0*theta-2.0*self.beta) - self.m*(self.m*self.m-2.0*np.cos(2.0*self.beta)*self.m+1.0)*rho2*np.sin(-2.0*theta) - self.m*np.sin(-2.0*self.beta))
        
        S2nr = rho6*np.cos(6.0*theta)-3.0*self.m*rho4*np.cos(4.0*theta)+3.0*self.m*self.m*rho2*np.cos(2.0*theta)-self.m*self.m*self.m
        S2ni = -rho6*np.sin(6.0*theta)+3.0*self.m*rho4*np.sin(4.0*theta)-3.0*self.m*self.m*rho2*np.sin(2.0*theta)
        
        ReS2 = (S2zr*S2nr+S2zi*S2ni)/(S2nr*S2nr+S2ni*S2ni)
        ImS2 = (S2zi*S2nr-S2zr*S2ni)/(S2nr*S2nr+S2ni*S2ni)
        
        sigmaxx = self.p*0.5*(S1+ReS2)
        sigmayy = self.p*0.5*(S1-ReS2)
        sigmaxy = self.p*0.5*(ImS2)
        
        value[0] = np.cos(self.beta)*np.cos(self.beta)*sigmaxx+np.sin(self.beta)*np.sin(self.beta)*sigmayy+2.0*np.cos(self.beta)*np.sin(self.beta)*sigmaxy
        value[1] = np.sin(self.beta)*np.sin(self.beta)*sigmaxx+np.cos(self.beta)*np.cos(self.beta)*sigmayy-2.0*np.cos(self.beta)*np.sin(self.beta)*sigmaxy
        value[2] = -np.cos(self.beta)*np.sin(self.beta)*sigmaxx+np.sin(self.beta)*np.cos(self.beta)*sigmayy+(np.cos(self.beta)*np.cos(self.beta)-np.sin(self.beta)*np.sin(self.beta))*sigmaxy
        
    


################################
#### PROBLEM DEFINITION ########
################################

# Define geometry and mesh
d = 2

# Parameters
R = 1
L = 3

# Create geometry
plate = Rectangle(Point(-L, 0), Point(0, L))
hole = Circle(Point(0, 0), R)
geometry = plate - hole

# Create mesh
mesh = generate_mesh(geometry, 16)
#mesh = refine(mesh)

#margin = 0.5
#for i in range (0, 4):
#	margin = margin*2.0/3.0
#	markers = CellFunction("bool", mesh)
#	markers.set_all(False)
#	CompiledSubDomain("x[1]>(0.5-self.m) && x[1]<(0.5+self.m)", self.m=margin).mark(markers, True)
#	mesh = refine(mesh, markers)

##for i in range (0, 2):
##	margin = 0.1
##	markers = CellFunction("bool", mesh)
##	markers.set_all(False)
##	CompiledSubDomain("x[0]>(0.5-self.m) && x[0]<(0.5+self.m) && x[1]>(0.5-self.m) && x[1]<(0.5+self.m)", self.m=margin).mark(markers, True)
##	mesh = refine(mesh, markers)


# Define boundaries
boundaries = MeshFunction("size_t", mesh, d-1)
boundaries.set_all(0)
left, right, bottom, top = 1, 2, 3, 4
CompiledSubDomain("near(x[0], side) && on_boundary", side = -L).mark(boundaries, left)
CompiledSubDomain("near(x[0], side) && on_boundary", side = 0).mark(boundaries, right)
CompiledSubDomain("near(x[1], side) && on_boundary", side = 0).mark(boundaries, bottom)
CompiledSubDomain("near(x[1], side) && on_boundary", side = L).mark(boundaries, top)
#CompiledSubDomain("near(x[1]*x[1]+x[2]*x[2], R*R) && on_boundary").mark(boundaries, hole)

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
s = Stresses(degree=p)
t_p_left = as_vector((-s[0], -s[2]))
t_p_top = as_vector((s[2], s[1]))

# Dirichlet boundary conditions
bcs = [DirichletBC(V.sub(0), Constant((0.0)), boundaries, right),
	   DirichletBC(V.sub(1), Constant((0.0)), boundaries, bottom)
       ]

# Stress tensor (linear isotropic elasticity)
def sigma(u):
    return lmbda*tr(sym(grad(u)))*Identity(d) + 2.0*mu*sym(grad(u))

# Weak form a==l
a = inner(grad(delta_u), sigma(u))*dx
l = dot(b, delta_u)*dx + dot(t_p_left, delta_u)*ds(left) + dot(t_p_top, delta_u)*ds(top)

#self.m = 0.5*inner(grad(u), sigma(u))*dx - dot(b, u)*dx - dot(t_p, u)*ds(right)

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
#	  tol=1.e-1, self.m=self.m)

e = as_matrix(((s[0],s[2]),(s[2],s[1])))-sigma(u)
print(sqrt(abs(assemble(inner(e,e)*dx))))
#errornorm(Expression((("0.0", "0.0"),("0.0", "0.0")), element=V.ufl_element()), sigma(u))

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

