from fenics import *
from mshr import *
import numpy as np


set_log_level(WARNING)

class Stress(Expression):
    def __init__(self, a, b, **kwargs):
        #Expression.__init__(self)
        self.a = A
        self.b = B
        self.m = (self.a-self.b)/(self.a+self.b)
        self.beta = 0
        self.t = 1
        
    def value_shape(self):
        return (2,2,)
        
    def eval(self, values, x):
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
        
        sigmaxx = self.t*0.5*(S1+ReS2)
        sigmayy = self.t*0.5*(S1-ReS2)
        sigmaxy = self.t*0.5*(ImS2)
        
        #sigma_xx
        values[0] = np.cos(self.beta)*np.cos(self.beta)*sigmaxx+np.sin(self.beta)*np.sin(self.beta)*sigmayy+2.0*np.cos(self.beta)*np.sin(self.beta)*sigmaxy
        #sigma_yy
        values[3] = np.sin(self.beta)*np.sin(self.beta)*sigmaxx+np.cos(self.beta)*np.cos(self.beta)*sigmayy-2.0*np.cos(self.beta)*np.sin(self.beta)*sigmaxy
        #sigma_xy
        values[1] = -np.cos(self.beta)*np.sin(self.beta)*sigmaxx+np.sin(self.beta)*np.cos(self.beta)*sigmayy+(np.cos(self.beta)*np.cos(self.beta)-np.sin(self.beta)*np.sin(self.beta))*sigmaxy
        #sigma_yx
        values[2] = values[1]
    


################################
#### PROBLEM DEFINITION ########
################################

d = 2
p = 3

# Parameters
A = 1
B = 1
L = 4

# Create geometry
plate = Rectangle(Point(-L, 0), Point(0, L))
hole = Ellipse(Point(0, 0), A, B)
geometry = plate - hole

file_u = File("displacement.pvd", "compressed") 
file_s = File("stress.pvd", "compressed")

file_h = open("history{}.dat".format(p), "w")
file_h.write("elems\thmax\terror\tclock\n")

for i in range(1, 7):
    
    # Create mesh
    mesh = generate_mesh(geometry, 2**i)
    print("MESH {0} with n_elems={1} and h_max={2}".format(i, mesh.num_cells(), mesh.hmax()))
    
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

    # Surface integral element
    ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

    # Define function space
    V = VectorFunctionSpace(mesh, "Lagrange", p)

    # Define trial and test functions
    u = TrialFunction(V)
    delta_u = TestFunction(V)

    # Material parameters
    E = 200.e9
    nu = 0.25
    #mu    = E/(2.0*(1.0 + nu))
    #lmbda = E*nu/((1.0 + nu)*(1.0 - 2.0*nu))
    mu    = E/(2.0*(1.0 + nu)) #PLANE STRESS
    lmbda = E*nu/((1.0 + nu)*(1.0 - nu)) #PLANE STRESS

    #rho = 8.e3
    #g = 9.81

    # Volume force/ heat source and prescribed tractions/ prescribed heat fluxes
    b = Constant((0.0, 0.0))
    stress = Stress(A, B, degree=p+3)
    t_p_left = as_vector((-stress[0,0], -stress[0,1]))
    t_p_top = as_vector((stress[1,0], stress[1,1]))

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

    tic()

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

    clock = toc()

    #e = as_matrix(((stress[0],stress[2]),(stress[2],stress[1])))-sigma(u)
    #error = sqrt(abs(assemble(inner(e,e)*dx)))    
    E = TensorFunctionSpace(mesh, "Discontinuous Lagrange", p+3)
    stressE = interpolate(stress, E)
    sigmaE = project(sigma(u), E)
    eE = Function(E)
    eE.assign(stressE)
    eE.vector().axpy(-1.0, sigmaE.vector())
    error = norm(eE, norm_type="L2", mesh=mesh)
    
    
    print("ERROR:", error)
    print("CLOCK:",clock,"\n")
    
    #errornorm(Expression((("0.0", "0.0"),("0.0", "0.0")), element=V.ufl_element()), sigma(u))

    ################################
    #### POST-PROCESSING ###########
    ################################

    # Create displacement and temperature file
    u.rename("u", "displacement")
    file_u << u

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
    stress_proj = project(sigma(u)/1.e6, T, solver_type="mumps")
    stress_proj.rename("sigma", "stress")

    file_s << stress_proj

    
    file_h.write("{}\t{}\t{}\t{}\n".format(mesh.num_cells(), mesh.hmax(), error, clock))
    
file_h.close()
