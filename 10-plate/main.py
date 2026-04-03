from ast import Constant

from fenics import *
from mshr import *
import numpy as np

# Set log level to only show warnings and errors (suppress info messages)
set_log_level(WARNING)

"""-------------------------2D Error Convergence Study: Stress Concentration Around an Elliptical Hole-----------------------------------------------------------------------------------------

Problem Overview:
This script performs a 2D plane-stress analysis of a rectangular plate with an 
elliptical cutout. It uses a C++ Expression to implement the exact analytical 
solution for stress, allowing for a rigorous L2-norm error convergence study 
across multiple mesh densities.

Geometry:
- Domain: A quadrant (Rectangle - Ellipse) representing a symmetric portion of a plate.
- Hole Geometry: Elliptical cutout at the origin (Semi-major A, Semi-minor B).
- Refinement: An automated loop generates increasingly finer meshes (from 2^2 to 2^9 
  cells) using the 'mshr' generator to track error reduction.

Material Model:
- Linear Isotropic Elasticity (Plane Stress formulation).
- Young's Modulus (E): 200 GPa; Poisson's Ratio (nu): 0.25.
- Exact Stress Solution: Provided via a custom C++ class (Stress) linked to FEniCS.

Boundary Conditions:
- Dirichlet BCs: Symmetry conditions (fixed u_x on right, fixed u_y on bottom).
- Neumann BCs: Analytical tractions derived from the exact stress solution are 
  applied to the left and top boundaries to simulate an infinite plate response.

Main Learnings:
- Integrating complex C++ analytical expressions into FEniCS for validation.
- Performing automated L2-norm error analysis (Analytical vs. FEM).
- Benchmarking performance (timers for mesh generation, solving, and error calculation).
------------------------------------------------------------------------------------------------------------------------------"""
# c++ code for the analytical stress 
cpp_code = ''' 
class Stress : public Expression {

public:
    // Constructor for the Stress class, which initializes the parameters a, b, beta, and t with default values.
    Stress() : Expression(2, 2), a(1.0), b(1.0), beta(0.0), t(1.0) {}
        
    void eval(Array<double> &values, const Array<double> &x) const {
        
        // Geometry parameter for elliptical hole
        const double m = (a-b)/(a+b);
        // Sign of coordinates (used for branch selection)
        double signx = (x[0] > 0) ? 1 : ((x[0] < 0) ? -1 : 0); 
        double signy = (x[1] > 0) ? 1 : ((x[1] < 0) ? -1 : 0); 
        
        // --- Conformal Mapping: Ellipse to Circle ---
        // Transformation for elliptical coordinates

        // These variables perform the mapping from Cartesian (x,y) to the complex plane (zeta) where the elliptical hole becomes a unit circle.
        double ure = (x[0]*x[0]-x[1]*x[1])-(a*a-b*b); // Real part of the complex variable in elliptical coordinates
        double urho = sqrt(ure*ure+4.0*x[0]*x[0]*x[1]*x[1]); // Magnitude of the complex variable in elliptical coordinates

        // zetare and zetaim represent the Real and Imaginary parts of the mapped coordinate zeta in the complex plane, which corresponds to the elliptical coordinates of the point (x[0], x[1]).
        double zetare = (x[0]+signx*sqrt(0.5*max(urho+ure, 0.0)))/(a+b); // Real part of the transformed coordinate in elliptical coordinates
        double zetaim = (x[1]+signy*sqrt(0.5*max(urho-ure, 0.0)))/(a+b); // Imaginary part of the transformed coordinate in elliptical coordinates
        
        // Polar coordinates
        // Convert the mapped coordinate to Polar form (rho, theta
        double rho = sqrt(zetare*zetare+zetaim*zetaim); # Radial coordinate 
        double theta = atan2(zetaim, zetare); # Angular coordinate 
        
        
        // Powers of rho for stress computation
        double rho2 = rho*rho;
        double rho4 = rho2*rho2;
        double rho6 = rho4*rho2;
        
        // --- Compute stress function terms ---

        // S1 represents the sum of principal stresses (sigma_xx + sigma_yy)
        // This is derived from the first Kolosov-Muskhelishvili potential.
        double S1 = (rho4-2.0*rho2*cos(2.0*theta-2.0*beta)+2.0*m*cos(2.0*beta)-m*m)/(rho4-2.0*m*rho2*cos(2.0*theta)+m*m); 
	    
        // Complex S2 components (real and imaginary)
        // S2zr/S2zi and S2nr/S2ni represent the Numerator and Denominator of the complex potential for the stress difference (sigma_yy - sigma_xx + 2i*sigma_xy)
        double S2zr = (rho6*cos(-6.0*theta+2.0*beta) - 2.0*m*rho2*(cos(-4.0*theta-2.0*beta) - m*cos(-4.0*theta)) - 3.0*m*rho4*cos(-4.0*theta+2.0*beta) - (m*m-2.0*cos(2.0*beta)*m+1.0)*rho4*cos(4.0*theta)	- 2.0*rho4*(cos(-2.0*theta-2.0*beta) - m*cos(-2.0*theta)) + 3.0*rho2*cos(2.0*theta+2.0*beta) - m*(m*m-2.0*cos(2.0*beta)*m+1.0)*rho2*cos(2.0*theta) - m*cos(2.0*beta));
        double S2zi = (rho6*sin(-6.0*theta+2.0*beta) - 2.0*m*rho2*(sin(-4.0*theta-2.0*beta) - m*sin(-4.0*theta)) - 3.0*m*rho4*sin(-4.0*theta+2.0*beta) - (m*m-2.0*cos(2.0*beta)*m+1.0)*rho4*sin(-4.0*theta) - 2.0*rho4*(sin(-2.0*theta-2.0*beta) - m*sin(-2.0*theta)) + 3.0*rho2*sin(-2.0*theta-2.0*beta) - m*(m*m-2.0*cos(2.0*beta)*m+1.0)*rho2*sin(-2.0*theta) - m*sin(-2.0*beta));
        
        double S2nr = rho6*cos(6.0*theta)-3.0*m*rho4*cos(4.0*theta)+3.0*m*m*rho2*cos(2.0*theta)-m*m*m;
        double S2ni = -rho6*sin(6.0*theta)+3.0*m*rho4*sin(4.0*theta)-3.0*m*m*rho2*sin(2.0*theta);
        
        // Perform complex division to get Real (ReS2) and Imaginary (ImS2) parts
        double ReS2 = (S2zr*S2nr+S2zi*S2ni)/(S2nr*S2nr+S2ni*S2ni);
        double ImS2 = (S2zi*S2nr-S2zr*S2ni)/(S2nr*S2nr+S2ni*S2ni);
        
         // Compute final stress components in local axes
        double sigmaxx = t*0.5*(S1+ReS2); // t is the far-field tension magnitude
        double sigmayy = t*0.5*(S1-ReS2);
        double sigmaxy = t*0.5*(ImS2);
        

        // --- Tensor Rotation ---
        // If the far-field stress is applied at an angle 'beta', we rotate the stress tensor back to the global coordinate system.
        
        //sigma_xx (Row 0, Col 0)
        values[0] = cos(beta)*cos(beta)*sigmaxx+sin(beta)*sin(beta)*sigmayy+2.0*cos(beta)*sin(beta)*sigmaxy;
        
        //sigma_yy (Row 1, Col 1
        values[3] = sin(beta)*sin(beta)*sigmaxx+cos(beta)*cos(beta)*sigmayy-2.0*cos(beta)*sin(beta)*sigmaxy;

        //sigma_xy (Row 0, Col 1)
        values[1] = -cos(beta)*sin(beta)*sigmaxx+sin(beta)*cos(beta)*sigmayy+(cos(beta)*cos(beta)-sin(beta)*sin(beta))*sigmaxy;

        // sigma_yx (Row 1, Col 0) - for symmetric stress tensor, sigma_yx = sigma_xy
        values[2] = values[1];
    }

public:
    double a, b, beta, t;  // ellipse geometry & stress scaling
}; 
'''

#---------------------------------------------------------------------------------------------------------
# Geometry and mesh generation
#---------------------------------------------------------------------------------------------------------


dim = 2 # spatial dimension
p = 1 # polynomial degree for finite element

# Parameters
A = 1 # semi-major axis of the elliptical hole
B = 1 # semi-minor axis of the elliptical hole
L = 4 # length of the plate in x and y directions

# Create geometry
plate = Rectangle(Point(-L, 0), Point(0, L)) # rectangular plate 
hole = Ellipse(Point(0, 0), A, B) # elliptical area at the origin
geometry = plate - hole # subtract elliptical area from plate to create the final geometry

#create output files
file_u = File("displacement_in_meters.pvd", "compressed") 
file_s = File("stress_in_MPa.pvd", "compressed")

# Stores convergence data.
file_h = open("history{}.dat".format(p), "w")
file_h.write("elems\thmax\terror\tclock_mesh\tclock_solve\tclock_error\n") # header for convergence data file

# mesh refinement loop
for i in range(2, 10):
    
    print("MESH", i)
    print("Generate: ", end="")
    # Create mesh
    tic() # start timer for mesh generation
    mesh = generate_mesh(geometry, 2**i) # generate mesh with 2^i cells in the longest direction (x or y) 
    clock_mesh = toc() # stop timer for mesh generation and store elapsed time
    print("n_elems={}, h_max={} ({})".format(mesh.num_cells(), mesh.hmax(), clock_mesh), flush=True) # print number of elements, maximum cell size, and time taken for mesh generation
    
    #mesh = refine(mesh) # uniform mesh refinement

    #margin = 0.5 # initial margin for adaptive mesh refinement around the hole
    
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

#---------------------------------------------------------------------------------------------------------
# Boundary identification and marking
#---------------------------------------------------------------------------------------------------------

    # Define boundaries
    boundaries = MeshFunction("size_t", mesh, dim-1) # MeshFunction to store boundary IDs on facets
    boundaries.set_all(0) # initialize all boundaries to 0
    
    left, right, bottom, top = 1, 2, 3, 4 # boundary IDs
    CompiledSubDomain("near(x[0], side) && on_boundary", side = -L).mark(boundaries, left) # mark left boundary with ID 1
    CompiledSubDomain("near(x[0], side) && on_boundary", side = 0).mark(boundaries, right) # mark right boundary with ID 2
    CompiledSubDomain("near(x[1], side) && on_boundary", side = 0).mark(boundaries, bottom) # mark bottom boundary with ID 3
    CompiledSubDomain("near(x[1], side) && on_boundary", side = L).mark(boundaries, top) #  mark top boundary with ID 4

    # Surface integral element
    ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

#-------------------------------------------------------------------------------------------------------
# Function spaces and variational functions
#-------------------------------------------------------------------------------------------------------

    # Define function space
    V = VectorFunctionSpace(mesh, "Lagrange", p) # vector function space for displacement field

    # Define trial and test functions
    u = TrialFunction(V) # trial function for unknown displacement field
    delta_u = TestFunction(V) # test function for virtual displacement field

#--------------------------------------------------------------------------------------------------------
# Material properties
#--------------------------------------------------------------------------------------------------------

    E = 200.e9 # Young's modulus in Pa
    nu = 0.25 # Poisson's ratio
    #mu    = E/(2.0*(1.0 + nu))
    #lmbda = E*nu/((1.0 + nu)*(1.0 - 2.0*nu))
    mu    = E/(2.0*(1.0 + nu)) #PLANE STRESS # Shear modulus # Lamé's second parameter
    lmbda = E*nu/((1.0 + nu)*(1.0 - nu)) #PLANE STRESS # Lamé's first parameter

    #rho = 8.e3 # density in kg/m^3
    #g = 9.81 # gravitational acceleration in m/s^2

    # Volume force/ heat source and prescribed tractions/ prescribed heat fluxes
    b = Constant((0.0, 0.0)) # body force in N/m^3
    
    #stress = Stress(A, B, degree=p+2)

    stress = Expression(cpp_code, degree=p+2) # Analytical stress expression
    # This tells FEniCS that we want to create a user-defined function (in this case, the analytical stress) to evaluate at any point in the domain.
    # cpp_code is the C++ code defining the class Stress we wrote earlier, which implements the exact stress solution around the elliptical hole.
    # degree=p+2 specifies the polynomial degree used for interpolation in FEniCS. Even though this is a C++ expression, 
    # FEniCS needs to know a degree for numerical integration inside the FEM assembly. Using p+2 ensures enough accuracy when comparing to the FEM solution

    stress.a, stress.b = 1.0, 1.0  # Parameters for the analytical stress solution (a and b are the semi-major and semi-minor axes of the elliptical hole)
    t_p_left = as_vector((-stress[0,0], -stress[0,1])) # Prescribed traction on the left boundary in N/m^2 (Pa) 
    t_p_top = as_vector((stress[1,0], stress[1,1])) # Prescribed traction on the top boundary in N/m^2 (Pa) 

    # Dirichlet boundary conditions
    bcs = [DirichletBC(V.sub(0), Constant((0.0)), boundaries, right), # prescribe zero displacement in x-direction on the right boundary
	       DirichletBC(V.sub(1), Constant((0.0)), boundaries, bottom) # prescribe zero displacement in y-direction on the bottom boundary
           ]

    #---------------------------------------------------------------------------------------------------------
    #  Variational formulationn (weak form)
    #---------------------------------------------------------------------------------------------------------

    # Strain tensor
    def epsilon(u):
	    return sym(grad(u))
        
    # Stress tensor (linear isotropic elasticity)
    def sigma(u):
        return lmbda*tr(epsilon(u))*Identity(dim) + 2.0*mu*epsilon(u)

    # Weak form a==l
    a = inner(grad(delta_u), sigma(u))*dx # bilinear form representing the internal virtual work (stiffness matrix)
    l = dot(b, delta_u)*dx + dot(t_p_left, delta_u)*ds(left) + dot(t_p_top, delta_u)*ds(top) # linear form representing the external virtual work 


#---------------------------------------------------------------------------------------------------------
# Solve the linear system
#---------------------------------------------------------------------------------------------------------

    print("Solve... ", end="")
    tic() # start timer for linear system solution
    u = Function(V) # function to hold the solution (displacement field) after solving the linear system
    #K = assemble(a)
    #F = assemble(l)
    #for bc in bcs:
    #	bc.apply(K, F)
    #U = u.vector()
    #solve(K, U, F)
    #solver = LUSolver(K, "mumps")
    #solver.solve(U, F)
    # solve Ku = F for u using the assembled stiffness matrix K and load vector F, with the specified boundary conditions applied
    solve(a == l, u, bcs=bcs, 
	      solver_parameters={"linear_solver": "mumps"},
	      form_compiler_parameters={"optimize": True})
    clock_solve = toc() # stop timer for linear system solution and store elapsed time
    print("({})".format(clock_solve), flush=True) # print time taken for linear system solution
    
    print("Error: ", end="")
    tic() # start timer for error calculation
    e = stress-sigma(u) # error in stress field (analytical stress - FEM stress)
    error = sqrt(abs(assemble(inner(e,e)*dx)))  # L2 norm of the error in the stress field, calculated by integrating the inner product of the error with itself over the domain and taking the square root.

    #E = TensorFunctionSpace(mesh, "Discontinuous Lagrange", p+2)
    #stressE = interpolate(stress, E)
    #sigmaE = project(sigma(u), E)
    #eE = Function(E)
    #eE.assign(stressE)
    #eE.vector().axpy(-1.0, sigmaE.vector())

    #tic()
    #stressE = interpolate(stress, E)
    #print("\nintp",toc())
    #tic()
    #sigmaE = Function(E)
    #sigmaE = project(sigma(u), E, solver_type="cg")
    #print("proj",toc())
    #tic()
    #eE = Function(E)
    #eE.assign(stressE)
    #eE.vector().axpy(-1.0, sigmaE.vector())
    #print("assg",toc())
    #tic()
    #error = norm(eE, norm_type="L2", mesh=mesh)
    #print("norm",toc())
    
    clock_error = toc() # stop timer for error calculation and store elapsed time
    print("{} ({})".format(error, clock_error), flush=True) # print error and time taken for error calculation
        
    #errornorm(Expression((("0.0", "0.0"),("0.0", "0.0")), element=V.ufl_element()), sigma(u))

    #---------------------------------------------------------------------------------------------------------
    # Post processing
    #---------------------------------------------------------------------------------------------------------


    # Create displacement file
    u.rename("u", "displacement") # rename displacement for output
    file_u << u # save displacement to file

    # Project stress field and create stress file
    S = FunctionSpace(mesh, "Lagrange", p) # scalar function space for stress component sigma_xx
    stress_proj = project(sigma(u)[0,0], S, solver_type="mumps") #
    stress_proj.rename("sigma_xx", "stress") # rename stress for output

    file_s << stress_proj  # save stres

    # record convergence and timings
    file_h.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(mesh.num_cells(), mesh.hmax(), error, clock_mesh, clock_solve, clock_error))
    
file_h.close() #close the history file
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------