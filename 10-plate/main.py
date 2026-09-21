from fenics import *
from dolfin import *
import numpy as np
import gmsh # Gmsh: geometry + mesh generator
import meshio # meshio: converts gmsh meshes to FEniCS-readable format
import time


"""
----------------------------------Static Analysis of a 2D Isotropic Plate with an Elliptical Hole----------------------------------

Problem Overview:
This script performs a 2D static linear elasticity analysis of a rectangular isotropic plate containing a central elliptical hole. The main objective is to compute the displacement and stress fields under mechanical loading and validate the numerical FEM solution against an analytical stress field.

Geometry & Mesh:

- 2D Rectangular Plate: L × L domain
- Central Elliptical Hole: semi-axes A and B
- Mesh generated using Gmsh with uniform refinement
- Mesh converted to FEniCS format using meshio (XDMF)

Material Model:

- Linear isotropic elasticity (Hooke’s law)
- Young’s Modulus: E = 200 GPa
- Poisson’s Ratio: ν = 0.25
- Plane stress assumption
- Lamé parameters λ and μ computed from E and ν


Boundary Conditions (BCs):

 Dirichlet BCs:
    - Right boundary: prescribed displacement to model symmetry boundary conditions
    - Bottom boundary: prescribed displacement to model symmetry boundary conditions
 Neumann BCs:
    - Left and top boundaries: traction applied using analytical stress field
 Elliptical hole boundary: traction-free condition

Analytical Solution:

Kirsch-type analytical stress solution implemented in compiled C++ expression
Used to:
- Define consistent boundary tractions
- Compute stress error norm
- Validate FEM accuracy

Mesh Strategy:
- Mesh refinement loop with decreasing element size
- Used for convergence study and error analysis

Post-processing:
- Displacement field exported for visualization
- Stress field projected and saved in MPa
- L2 stress error computed between FEM and analytical solution
- Convergence results stored in text file

Main Objective:
- Study stress concentration around elliptical holes
- Validate FEM solution against analytical benchmark
- Analyze convergence behavior under mesh refinement

Main Learnings:
- Integrating Gmsh + meshio workflow to generate complex geometries (plate with elliptical hole) and import them into FEniCS
- Using physical groups in Gmsh to correctly define and transfer boundary markers for applying boundary conditions
- Embedding analytical solutions via compiled C++ expressions (pybind11) for efficient evaluation inside FEM
- Performing coordinate transformations (Cartesian → elliptical → polar) to evaluate exact stress solutions
- Performing mesh refinement studies to investigate convergence behavior

---------------------------------------------------------------------------------------------------------------------------"""
def tic(): # store start time globally
    global _start_time
    _start_time = time.time()

def toc(): # compute elapsed time since last tic()
    elapsed = time.time() - _start_time
    print("Elapsed time:", elapsed)
    return elapsed

# -------------------------
# C++ compiled stress field
# -------------------------

cpp_code = '''
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <dolfin/function/Expression.h>
#include <cmath>

namespace py = pybind11;
using namespace dolfin;

class Stress : public Expression  //This defines a FEM field (stress tensor) that can be evaluated at any point in the domain. We will use this to compute the analytical stress field for error calculation and to apply boundary tractions.
{
public:
    double a, b, beta, t; // parameters

    Stress() : Expression(2,2), a(1.0), b(1.0), beta(0.0), t(1.0) {} // creates a 2 x 2 tensor-valued expression and initializes parameters

    void eval(Eigen::Ref<Eigen::VectorXd> values,
              Eigen::Ref<const Eigen::VectorXd> x) const override // this function is called by FEniCS to evaluate the expression at a given point x. The computed stress tensor components are stored in values.
    {
        const double m = (a-b)/(a+b); // shape parameter for ellipse

        double signx = (x[0] > 0) ? 1 : ((x[0] < 0) ? -1 : 0); // sign function for x-coordinate
        double signy = (x[1] > 0) ? 1 : ((x[1] < 0) ? -1 : 0); // sign function for y-coordinate

        double ure = (x[0]*x[0]-x[1]*x[1])-(a*a-b*b);  // transforms from Cartesian to elliptical coordinates
        double urho = sqrt(ure*ure+4.0*x[0]*x[0]*x[1]*x[1]); 

        // normalized elliptical coordinates (zeta = re/(a+b), zetaim = rho/(a+b))
        double zetare = (x[0]+signx*sqrt(0.5*std::max(urho+ure, 0.0)))/(a+b);
        double zetaim = (x[1]+signy*sqrt(0.5*std::max(urho-ure, 0.0)))/(a+b);

        // convert to polar coordinates (rho, theta) for use in stress equations
        double rho = sqrt(zetare*zetare+zetaim*zetaim);
        double theta = atan2(zetaim, zetare);

        // used in analytical Kirsch stress solution
        double rho2 = rho*rho;
        double rho4 = rho2*rho2;
        double rho6 = rho4*rho2;

        // represent complex stress function solution around hole
        double S1 = (rho4-2.0*rho2*cos(2.0*theta-2.0*beta)+2.0*m*cos(2.0*beta)-m*m)/(rho4-2.0*m*rho2*cos(2.0*theta)+m*m);
	    
        double S2zr = (rho6*cos(-6.0*theta+2.0*beta) - 2.0*m*rho2*(cos(-4.0*theta-2.0*beta) - m*cos(-4.0*theta)) - 3.0*m*rho4*cos(-4.0*theta+2.0*beta) - (m*m-2.0*cos(2.0*beta)*m+1.0)*rho4*cos(4.0*theta)	- 2.0*rho4*(cos(-2.0*theta-2.0*beta) - m*cos(-2.0*theta)) + 3.0*rho2*cos(2.0*theta+2.0*beta) - m*(m*m-2.0*cos(2.0*beta)*m+1.0)*rho2*cos(2.0*theta) - m*cos(2.0*beta));
        double S2zi = (rho6*sin(-6.0*theta+2.0*beta) - 2.0*m*rho2*(sin(-4.0*theta-2.0*beta) - m*sin(-4.0*theta)) - 3.0*m*rho4*sin(-4.0*theta+2.0*beta) - (m*m-2.0*cos(2.0*beta)*m+1.0)*rho4*sin(-4.0*theta) - 2.0*rho4*(sin(-2.0*theta-2.0*beta) - m*sin(-2.0*theta)) + 3.0*rho2*sin(-2.0*theta-2.0*beta) - m*(m*m-2.0*cos(2.0*beta)*m+1.0)*rho2*sin(-2.0*theta) - m*sin(-2.0*beta));
        
        double S2nr = rho6*cos(6.0*theta)-3.0*m*rho4*cos(4.0*theta)+3.0*m*m*rho2*cos(2.0*theta)-m*m*m;
        double S2ni = -rho6*sin(6.0*theta)+3.0*m*rho4*sin(4.0*theta)-3.0*m*m*rho2*sin(2.0*theta);
        
        // convert complex stress function to Cartesian components
        double ReS2 = (S2zr*S2nr+S2zi*S2ni)/(S2nr*S2nr+S2ni*S2ni);
        double ImS2 = (S2zi*S2nr-S2zr*S2ni)/(S2nr*S2nr+S2ni*S2ni);

        // actual stress components in Cartesian coordinates (Kirsch solution)
        double sigmaxx = t*0.5*(S1+ReS2);
        double sigmayy = t*0.5*(S1-ReS2);
        double sigmaxy = t*0.5*(ImS2);

        // rotate stress tensor by angle beta to align with global coordinates
        //sigma_xx
        values[0] = cos(beta)*cos(beta)*sigmaxx+sin(beta)*sin(beta)*sigmayy+2.0*cos(beta)*sin(beta)*sigmaxy;
        //sigma_yy
        values[3] = sin(beta)*sin(beta)*sigmaxx+cos(beta)*cos(beta)*sigmayy-2.0*cos(beta)*sin(beta)*sigmaxy;
        //sigma_xy
        values[1] = -cos(beta)*sin(beta)*sigmaxx+sin(beta)*cos(beta)*sigmayy+(cos(beta)*cos(beta)-sin(beta)*sin(beta))*sigmaxy;
        //sigma_yx
        values[2] = values[1];
    }
};

PYBIND11_MODULE(SIGNATURE, m) //creates a Python module, SIGNATUre -> module name (FEniCS replaces this internally), m -> the module object you add things to
{
    py::class_<Stress, std::shared_ptr<Stress>, Expression>(m, "Stress")
        .def(py::init<>()) //Allow Python to call Stress() (default constructor)
        .def_readwrite("a", &Stress::a) //exposes C++ variables to Python
        .def_readwrite("b", &Stress::b)
        .def_readwrite("beta", &Stress::beta)
        .def_readwrite("t", &Stress::t);
}
'''

#---------------------------------------------------------------------------------------------------------
# Geometry and mesh generation
#---------------------------------------------------------------------------------------------------------

d = 2      # spatial dimension (2D)
p = 1      # polynomial degree of FEM

A = 1      # ellipse semi-axis in x
B = 1      # ellipse semi-axis in y
L = 4      # plate half-size
base_lc = 0.5  # base mesh element size

# -----------------------------------------------------------------------------------------------------------------
# Initialize Gmsh
# ----------------------------------------------------------------------------------------------------------

gmsh.initialize() # start gmsh kernel
gmsh.model.add("ellipse_refine") # create new model

# Geometry 
rect    = gmsh.model.occ.addRectangle(-L, 0, 0, L, L)
ellipse = gmsh.model.occ.addDisk(0, 0, 0, A, B)

# Cut hole
result, _ = gmsh.model.occ.cut([(2, rect)], [(2, ellipse)])
gmsh.model.occ.synchronize()

#  Define physical domain (solid region)
domain = [result[0][1]]
gmsh.model.addPhysicalGroup(2, domain, 1)

boundary = gmsh.model.getBoundary(result, oriented=False) # Extract boundary edges of geometry

#  containers for boundary tags
left   = [] 
right  = []
bottom = []
top    = []

# boundary marker IDs
LEFT, RIGHT, BOTTOM, TOP = 1, 2, 3, 4

# classify boundary edges by centroid location
for dim, tag in boundary:
    com = gmsh.model.occ.getCenterOfMass(dim, tag) # compute centroid of boundary edge; this is used to classify which edge belongs to which side of the plate for applying boundary conditions
    if np.isclose(com[0], -L): # if x-coordinate of centroid is close to -L, it's a left edge
        left.append(tag) # append the edge tag to the left boundary group
    elif np.isclose(com[0], 0):
        right.append(tag)
    elif np.isclose(com[1], 0):
        bottom.append(tag)
    elif np.isclose(com[1], L):
        top.append(tag)

# assign physical groups for boundary conditions
gmsh.model.addPhysicalGroup(1, left,   LEFT) # assign left edges to physical group LEFT, which we will use to apply Dirichlet BCs in FEniCS; the number 1 indicates that these are 1D entities (lines) in a 2D geometry; the variable LEFT is just an integer tag for this group
gmsh.model.addPhysicalGroup(1, right,  RIGHT)
gmsh.model.addPhysicalGroup(1, bottom, BOTTOM)
gmsh.model.addPhysicalGroup(1, top,    TOP)

# Force MSH2 ASCII format — meshio reads this reliably at all mesh sizes
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2) # Gmsh 4.x defaults to MSH4 format which meshio cannot read; this forces it back to MSH2
gmsh.option.setNumber("Mesh.Binary", 0) # Gmsh 4.x defaults to binary MSH2 which meshio cannot read; this forces it to ASCII MSH2

#---------------------------------------------------------------------------------------------------------
# Output files
#---------------------------------------------------------------------------------------------------------

file_u = File("displacement_in_meters.pvd", "compressed") #displacement output file
file_s = File("stress_in_MPa.pvd",       "compressed") #stress output file

file_h = open("convergence.txt", "w") #convergence data output file
file_h.write("n_cells\th_max\terror\tt_mesh\tt_solve\tt_error\n")

prev_error = None # store previous error for convergence ratio calculation

#---------------------------------------------------------------------------------------------------------
# Refinement loop
#---------------------------------------------------------------------------------------------------------

for i in range(2, 6):   # i=0 is the coarsest mesh, we start from i=2 

    print("\n========== MESH {} ==========".format(i))
    print("Generate: ", end="", flush=True)
    tic() # start timer for mesh generation

    # Refine mesh by halving element size in each iteration
    lc = base_lc / (2**i)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", lc) # set minimum element size in Gmsh; this controls how fine the mesh is, especially around the hole where we expect high stress gradients; halving lc in each iteration creates a finer mesh for convergence testing
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", lc) # set Maximum element size to the same value to create a uniform mesh; in practice, you might want to allow larger elements away from the hole for efficiency, but here we keep it uniform for simplicity 
    gmsh.model.mesh.clear() # clear previous mesh (keep geometry)
    gmsh.model.mesh.generate(2) # generate 2D mesh
    gmsh.write("geom.msh") # write mesh to file (Gmsh format)

    msh = meshio.read("geom.msh") # read mesh using meshio (creates a mesh object with points, cells, and cell data)

     # extract triangle elements (domain) and their physical tags
    triangle_cells = msh.get_cells_type("triangle") 
    triangle_data  = msh.get_cell_data("gmsh:physical", "triangle")

    # extract boundary line elements and their physical tags
    line_cells = msh.get_cells_type("line")
    line_data  = msh.get_cell_data("gmsh:physical", "line")

    # Strip z-coordinate to keep mesh genuinely 2D
    points_2d = msh.points[:, :2]

    # write volume mesh (triangles) and boundary mesh (lines) in XDMF format for FEniCS
    meshio.write("mesh.xdmf",
        meshio.Mesh(points=points_2d,
                    cells={"triangle": triangle_cells},
                    cell_data={"name_to_read": [triangle_data]}))

    # write boundary mesh (lines) and their physical tags in XDMF format for FEniCS
    meshio.write("mf.xdmf",
        meshio.Mesh(points=points_2d,
                    cells={"line": line_cells},
                    cell_data={"name_to_read": [line_data]}))

#---------------------------------------------------------------------------------------------------------
# Load mesh into FEniCS
#---------------------------------------------------------------------------------------------------------

    mesh = Mesh() # create empty FEniCS mesh object
    with XDMFFile("mesh.xdmf") as f: # read mesh geometry and topology from XDMF file created by meshio; this populates the mesh object with vertices and cells (triangles) that define the computational domain
        f.read(mesh)

    mvc = MeshValueCollection("size_t", mesh, 1) # create a MeshValueCollection to hold the boundary markers; "size_t" indicates that the markers are integers, mesh is the mesh we just loaded, and 1 indicates that these are defined on facets (edges in 2D)
    with XDMFFile("mf.xdmf") as f: # read boundary markers from XDMF file created by meshio; this populates the MeshValueCollection with the physical tags for each boundary edge, which we will use to apply boundary conditions and integrate over boundaries in FEniCS
        f.read(mvc, "name_to_read")

    clock_mesh = toc()
    print("n_elems={}, h_max={:.6f} ({:.3f}s)".format(
        mesh.num_cells(), mesh.hmax(), clock_mesh), flush=True)

    # convert mesh markers into FEniCS boundary function for applying boundary conditions and integrating over boundaries
    boundaries = MeshFunction("size_t", mesh, mvc)
    print("Unique boundary markers:", np.unique(boundaries.array()))

    # define measure for integrating over boundaries, using the boundary markers from Gmsh
    ds = Measure("ds", domain=mesh, subdomain_data=boundaries)
    
#-------------------------------------------------------------------------------------------------------
# Function spaces
#-------------------------------------------------------------------------------------------------------

    V       = VectorFunctionSpace(mesh, "Lagrange", p) # vector function space for displacement field
    u       = TrialFunction(V) # trial function (unknown solution we want to solve for)
    delta_u = TestFunction(V) # test function (used in weak form; think of it as a "virtual displacement" in mechanics)

#--------------------------------------------------------------------------------------------------------
# Material properties
#--------------------------------------------------------------------------------------------------------

    E  = 200.e9 # Young's modulus in Pascals (200 GPa)
    nu = 0.25 # Poisson's ratio (dimensionless)
    mu    = E / (2.0*(1.0 + nu))           # plane stress
    lmbda = E*nu / ((1.0 + nu)*(1.0 - nu)) # plane stress

    b = Constant((0.0, 0.0)) # body force (zero in this case since we're only applying boundary tractions)

     # analytical stress field from compiled C++ expression
    stress = CompiledExpression(compile_cpp_code(cpp_code).Stress(), degree=p+2) # create an instance of the compiled C++ expression class Stress and wrap it as a FEniCS expression; degree=p+2 indicates that we want to use a higher degree for evaluating this expression to capture the stress gradients accurately, especially near the hole
    stress.a, stress.b, stress.beta, stress.t = float(A), float(B), 0.0, 1.0 

    # traction vectors on boundaries
    t_p_left   = as_vector((-stress[0,0], -stress[0,1])) 
    t_p_top    = as_vector(( stress[1,0],  stress[1,1]))
    
    # Dirichlet boundary conditions
    bcs = [DirichletBC(V.sub(0), Constant((0.0)), boundaries, RIGHT), # fix x-displacement on right edge to prevent rigid body motion in x
	       DirichletBC(V.sub(1), Constant((0.0)), boundaries, BOTTOM) # fix y-displacement on bottom edge to prevent rigid body motion in y
           ]

#---------------------------------------------------------------------------------------------------------
# Variational formulation (weak form)
#---------------------------------------------------------------------------------------------------------
    # Strain tensor
    def epsilon(u):
        return sym(grad(u))
    
    def sigma(u): # stress tensor
        return lmbda*tr(epsilon(u))*Identity(2) + 2.0*mu*epsilon(u)

    

     # weak form of linear elasticity: integrate the dot product of the stress tensor (sigma) and the symmetric gradient of the test function (delta_u) over the domain; this represents the internal virtual work; on the right-hand side, we have contributions from body forces and boundary tractions
    a = inner(sigma(u), grad(delta_u)) * dx

    l = (dot(b,          delta_u)*dx # body force contribution (zero in this case)
       + dot(t_p_left,   delta_u)*ds(LEFT) # traction contribution on left edge (note the negative sign in t_p_left to convert from stress to traction)
       + dot(t_p_top,    delta_u)*ds(TOP)) # traction contribution on top edge (note the positive sign in t_p_top to convert from stress to traction); no tractions on right and bottom edges since they are Dirichlet boundaries

#---------------------------------------------------------------------------------------------------------
# Solve system
#---------------------------------------------------------------------------------------------------------

    print("Solve... ", end="", flush=True)
    tic()
    u_sol = Function(V) # displacement solution function to be solved for
    solve(a == l, u_sol, bcs=bcs,
          solver_parameters={"linear_solver": "mumps"},
          form_compiler_parameters={"optimize": True})
    clock_solve = toc()
    print("({:.3f}s)".format(clock_solve), flush=True)

#---------------------------------------------------------------------------------------------------------
# Error
#---------------------------------------------------------------------------------------------------------

    print("Error: ", end="", flush=True)
    tic()
    e     = stress - sigma(u_sol)  # stress error tensor
    error = sqrt(abs(assemble(inner(e, e)*dx))) # compute L2 norm of stress error over the domain
    clock_error = toc()
    print("{:.6f}  ({:.3f}s)".format(error, clock_error), flush=True)

    prev_error = error

#---------------------------------------------------------------------------------------------------------
# Post processing
#---------------------------------------------------------------------------------------------------------

    u_sol.rename("u", "displacement")
    file_u << u_sol # write displacement solution to file for visualization in Paraview

    S_scalar  = FunctionSpace(mesh, "Lagrange", p) # scalar function space for projecting stress components

    stress_xx = project(sigma(u_sol)[0, 0], S_scalar, solver_type="mumps")
    stress_xx.rename("sigma_xx", "stress")
    file_s << stress_xx

    file_h.write("{}\t{:.8e}\t{:.8e}\t{:.4f}\t{:.4f}\t{:.4f}\n".format(
        mesh.num_cells(), mesh.hmax(), error,
        clock_mesh, clock_solve, clock_error))
# ---------------------------------------------------------------------------------------------------------------
# Cleanup
# ---------------------------------------------------------------------------------------------------------------

file_h.close() # close convergence data file
gmsh.finalize() # cleanly shut down gmsh kernel and free resources
#--------------------------------------------------------------------------------------------------------------------------------
