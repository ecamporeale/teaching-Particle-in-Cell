# Define some quantities

L=2.0 * np.pi # box length
DT=0.5 # timestep
NT=500 # total number of timesteps
NTOUT=25 # save output every other NTOUT
NG=32 # number of gridpoints DECIDE: is this the number of cells or number of edges?
N=10000 # number of particles
WP=1  
QM=- 1 # charge over mass (normalized). For electrons = -1
V0=0.2  # mean initial velocity (one for each species)
T=0.0  # thermal velocity (standard deviation for a Maxwellian distribution)
XP1=0.001 # initial perturbation
# Some more?

# Initialize particles
# position-> uniformly distributed in the box
# velocity -> Normal distribution with std VT and mean V0

# Add a perturbation in the position (small displacement)
mode=1
xp=xp + XP1 * np.sin(2 * np.pi * xp / L * mode)

# check for periodic boundary conditions

# Main Loop

    # Mover --> Use leapfrog to update positions

    # check for periodic boundary conditions

    # Accumulate density from Particles to Grid	(mind the BC)
    # and calculate charge density

    # Solve Poisson's equation with the new charge density	
    # (solve first laplacian(phi) = -rho
    # and then E = - grad(phi)
		
    # Once again: are your boundary conditions satisfied?
    # Quick test: is mean(E)=0	satisfied?

    # Calculate Electric field on the position of the particles
    # (interpolate particle --> grid)	
    
    # Update particle velocity with leapfrog
    vp=vp + QM * Ep * DT
    
    # Put your diagnostics here: calculate energy and momentum 

    # Make plots, movies, etc.



