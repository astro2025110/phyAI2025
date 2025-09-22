import numpy as np


### Dictionary for Constants in Table 1:

constant = {
    "c": 299792458, # speed of light, m/s
    "e": 1.602176634e-19, # electron charge, A s
    "m_e": 9.1093837015e-31, # electron mass, kg
    "mu_0": 1.25663706212e-6, # magnetic constant, N/(A^2)
    "epsilon_0": 8.8541878128e-12 # electric constant, F/m
}

### Dictionary for Values in Table 2 + alternative values:

lens = {
    "z_0": -0.3, # initial z position, m

    "z_1": 0, # lens opening, m

    "L": 0.2, # lens length, m
    "L_alt": 0.75, # alternative length, m

    "R": 0.1, # lens radius, m
    "R_alt": 0.05, # alternative radius, m

    "B_0": 1e-3, # lens magnetic field, T

    "v_0": 0.06 * constant["c"] , # initial velocity, m/s
}

### Function for Finding B_z

def B_z(z, zL_1 = lens["z_1"], z_2 = lens["L"], R = lens["R"], B_0=lens["B_0"]):
    ### returns tuple:  (B_z, B_z / B_0)
    
    ## Finding delta_z values
    dz_1 = z - z_1
    dz_2 = z - z_2

    ## Finding H (hypothenuse) values
    H_1 = np.sqrt(dz_1 ** 2 + R ** 2)
    H_2 = np.sqrt(dz_2 ** 2 + R ** 2)
    
    ## Equation 3.1 
    B_z = 1/2 * B_0 * (dz_1 / H_1 - dz_2 / H_2)

    ## Bz/B0 ratio
    B_ratio = B_z / B_0

    return B_z, B_ratio

### Function for Finding B_r

def B_r(r, z, z_1 = lens["z_1"], z_2 = lens["L"], R = lens["R"], B_0=lens["B_0"]):
    ### returns tuple:  (B_r, B_z / B_0)

    ## Finding delta_z values
    dz_1 = z - z_1
    dz_2 = z - z_2

    ## Finding H (hypothenuse) values
    H_1 = np.sqrt(dz_1 ** 2 + R ** 2)
    H_2 = np.sqrt(dz_2 ** 2 + R ** 2)

    ## Equation 5.2
    B_r = -B_0 * r * (R ** 2) / 4 * (1/(H_1 ** 3) - 1/(H_2 ** 3))

    ## Br/B0 ratio
    B_ratio = B_r / B_0

    return B_r, B_ratio


#### Function that finds Bx, By, Bz, and B

def B(position, z_1 = lens["z_1"], z_2 = lens["L"], R = lens["R"], B_0=lens["B_0"]):
    """
    Function that returns magnetic field in cartesian components.
    
    -Position argument should be array [x, y, z]
    
    -Uses B_r and B_z functions defined above.
    
    -Defines r for B_r function from x and y arguments.
    
    -Returns tuple (Bx, By, Bz, B)"""

    ## Establishing position components
    x, y, z = position[0], position[1], position[2]
        
    ## Defining r
    r = np.sqrt(x**2 + y**2)
    
    ## Finding Br and its cartesian components from equation 9
    Br = B_r(r, z, z_1, z_2, R, B_0)
    Br = Br[0]

    if r == 0:
        Bx = 0
        By = 0
    else:
        Bx = Br * x / r
        By = Br * y / r

    ## Finding B_z from function
    Bz = B_z(z, z_1, z_2, R, B_0)
    Bz = Bz

    ## Finding B from sum of vectors
    B = np.sqrt(Br ** 2 + Bz **2)

    return np.array([Bx, By, Bz]), B

#### Function that Lorentz Force

def F_lorentz(v, B, q = constant["e"]):
    """
    v and B are 3d arrays
    
    returns 3D array"""

    #Fx = q * ((vy * Bz) - (vz * By))
    #Fy = q * ((vx * Bz) - (vz * Bx))
    #Fz = q * ((vx * By) - (vy * Bx))


    return q * np.cross(v, B)

def a(F_lorentz, m=constant["m_e"]):
    """
    F_lorentz is a 3D array
    
    returns 3D array"""

    return F_lorentz / m




