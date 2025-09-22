import numpy as np
import matplotlib.pyplot as plt



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
    "z_0": -0.3, # initial position, m

    "z_1": 0, # lens opening, m

    "L": 0.2, # lens length, m
    "L_alt": 0.75, # alternative length, m

    "R": 0.1, # lens radius, m
    "R_alt": 0.05, # alternative radius, m

    "B_0": 1e-3, # lens magnetic field, T

    "v_0": 0.06 * constant["c"] , # initial velocity, m/s
}

### Function for Finding B_z

def B_z(z, z_1 = lens["z_1"], z_2 = lens["L"], R = lens["R"], B_0=lens["B_0"]):
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
    ### returns tuple:  (B_z, B_z / B_0)

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

## All of the above is from angels code. Can be imported in a separate file
## Function for finding the magnetic force on the particle

def F_B_r(t, y):
    #constants
    q = constant["e"]
    m = constant["m_e"]
    #y = [rx, ry, rz, vx, vy, vz, L, R, B_0]
    r = np.array([y[0], y[1], y[2]])
    v = y[3:]
    B = B_r(r, y[2])[0]
    ## Equations of motion
    a = (q/m)*np.cross(v,B)
    return np.array([a[0], a[1], a[2], a[0], a[1], a[2]])

## rk4 method
def rk4(f, y0, tmax, dt):
    """
    Applies the 4th order Runge-Kutta method to solve an ODE.

    Args:
        f (function): The ODE function, must be in the form f(y, t, *args)
        y0 (float or array): The initial value of y
        t (array): The time values at which to solve for y

    Returns:
        array: The approximate solution y at the given time values t
    """
    t = np.arange(0, tmax+dt, dt)
    y = np.zeros([len(t), len(y0)])
    #positions = [[0, y[0], y[1], y[2]]]
    #velocities = [[0, y[3], y[4], y[5]]]

    for i in range(len(t)-1):
        k1 = dt*f(t[i], y[i])
        k2 = dt*f(t[i] + 0.5*dt, y[i] + 0.5*dt*k1)
        k3 = dt*f(t[i] + 0.5*dt, y[i] + 0.5*dt*k2)
        k4 = dt*f(t[i] + dt, y[i] + k3)
        y[i+1] = y[i] + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
        #positions.append([t, y2[0], y2[1], y2[2]])
        #velocities.append([t, y[3], y[4], y[5]])

    return t, y


## Plots
# Plot x(t), y(t), r(t) (Eq. 4) in one graph, z(t) in a second graph.
# Describe the motions.

#z = np.linspace(-0.5, 1, 1501)
r0_lens = np.array([0.75*lens["R"],0, lens["z_0"]])
v0_lens = np.array([0, 0, lens["v_0"]])
y0_lens = np.array([r0_lens[0], r0_lens[1], r0_lens[2], v0_lens[0], v0_lens[1], v0_lens[2]])

t, y = rk4(f=F_B_r, y0=y0_lens, tmax=1e-7, dt=1e-8)
#t = vals[:,0]
#x_t = vals[:,1]
#y_t = vals[:,2]
#z_t = vals[:,3]
#r_t = np.sqrt(x_t**2 + y_t**2)

#plt.plot(t, x_t)
#plt.xlabel('t [s]')
#plt.ylabel('position [m]')
#plt.legend()
#plt.savefig("positions.png")
# Show the plot
#plt.show()
