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
def F_B_rx(v=lens["v_0"], z=1, q=constant["e"]):

    ## Equation 16
    return q*(v[1]*B_z(z) - v[2]*B_r(r, z)[1])

def F_B_ry(v=lens["v_0"], z=1, q=constant["e"]):

    ## Equation 16
    return q*(v[0]*B_z(z) - v[2]*B_r(r, z)[0])

def F_B_rz(v=lens["v_0"], z=1, q=constant["e"]):

    ## Equation 16
    return q*(v[0]*B_r(r, z)[1] - v[1]*B_z(z))
## rk4 method

def rk4(y, f, t, h):
    k1 = f(t, y)
    k2 = f(t+0.5*h, y+0.5*h+k1)
    k3 = f(t+0.5*h, y+0.5*h+k2)
    k4 = f(t+h, y+h+k3)
    return y + (1/6)*h*(k1 + 2*k2 + 2*k3 + k4)

def f_standard(t, y, force, m=1):
    """Force vector in standard ODE form (n=2)

    Arguments
    ---------
    t : float
        time
    y : array
        dependent variables in ODE standard form (2d)
    force : function
        `force(y[0])` returns force (note: will not be
        able to handle velocity dependent forces)
    m : float
        mass
    """
    return np.array([y[1], force(y[0])/m])

def integrate_newton(x0=0, v0=1, t_max=100, h=0.001, mass=1, force=1, integrator=rk4):
    """Integrate Newton's equations of motions.

    Note that all problem parameters such as spring constant k must be
    set consistently in the force function.

    Arguments
    ---------
    x0 : float
       initial position
    v0 : float
       initial velocity
    t_max : float
       time to integrate out to
    h : float (default 0.001)
       integration time step
    mass : float (default 1)
       mass of the particle
    force : function `f(x)`
       function that returns the force when particle is
       at position `x`
    integrator : function `I(y, f, t, h)`
       function that takes the ODE standard form vectors y and f
       together with the current time and the step `h` and returns
       y at time t+h.

    Returns
    -------
    Tuple ``(t, y)`` with times and the ODE standard form vector.
    `y[:, 0]` is position and `y[:, 1]` velocity.

    """

    Nsteps = t_max/h
    t_range = h * np.arange(Nsteps)
    y = np.zeros((len(t_range), 2))

    # initial conditions
    y[0, :] = x0, v0

    # build a function with "our" force
    def f(t, y):
        """ODE force vector"""
        return f_standard(t, y, force, m=mass)

    for i, t in enumerate(t_range[:-1]):
        # copy is necessary to avoid changing the full y trajectory
        y[i+1, :] = integrator(y[i].copy(), f, t, h)

    return t_range, y
## Plots
# Plot x(t), y(t), r(t) (Eq. 4) in one graph, z(t) in a second graph.
# Describe the motions.

z = np.linspace(-0.5, 1, 1501)
r0_lens = np.array([0.75*lens["R"],0, lens["z_0"]])
v0_lens = np.array([0, 0, lens["v_0"]])

x0_lens = r0_lens[0]

t, y = integrate_newton(x0=r0_lens[0], v0=v0_lens[0], t_max=10, h=0.001, mass=constant["m_e"], force=F_B_rx, integrator=rk4)

#plt.plot(times, x_lens)
#plt.xlabel('t [s]')
#plt.ylabel('position [m]')
#plt.legend()
#plt.savefig("positions.png")
# Show the plot
plt.show()
