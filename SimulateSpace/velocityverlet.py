import numpy as np
#%matplotlib inline 
import matplotlib.pyplot as plt

m_earth = 3.00346e-6  # in solar masses
M_sun = 1.
G = 4*np.pi**2

# reduced mass
mu = m_earth * M_sun / (m_earth + M_sun)

v0y = 6.179 # estim: 2*np.pi       # initial velocity AU/yr
r0x = 1.0           # AU

r0 = np.array([r0x, 0])
v0 = np.array([0, v0y])

dt = 1e-2   # in years (try a coarse time step)
t_max = 10. # year
N_steps = int(t_max / dt)

r = np.zeros(2)  # (x, y)
v = np.zeros(2)  # (vx, vy)
t = 0

def F_gravity(r, m, M):
    rlength = np.sqrt(np.sum(r*r))
    return -G*m*M/rlength**3 * r

# initial conditions
r[:] = r0
v[:] = v0

r_values = [r0]
v_values = [v0]
t_values = [t]

# need to start velocity Verlet with current acceleration
a = F_gravity(r, m_earth, M_sun)/mu

for i in range(N_steps):
    t += dt
    
    # velocity Verlet
    v_half = v + 0.5*dt*a        # half-step velocity
    r_new = r + dt*v_half
    a_new = F_gravity(r_new, m_earth, M_sun)/mu
    v_new = v_half + 0.5*dt*a_new
    
    r_values.append(r_new)
    v_values.append(v_new)
    t_values.append(t)

    r[:] = r_new
    v[:] = v_new
    a[:] = a_new    # important: use force/acceleration for next step!

# turn lists of results into arrays for easier processing
r_values = np.array(r_values)
v_values = np.array(v_values)
t_values = np.array(t_values)

ax = plt.subplot(1,1,1)
ax.plot(r_values[:, 0], r_values[:, 1])
ax.set_aspect(1)