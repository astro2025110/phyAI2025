'''
This file contains code used for to generate the plots needed for our solar system simulator.

1. Create a function that computes the positions, velocities of the planets and 10 asteroids in the solar system given:
    - The gravity from the Sun and other bodies
    - The radial distance between the object and others
    - The resultant acceleration

'''
# Import all needed packages here
import numpy as np
#%matplotlib inline
#%matplotlib widget
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

# Insert all constants from NASA Horizons data here in a dictionary
const = {
    "G": 6.6743e-11, # m^3 / kg*yr^2
    "AU": 149597870700.0, #m
}

convert = {
    "years to days": 365.25, #average number of days in an Earth year
}
# k = GM, F_grav = -k/r^2
mass = {
    
    "Sun": 1.99e+30,

    #Planets 
    "Mercury": 0.330103e+24, "Venus": 4.86731e+24, 
    "Earth": 5.97217e+24, "Mars": 0.641691e+24,
    "Jupiter": 1898.125e+24, "Saturn": 568.317e+24, 
    "Uranus": 86.8099e+24, "Neptune": 102.4092e+24,

    #Asteroids
    "Ceres": 938.416e+18, "Vesta": 259.076e+18, 
    "Pallas": 204e+18, "Hygiea": 87e+18, 
    "Interamnia": 35e+18, "Eunomia": 30e+18,
    "Juno": 27e+18, "Davida": 27e+18,
    "Europa": 24e+18, "Psyche": 23e+18,
}

r_j0 = {
    
    "Sun": 0,

    #Planets 
    "Mercury": 0.38709843*const["AU"],
    "Venus": 0.72332102*const["AU"], 
    "Earth": const["AU"],
    "Mars": 1.52371243*const["AU"],
    "Jupiter": 5.20248019*const["AU"], 
    "Saturn": 9.54149883*const["AU"], 
    "Uranus": 19.18797948*const["AU"], 
    "Neptune": 30.06952752*const["AU"],

    #Asteroids
    "Ceres": 2.77*const["AU"], 
    "Vesta": 2.36*const["AU"], 
    "Pallas": 2.77*const["AU"],
    "Hygiea": 3.14*const["AU"], 
    "Interamnia": 3.06*const["AU"], 
    "Eunomia": 2.64*const["AU"],
    "Juno": 2.67*const["AU"], 
    "Davida": 3.17*const["AU"], 
    "Europa": 3.10*const["AU"],
    "Psyche": 2.92*const["AU"],
} 


# Define functions

# Universal Gravitational Force
def F_gravity(ri=1, rj=1, mi=1, mj=1):
    k = const["G"] * mi * mj
    r = rj - ri
    rlength = np.sqrt(np.sum(r*r))
    return k*(1/rlength**3) * r

# Net gravity
def a_net(n, ri=1, rj=1, mi=1, mj=1):
    F = []
    a_net = []
    if n > 2:
        for i in range(n-1):
            #print(f"i={i}")
            Fnet = 0
            for j in range(n):
                if np.array_equal(ri[i], rj[j]):
                    #print(f"j={j}, ri={ri[i]}, rj={rj[j]} Fnet={Fnet}")
                    continue
                else:
                    k = const["G"] * mi[i] * mj[j]
                    r = rj[j] - ri[i]
                    rlength = np.sqrt(np.sum(r*r))
                    Fnet += k*(1/rlength**3) * r
                    #print(f"j={j}, ri={ri[i]}, rj={rj[j]} Fnet={Fnet}")
            F.append(Fnet)
            a_net.append(Fnet/mi[i])
    else:
        k = const["G"] * mi * mj
        r = rj - ri
        rlength = np.sqrt(np.sum(r*r))
        a_net = (k*(1/rlength**3)*r)/mi
    return a_net
 
# Velocity verlet integrator
def integrate(n=1, ri=1, rj=1, mi=1, mj=1, vi=1, dt=1, t=0, t_max=1):
    '''
    This function computes r(t) and v(t) given initial conditions r(0), v(0)
    '''
    if n <= 2:
        r = np.zeros(3)
        v = np.zeros(3)
        a = np.zeros(3)
    else:
        r = np.zeros((n-1, 3))  # (x, y, z)
        v = np.zeros((n-1, 3))  # (vx, vy, vz)
        a = np.zeros((n-1, 3)) # (ax, ay, az)
    
    #initial conditions
    
    ri_0 = ri
    vi_0 = vi
    rj_0 = rj
    r[:] = ri_0
    v[:] = vi_0

    r_values = [ri_0]
    
    v_values = [vi_0]    
    t_values = [t]
    
    N_steps = int(t_max / dt)
    a = a_net(n=n, ri=ri_0, rj=rj_0, mi=mi, mj=mj)
    a_vals = [a]

    for i in range(N_steps):
        t += dt
    
        # velocity Verlet
        if n <= 2:
            v_half = v + 0.5*dt*a    # half-step velocity
            r_new = r + dt*v_half
            a_new = a_net(n=n, ri=r_new, rj=rj_0, mi=mi, mj=mj)
            v_new = v_half + 0.5*dt*a_new
        else:
            v_half = np.zeros((n-1, 3))
            ri_new = np.zeros((n-1, 3))
            vi_new = np.zeros((n-1, 3))
            rj_new = np.zeros((n, 3))
            for i in range(len(a)):
                v_half[i] = v[i] + 0.5*dt*a[i]    # half-step velocity
                ri_new[i] = r[i] + dt*v_half[i]
                rj_new[i+1] = ri_new[i]    
            
            a_new = a_net(n=n, ri=ri_new, rj=rj_new, mi=mi, mj=mj)
            for i in range(len(a_new)):
                vi_new[i] = v_half[i] + 0.5*dt*a_new[i]
    
        r_values.append(ri_new)
        v_values.append(vi_new)
        t_values.append(t)
        a_vals.append(a_new)

        r[:] = ri_new
        v[:] = vi_new
        a = a_new    # important: use force/acceleration for next step!

        # turn lists of results into arrays for easier processing
    r_values = np.array(r_values)
    v_values = np.array(v_values)
    t_values = np.array(t_values)
        
    
    return r_values, v_values, t_values

    
    
# Earth-Sun Trajectory for 1 year
n = 2
dt = 1e-5
ri = np.array([r_j0["Earth"], 0,0])
rj = np.array([r_j0["Sun"], 0, 0])
mi = mass["Earth"]
mj = mass["Sun"]
vi = np.array([0, 6.179*const["AU"], 0])
# 3-body initial conditions
#n = 3
vi_V = np.array([0, (6.179/0.62)*r_j0["Venus"],0])
vi2 = np.array([vi, vi_V])
ri_Sun = np.array([r_j0["Sun"],0,0])
ri_E = np.array([r_j0["Earth"],0,0])
ri_V = np.array([r_j0["Venus"],0,0]) 
ri2 = np.array([ri_E, ri_V])
mi2 = np.array([mass["Earth"], mass["Venus"]])
rj3 = np.array([ri_Sun, ri_E, ri_V])
mj3 = np.array([mass["Sun"], mass["Earth"], mass["Venus"]])
#4-body
vi_Me = np.array([0, (6.179/0.24)*r_j0["Mercury"],0])
ri_Me = np.array([r_j0["Mercury"],0,0])
ri3 = np.array([ri_E, ri_V, ri_Me])
rj4 = np.array([ri_Sun, ri_E, ri_V, ri_Me])
mi3 = np.array([mass["Earth"], mass["Venus"], mass["Mercury"]])
mj4 = np.array([mass["Sun"], mass["Earth"], mass["Venus"], mass["Mercury"]])
vi3 = vi2 = np.array([vi, vi_V, vi_Me])

'''
t = 0
t_max = 1

#Optimization:
#Will need a function to sort the calculated data(r) for each object
#Should take a list of ri and convert each ri to a 3d array
#Should make a mi, rj, and mj based off of ri and the sun
r, v, t = integrate(n=4, ri=ri3, rj=rj4, mi=mi3, mj=mj4, vi=vi3, dt=dt, t_max=t_max)
r_Earth = r[:,0,:]
r_Venus = r[:,1,:]
r_Merc = r[:,2,:]
x = r_Earth[:, 0]/const["AU"]
y = r_Earth[:, 1]/const["AU"]
z = r_Earth[:, 2]/const["AU"]
x_V = r_Venus[:, 0]/const["AU"]
y_V = r_Venus[:, 1]/const["AU"]
z_V = r_Venus[:, 2]/const["AU"]
x_Me = r_Merc[:, 0]/const["AU"]
y_Me = r_Merc[:, 1]/const["AU"]
z_Me = r_Merc[:, 2]/const["AU"]

#t_days = t*convert["years to days"]

#Plot of Earth-Sun Trajectory
plt.plot(x, y, x_V, y_V, x_Me, y_Me)
plt.scatter(x[5], y[5], s=100, color='green', marker='o') # Adds earth on the curve
plt.scatter(x_V[5], y_V[5], s=100, color='orange', marker='o') # Adds venus on the curve
plt.scatter(x_Me[5], y_Me[5], s=100, color='brown', marker='o') # Adds mercury on the curve
plt.scatter(0, 0, s=100, color='yellow', marker='o') # Adds Sun on the curve
plt.title(f"Solar System Simulation for {t_max} years")
plt.xlabel("X (AU)")
plt.ylabel("Y (AU)")
#plt.savefig("Earth-Sun System(1 year)")
plt.show()
plt.close()
'''
'''
#3d plot
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax = plt.axes(projection='3d')
ax.plot3D(x, y, z, x_V, y_V, z_V, x_Me, y_Me, z_Me)
ax.set_xlabel('X (AU)')
ax.set_ylabel('Y (AU)')
ax.set_zlabel('Z (AU)')
plt.show()
'''
# Plot of radial distance between Earth and Sun
#distance_AU = np.sqrt(np.sum(r**2, axis=1))/const["AU"]

#plt.plot(t, distance_AU)
#plt.show()
#'''
'''
#Animation of plot
# Create a figure and axes
fig, ax = plt.subplots()
x_data, y_data = t, distance_AU
line, = ax.plot([], [], lw=2)

# Initialization function
def init():
    ax.set_xlim(min(x_data), x_data[2000]) 
    ax.set_ylim(min(y_data), 1.01)
    return line,

# Animation function
def animate(i, ax=plt.gca()):
    line.set_data(x_data[:i], y_data[:i])
    print(i)

# Create the animation
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=2000, interval=25)

#'''