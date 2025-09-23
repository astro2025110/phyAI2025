# Import all needed packages here
import numpy as np
#%matplotlib inline
#%matplotlib widget
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

#import data
from data import *
#import functions
import functions
from functions import Positions
import code_chris as chris
#run code

planets = ["Earth", "Venus", "Mercury", "Mars", "Jupiter", "Saturn", "Neptune"]
n = len(planets) + 1
ri = np.array([Earth.position, Venus.position, Mercury.position, Mars.position, Jupiter.position, Saturn.position, Neptune.position])

rj = np.array([Sun.position, Earth.position, Venus.position, Mercury.position, Mars.position, Jupiter.position, Saturn.position, Neptune.position])

mi = np.array([Earth.mass, Venus.mass, Mercury.mass, Mars.mass, Jupiter.mass, Saturn.mass, Neptune.mass])

mj = np.array([Sun.mass, Earth.mass, Venus.mass, Mercury.mass, Mars.mass, Jupiter.mass, Saturn.mass, Neptune.mass])

vi = np.array([Earth.velocity, Venus.velocity, Mercury.velocity, Mars.velocity, Jupiter.velocity, Saturn.velocity, Neptune.velocity])


#CM = functions.pos_CM(rj=rj, mj=mj)
#print(CM)
#a_g = functions.a_cg(Earth.position, rj, mj)
#print(Earth.position)
#print(a_g)
t_max = 20
r, v, t = functions.integrate(n=n,ri=ri,vi=vi,mi=mi,rj=rj,mj=mj,dt=1e-1,t=0,t_max=t_max)

for i in range(n-1):
        planets[i] = Positions(planets[i], r[:,i,:])
        
#2D Plot
plt.figure(figsize=(8, 6))
for i in range(n-1):
    plt.plot(planets[i].x, planets[i].y, label=planets[i].name)
    plt.scatter(planets[i].x[15], planets[i].y[15], s=100, marker='o')
plt.scatter(0, 0, s=100, color='yellow', marker='o', label='Sun') # Adds Sun on the curve
plt.title(f"Solar System Simulation for {t_max} years")
plt.xlabel("X (AU)")
plt.ylabel("Y (AU)")
plt.xlim(-33, 33)
plt.ylim(-33, 33)
plt.legend()
plt.savefig(f"{t_max} year Solar System Orbits.png")
plt.show()
plt.close()

#functions.plot_3D(n, planets)
#plt.xlim(-33, 33)
#plt.ylim(-33, 33)
#plt.zlim(-33,33)