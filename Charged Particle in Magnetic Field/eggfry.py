
import numpy as np
import matplotlib.pyplot as plt

#global constants

c=299792458
e=-1.602176634e-19
m_e= 9.1093837015e-31 # electron mass, kg



# Local variables Initial Conditions
#Solenoid Parameters
B0 = 0.0035     # Magnetic field strength (T)
R = 0.1         # Radius of solenoid (m)
L = 0.75        # Length of solenoid (m)

#Spatial Initialization, Velocity of electron

x0 = (3/4) * R  # Initial x position (m)
y0 = 0          # Initial y position (m)
z0 = -0.3       # Initial z position (m)
vx0 = 0         # Initial velocity in x (m/s)
vy0 = 0         # Initial velocity in y (m/s)
vz0 = 0.06*c       # Initial velocity in z (m/s)


def B_func(position,B0=0.0035,L=0.75,R=0.1):
  x,y,z=position[0],position[1],position[2]
  Bz = (B0/2) * ((z / np.sqrt(z**2 + R**2)) - ((z-L) / np.sqrt((z-L)**2 + R**2)))

  denom = np.sqrt(x**2 + y**2)
  if denom == 0:
    Bx=0
    By=0
  else:
    Bx = -((B0*R**2)/4) * (x**2) * (((z**2 + R**2)**(-3/2)) - (((z-L)**2 + R**2)**(-3/2))) / denom
    By = -((B0*R**2)/4) * (y**2) * (((z**2 + R**2)**(-3/2)) - (((z-L)**2 + R**2)**(-3/2))) / denom
  return np.array([Bx, By, Bz])


def lorentz_force(e, v, B):
  return e * (np.cross(v, B) )

r0 = np.array([x0,y0,z0])
v0 = np.array([vx0,vy0,vz0])
Bo=B_func(r0)


dt = 3*10e-10  #  time step  # taken care based on cyclotron frequecy limit for B0
t_max = 1/1000  #  heads on! dont increase it is time consuming!
N_steps = int(t_max / dt)

r = np.zeros(3)  # (x, y,z)
v = np.zeros(3)  # (vx, vy,vz)
t = 0



# initial conditions
r[:] = r0
v[:] = v0

r_values = [r0]
v_values = [v0]
t_values = [t]

# need to start velocity Verlet with current acceleration
a = lorentz_force(e, v0,Bo )/m_e

for i in range(N_steps):
    t += dt

    # velocity Verlet
    v_half = v + 0.5*dt*a        # half-step velocity
    r_new = r + dt*v_half

    B_new=B_func(r_new)
    a_new = lorentz_force(e, v_half,B_new )/m_e  #egg chicken problem resolved, we take from vhalf here. justification necessary
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

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Extract x, y, z coordinates of trajectory
x_traj = r_values[:, 0]
y_traj = r_values[:, 1]
z_traj = r_values[:, 2]

# Cylinder parameters
R = 0.1   # Radius of the cylinder
L = 0.75   # Length of the cylinder

# Generate cylinder
theta = np.linspace(0, 2 * np.pi, 100)  # Circular base angle
z = np.linspace(0, L, 100)  # Height along the z-axis
Theta, Z = np.meshgrid(theta, z)  # Create a grid

# Convert to Cartesian coordinates
X = R * np.cos(Theta)
Y = R * np.sin(Theta)

# Create figure
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

# Plot cylinder surface
ax.plot_surface(X, Y, Z, color='c', alpha=0.3, edgecolor='none')

# Plot the trajectory
ax.plot(x_traj, y_traj, z_traj, marker='o', linestyle='-', color='r', label='Trajectory', linewidth=1)

# Labels and title
ax.set_xlabel("X-axis")
ax.set_ylabel("Y-axis")
ax.set_zlabel("Z-axis")
ax.set_title("3D Trajectory inside a Cylinder")

# Set limits
ax.set_xlim([-R, R])
ax.set_ylim([-R, R])
ax.set_zlim([0, L])

plt.legend()
plt.savefig("trjesol.png")

plt.figure(figsize=(10,2))
plt.plot(r_values[:,0],t_values,color="grey", label="x-coordinate")
plt.plot(r_values[:,1],t_values,linestyle='--',label="y-coordinate")
plt.plot(   np.sqrt((r_values[:,1])**2+(r_values[:,0])**2)   ,t_values,label=r"$r = {\sqrt{x^2 + y^2}}$")
plt.xlabel("time (s)")
plt.ylabel("coordinate component(m)")
plt.legend()
plt.savefig("xyrvst.png")

plt.figure(figsize=(10,2))
plt.plot(r_values[:,2],t_values,color="grey", label="z-coordinate")

plt.xlabel("time (s)")
plt.ylabel("Z coordinate (m)")
plt.legend()
plt.savefig("Zvst.png")



plt.figure(figsize=(10,2))
plt.title("Plot of r as function of z")
plt.plot(  np.sqrt((r_values[:,1])**2+(r_values[:,0])**2) ,r_values[:,2])
plt.xlabel("z")
plt.ylabel(r"$r = {\sqrt{x^2 + y^2}} $")
plt.legend()
plt.savefig("rvsz.png")

#Kinetic Energy
KE_batta=[]  #empty list
KE_initial = 0.5 * m_e * (vx0**2 + vy0**2 + vz0**2)
for count in range(0,len(v_values[:, 0])):
  v1=v_values[count, 0]
  v2=v_values[count, 1]
  v3=v_values[count, 2]
  KE_batta.append(0.5 * m_e * (v1**2 + v2**2 + v3**2))
plt.figure(figsize=(6,7)) 
KE_error = np.abs(np.array(KE_batta) / KE_initial - 1)
plt.plot(        np.arange(1,len(v_values[:, 0])+1),np.log(KE_error)          )
plt.xlabel("integration step")
plt.ylabel(f"$ ln \Delta K /K_0$")
plt.savefig("KErel.png")

