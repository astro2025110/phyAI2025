#global constants
#global constants
import numpy as np
c=299792458
e=-1.602176634e-19
m_e= 9.1093837015e-31 # electron mass, kg

B0 = 0.0035     # Magnetic field strength (T)
R = 0.1         # Radius of solenoid (m)
L = 0.75        # Length of solenoid (m)

# Initial Conditions
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


dt = 3*10e-10  #  time step
t_max = 1/1000
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
    a_new = lorentz_force(e, v_half,B_new )/m_e  #egg chicken problem
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




'''c=299792458
e=-1.602176634e-19
m_e= 9.1093837015e-31 # electron mass, kg
 

def B_func(position):
  x,y,z=position[0],position[1],position[2]
  Bz = (1 / 2000) * ((z / np.sqrt(z**2 + 0.01)) - ((z-0.2) / np.sqrt((z-0.2)**2 + 0.01)))
    
  denom = np.sqrt(x**2 + y**2)
  if denom == 0:
    Bx=0
    By=0
  else:
    Bx = (-0.01 / 4000) * (x**2) * (((z**2 + 0.01)**(-3/2)) - (((z-0.2)**2 + 0.01)**(-3/2))) / denom
    By = (-0.01 / 4000) * (y**2) * (((z**2 + 0.01)**(-3/2)) - (((z-0.2)**2 + 0.01)**(-3/2))) / denom
  return np.array([Bx, By, Bz])
    

def lorentz_force(e, v, B):
  return e * (np.cross(v, B) )

r0 = np.array([0, 0,-0.3])
v0 = np.array([0,0,c*0.06])
Bo=B_func(r0)


dt = 0.5   #  time step
t_max = 10
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
    a_new = lorentz_force(e, v_new,B_new )/m_e  #egg chicken problem
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
t_values = np.array(t_values)'''