import numpy as np
import matplotlib.pyplot as plt

# Parameters
R = 0.1  # Define the radius R
sigma = R / 3  # Standard deviation for Gaussian distribution
num_points = 1000  # Number of points to generate

# Generate x, y coordinates from normal distribution
x = np.random.normal(0, sigma, num_points)
y = np.random.normal(0, sigma, num_points)

# Filter out points that are outside the radius R
mask = np.sqrt(x**2 + y**2) <= R
x_valid = x[mask]
y_valid = y[mask]

# Plotting the scatter plot
plt.figure(figsize=(6, 6))
plt.scatter(x_valid, y_valid, s=5, color='blue')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Scatter Plot of Gaussian Beam Initial Positions')
plt.gca().set_aspect('equal', adjustable='box')
plt.xlim(-R, R)
plt.ylim(-R, R)
plt.grid(True)
plt.savefig("part5gasuusian.png")
plt.close()
def simel(x0,y0):
  import numpy as np
  c=299792458
  e=-1.602176634e-19
  m_e= 9.1093837015e-31 # electron mass, kg
  B0 = 0.0035     # Magnetic field strength (T)
  R = 0.1         # Radius of solenoid (m)
  L = 0.75        # Length of solenoid (m)
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
  Rreq=np.sqrt(        (r_values[:,1])**2+(r_values[:,0])**2 )
  Zreq=r_values[:,2]
  return Zreq,Rreq

plt.figure(figsize=(8,10))
for i in range(1,21):
    ychose=y_valid[i]
    xchose=x_valid[i]
    alfa,beta=simel(xchose,ychose)
    plt.plot(alfa,beta)
plt.xlabel("z-coordinate initialization")
plt.ylabel("r(z)")
plt.savefig("rvsz_part5.png")
