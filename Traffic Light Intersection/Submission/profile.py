import numpy as np
import matplotlib.pyplot as plt
#This module computes v(t) and x(t) for a car depending on the
# decision to brake or to continue driving given an initial
# velocity V0 and position x0 for any t
#functions for v(t) and x(t)
def brake_v(v0, t, a=-3, delta=0.8):
#                           0   if t-delta  < 0
#heaviside(t-delta, t) =    t   if t-delta == 0
#                           1   if t-delta  > 0
    t1 = np.heaviside(delta-t, 1)
    v1 = v0*t1
    #Add deceleration component here
    delta2 = delta*np.heaviside(t-delta, 0)
    t2 = t*np.heaviside(t-delta, 0) - delta2
    v2 = v0*np.heaviside(t-delta, 0) + a*t2
    #put it all together
    v = v1+v2
    v = v*np.heaviside(v, v)
    return v

def brake_x(v0, t, x0, a=-3, delta=0.8):
    #Begin with t<delta. x should increase by the same increment here
    t1 = np.heaviside(delta-t, 1)
    x1 = x0*t1 + (v0*t1)*t
    v = brake_v(v0, t)
    del_v = v - v0*np.heaviside(t, 1)
    #Add deceleration component here. Start with time
    t2 = t*np.heaviside(t-delta, 0) - delta*np.heaviside(t-delta, 0)
    brake = np.where(t == delta)[0]
    x2 = x1[brake]*np.heaviside(t-delta, 0) + v0*t2 + 0.5*a*(t2**2)
    x3 = x1+x2
    #Insert stopping point. This is when V=0 or the car comes to a stop.
    #Replace the displacements beyond the stopping point with xstop
    stop = np.where(v == 0)[0]
    stop = stop[0]
    rest = x3[stop]
    xf1 =  x3*np.heaviside(t[stop]-t, 1)
    xf2 = rest*np.heaviside(t-t[stop], 0)
    xf = xf1+xf2
    return xf

def drive_x(v0, t, x0, a=-3, delta=0.8):
    t1 = np.heaviside(delta-t, 1)
    x1 = (x0 + v0*t)*t1
    t2 = t*np.heaviside(t-delta, 0) - delta*np.heaviside(t-delta, 0)
    i = np.where(t == delta)[0]
    x2 = x1[i]*np.heaviside(t-delta, 0) + v0*t2
    return x1+x2
#Plot v(t) and x(t). Should be 2 subplots for each x0
t = np.arange(0, 7, 0.1)

x0 = -30
v0 = 15.0
tau = 3
W = 30

v_brake = brake_v(v0, t)
v_drive = v0*np.heaviside(t, 1)
x_brake = brake_x(v0, t, x0)
x_drive = drive_x(v0, t, x0)


fig, axes = plt.subplots(1, 2, figsize=(10, 6))
fig.suptitle(f"Motion of car for x0 = {x0} m")

axes[0].plot(t, v_brake, 'b', label='brake') 
axes[1].plot(t, x_brake, 'b', label='brake')
axes[0].plot(t, v_drive, 'r', label='drive')
axes[1].plot(t, x_drive, 'r', label='drive')
axes[0].plot([tau, tau], [x0, 1.05*W], "--", color="red", lw=1)
axes[1].fill_between([0, 7], [W, W], color="black", alpha=0.3)
axes[1].plot([tau, tau], [x0, 1.05*W], "--", color="red", lw=1)

axes[0].set_xlim(0, 6.9)
axes[0].set_ylim(0, v0+10)
axes[0].set_xlabel("time t (s)")
axes[0].set_ylabel("speed v(t) (m/s)")
axes[0].legend(loc="best")   
axes[1].set_xlabel("time t (s)")
axes[1].set_ylabel("position x(t) (m)")
axes[1].set_ylim(x0, 80)
axes[1].set_xlim(0, 6.9)
axes[1].legend(loc="best")
plt.tight_layout()
plt.savefig("time_series1.png")
plt.close()

