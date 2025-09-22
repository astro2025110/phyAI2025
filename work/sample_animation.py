import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

# Create a figure and axes
fig, ax = plt.subplots()
line, = ax.plot([], [], lw=2)
ax.set_xlim(0, 10)
ax.set_ylim(-1, 1)

# Initialization function
def init():
    line.set_data([], [])
    return line,

# Animation function
def animate(i):
    x = np.linspace(0, 10, 1000)
    y = np.sin(x + i/10.0)
    line.set_data(x, y)
    return line,

# Create the animation
ani = animation.FuncAnimation(fig, animate, init_func=init, frames=200, interval=20, blit=True)

plt.show()
plt.close()