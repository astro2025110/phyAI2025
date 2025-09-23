import matplotlib.pyplot as plt
import numpy as np

# Sample Data
x = np.linspace(0, 10, 100)
y_values = [
    np.sin(x),
    np.cos(x),
    np.sin(x) * np.cos(x),
    np.sin(x)**2
]
labels = ["sin(x)", "cos(x)", "sin(x)cos(x)", "sin(x)^2"]

# Plotting in a loop
plt.figure(figsize=(8, 6))
for i in range(len(y_values)):
    plt.plot(x, y_values[i], label=labels[i])

# Customizations
plt.xlabel("x")
plt.ylabel("y")
plt.title("Multiple Curves")
plt.legend()
plt.grid(True)
plt.show()