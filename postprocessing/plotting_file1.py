import numpy as np
# import matplotlib
# matplotlib.use('TkAgg')  # Use an interactive backend
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pandas as pd
import json

with open("../code/configuration.json", "r") as file:
    config1 = json.load(file, )

for key, value in zip(config1.keys(), config1.values()):
    print(f"{key:<30s}", f"{value:>}")

colnames = ["time"] + [f"col{col_num:04d}" for col_num in range(1, config1["Nx"]+2)]
df_output_sol_u = pd.read_csv("../postprocessing/output_sol_u_002.csv", header=None, names=colnames)
df_output_sol_dens = pd.read_csv("../postprocessing/output_sol_dens_002.csv", header=None, names=colnames)
df_output_sol_temp = pd.read_csv("../postprocessing/output_sol_temp_002.csv", header=None, names=colnames)
df_output_sol_pres = pd.read_csv("../postprocessing/output_sol_pres_002.csv", header=None, names=colnames)

u_vel = df_output_sol_u.iloc[:, 1:(config1["Nx"]+1)].values
dens = df_output_sol_dens.iloc[:, 1:(config1["Nx"]+1)].values
temp = df_output_sol_temp.iloc[:, 1:(config1["Nx"]+1)].values
pres = df_output_sol_pres.iloc[:, 1:(config1["Nx"]+1)].values


# Assuming x, u_vel, dens, temp, pres, and config1["Nx"] are defined
x = np.linspace(0, 1, config1["Nx"], endpoint=False)  # Adjusted to use config1["Nx"] for generality

# Create figure and subplot axes (4 rows, 1 column)
fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(15, 20), sharex=True)

# Data and labels for each subplot
data_vars = [u_vel, dens, temp, pres]
titles = ['u velocity (time evolution)', 'Density (time evolution)', 'Temperature (time evolution)', 'Pressure (time evolution)']
y_labels = ['u velocity', 'Density', 'Temperature', 'Pressure']

# Plot initialization function (called at the start of the animation)
def init():
    for ax in axes:
        ax.clear()  # Clear the axes for the initial setup
    return []

# Update function to animate each time step
def update(frame):
    # Clear previous data
    for ax, data, title, ylabel in zip(axes, data_vars, titles, y_labels):
        ax.clear()
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(True)
        
        # Plot new data at the current time step
        ax.plot(x, data[frame, :], alpha=0.7, label=f"t = {frame * 5e-9:.2e}")
        ax.legend(loc="best")
        
    # Set the common x-label only once for the last subplot
    axes[-1].set_xlabel('Domain')
    
    return [ax for ax in axes]

# Create the animation
ani = FuncAnimation(fig, update, frames=range(0, 10000, 1), init_func=init, blit=True, interval=200)

# Show the animation
plt.tight_layout()
plt.show()