import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm

run_num = 111
grid_size = 129
colnames_df =  ["time"] + [f"col{col_num:04d}" for col_num in range(1, grid_size+2)]

data_folder = f"/media/ssudhakar/DATA_10TB/data_solver/case1_{run_num}/"
imgs_folder = f"/media/ssudhakar/DATA_10TB/data_solver/plots_case1/"

df_output_sol_u = pd.read_csv(f"{data_folder}output_sol_u_{run_num}.csv", header=None, names=colnames_df)
df_output_sol_dens = pd.read_csv(f"{data_folder}output_sol_dens_{run_num}.csv", header=None, names=colnames_df)
df_output_sol_temp = pd.read_csv(f"{data_folder}output_sol_temp_{run_num}.csv", header=None, names=colnames_df)
df_output_sol_pres = pd.read_csv(f"{data_folder}output_sol_pres_{run_num}.csv", header=None, names=colnames_df)

u_vel = df_output_sol_u.iloc[:, 1:(grid_size+1)].values
dens = df_output_sol_dens.iloc[:, 1:(grid_size+1)].values
temp = df_output_sol_temp.iloc[:, 1:(grid_size+1)].values
pres = df_output_sol_pres.iloc[:, 1:(grid_size+1)].values

# Assuming x, u_vel, dens, temp, pres, grid_size are defined
x = np.linspace(0, 1, grid_size, endpoint=False)  # Adjusted to use grid_size

# Create figure and subplot axes (4 rows, 1 column)
fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(15, 20), sharex=True)

# Data and labels for each subplot
data_vars = [u_vel, dens, temp, pres]
titles = ['u velocity (time evolution)', 'Density (time evolution)', 'Temperature (time evolution)', 'Pressure (time evolution)']
y_labels = ['u velocity', 'Density', 'Temperature', 'Pressure']

# Plot each variable in its respective subplot
for ax, data, title, ylabel in zip(axes, data_vars, titles, y_labels):
    for i in range(0, 200, 10):
        ax.plot(x, data[i, :], alpha=0.7, label=f"{i * 5e-9:.2e}")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True)
    ax.legend()

# Set common x-label for the bottom subplot
axes[-1].set_xlabel('Domain')

# Adjust layout to prevent overlap
plt.tight_layout()
plt.savefig(f"{imgs_folder}test_img.png", dpi=600)








