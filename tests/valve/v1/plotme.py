import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv('./work/plant_res.csv')

# Create the figure and axes
fig, axs = plt.subplots(nrows=2,ncols=1)
ax = axs[0]
ax.plot(df["time"], df["valve[1].inlet.m_flow"], label='Linear', color='blue', linestyle='-')
ax.plot(df["time"], df["valve[2].inlet.m_flow"], label='Equi Percent', color='red', linestyle='-')
ax.plot(df["time"], df["valve[3].inlet.m_flow"], label='Fast Acting', color='green', linestyle='-')
ax.set_xlabel('time')
ax.set_ylabel('Mass Flow')
ax.set_title('Three types of valves')
ax.legend()
ax.grid(True)

ax = axs[1]
ax.plot(df["time"], df["valve[1].po"], color='blue', linestyle='-')
ax.plot(df["time"], df["valve[2].po"], color='red', linestyle='-')
ax.plot(df["time"], df["valve[3].po"], color='green', linestyle='-')
ax.set_xlabel('time')
ax.set_ylabel('Valve Opening (%)')
ax.legend()
ax.grid(True) # Add a grid for better readability
plt.show()
