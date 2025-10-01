import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv('./work/plant_res.csv', sep=",")

# Create the figure and axes
fig, axs = plt.subplots(nrows=3,ncols=1)
ax = axs[0]
ax.plot(df["time"], df["v1.inlet.m_flow"], label='V1', color='blue', linestyle='-')
ax.plot(df["time"], df["v2.inlet.m_flow"], label='V2', color='red', linestyle='-')
ax.set_xlabel('time')
ax.set_ylabel('Mass Flow')
ax.set_title('Parallel valves from Port Mixer ')
ax.legend()
ax.grid(True) # Add a grid for better readability

ax = axs[1]
ax.plot(df["time"], df["Node.m"], color='blue', linestyle='-')
ax.set_xlabel('time')
ax.set_ylabel('Node Flow Rate (kg/s)')
ax.legend()
ax.grid(True) # Add a grid for better readability

ax = axs[2]
ax.plot(df["time"], df["Node.Q_in"], color='blue', linestyle='-')
ax.set_xlabel('time')
ax.set_ylabel('Node Heat Flow (J/s)')
ax.legend()
ax.grid(True) # Add a grid for better readability
plt.show()
