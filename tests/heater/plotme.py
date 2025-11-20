import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv('work/plant_res.csv', sep=",")

# Create the figure and axes
fig, axs = plt.subplots(nrows=3,ncols=1)
fig.suptitle('Heater/Cooler Test')
ax = axs[0]
ax.set_ylabel('Temp. (C)')
ax.plot(df["time"], df["htr.Tf"]-273, label='fluid', color='blue', linestyle='-')
ax.plot(df["time"], df["htr.Tw"]-273, label='wall', color='red', linestyle='-')
ax.plot(df["time"], df["Node.T"]-273, label='node', color='green', linestyle='-')
ax.plot(df["time"], df["Node2.T"]-273, label='node2', color='black', linestyle='-')
ax.legend()
ax.grid(True) # Add a grid for better readability

ax = axs[1]
ax.set_ylabel('Flow (lpm)')
ax.plot(df["time"], df["Node2.port[1].m_flow"] * 1000*60, label='Node2A', color='blue', linestyle='-')
ax.plot(df["time"], df["Node2.port[2].m_flow"] * 1000*60, label='Node2B', color='green', linestyle='-')
ax.plot(df["time"], df["Node2.port[3].m_flow"] * 1000*60, label='Node2C', color='red', linestyle='-')

ax = axs[2]
ax.set_ylabel('Flow (lpm)')
ax.plot(df["time"], df["htr.inlet.m_flow"] * 1000*60, label='inlet', color='blue', linestyle='-')
ax.plot(df["time"], -df["htr.outlet.m_flow"] * 1000*60, label='outlet', color='red', linestyle='-')
ax.set_xlabel('time')
ax.legend()
ax.grid(True) # Add a grid for better readability
plt.show()
