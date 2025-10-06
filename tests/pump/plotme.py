import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv('./work/plant_res.csv')

# Create the figure and axes
fig, axs = plt.subplots(nrows=5,ncols=1)
ax = axs[0]
ax.plot(df["time"], df["tank.pFull"], color='blue', linestyle='-')
ax.set_xlabel('time')
ax.set_ylabel('% full')
ax.set_title('Pump Tests')
ax.legend()
ax.grid(True)

ax = axs[1]
ax.plot(df["time"], df["tank.level"], color='blue', linestyle='-')
ax.set_xlabel('time')
ax.set_ylabel('Level (m)')
ax.legend()
ax.grid(True) 


ax = axs[2]
ax.plot(df["time"], df["tank.Tf"]-273, label = 'tank.Tf', color='blue', linestyle='-')
ax.plot(df["time"], df["pump.T"]-273, label = 'pump.T', color='red', linestyle='-')
ax.set_xlabel('time')
ax.set_ylabel('Temp. (C)')
ax.legend()
ax.grid(True) 


ax = axs[3]
ax.plot(df["time"], df["tank.inlet.m_flow"], label = 'tank', color='blue', linestyle='-')
ax.plot(df["time"], df["pump.inlet.m_flow"], label = 'pump', color='red', linestyle='-')
ax.set_xlabel('time')
ax.set_ylabel('inflow (kg/s)')
ax.legend()
ax.grid(True) 


ax = axs[4]
ax.plot(df["time"], df["pump.head"], label = 'pump head', color='blue', linestyle='-')
ax.plot(df["time"], df["dpvh"], label = 'valve pdrop', color='red', linestyle='-')
ax.set_xlabel('time')
ax.set_ylabel('head (m)')
ax.legend()
ax.grid(True) 

plt.show()
