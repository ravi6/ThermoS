import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv('work/plant_res.csv', sep=",")

# Create the figure and axes
fig, axs = plt.subplots(nrows=3,ncols=1)
fig.suptitle('Heater + Tank test')
ax = axs[0]
ax.set_ylabel('Flow (kg/s)')
ax.plot(df["time"], df["tank.inlet.m_flow"], label='tankIn', color='blue', linestyle='-')
ax.plot(df["time"], -df["tank.outlet.m_flow"], label='tankOut', color='red', linestyle='-')
ax.plot(df["time"], df["htr.inlet.m_flow"], label='htrIn', color='black', linestyle='-')
ax.plot(df["time"], -df["htr.outlet.m_flow"], label='htrOut', color='green', linestyle='-')
ax.legend()
ax.grid(True) # Add a grid for better readability
ax.set_xlabel('time')
ax.legend()
ax.grid(True) # Add a grid for better readability



ax = axs[1]
ax.set_ylabel('Temp (C)')
ax.plot(df["time"], df["tank.T"]-273, label='tank', color='blue', linestyle='-')
ax.plot(df["time"], df["htr.Tf"]-273,  label='htr', color='red', linestyle='-')
ax.legend()
ax.grid(True) # Add a grid for better readability
ax.set_xlabel('time')
ax.legend()
ax.grid(True) # Add a grid for better readability


ax = axs[2]
ax.set_ylabel('Pressure (bar)')
ax.plot(df["time"], df["tank.p"]*1e-5, label='tank', color='blue', linestyle='-')
ax.plot(df["time"], df["htr.p"]*1e-5, label='htr', color='red', linestyle='-')
ax.legend()
ax.grid(True) # Add a grid for better readability
ax.set_xlabel('time')
ax.legend()
ax.grid(True) # Add a grid for better readability


plt.show()
