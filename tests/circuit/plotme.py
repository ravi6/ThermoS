import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv('./work/plant_res.csv')

# Create the figure and axes
fig, axs = plt.subplots(nrows=5,ncols=1)
ax = axs[0]
ax.plot(df["time"], df["tank.level"], color='blue', linestyle='-')
ax.set_xlabel('time')
ax.set_ylabel('% full')
ax.set_title('Closed Loop Test')
ax.legend()
ax.grid(True)

ax = axs[1]
ax.plot(df["time"], df["tank.Tf"], color='red', linestyle='-')
ax.set_xlabel('time')
ax.set_ylabel('Tank Temp (C)')
ax.legend()
ax.grid(True) 


ax = axs[2]
ax.plot(df["time"], df["delp.pump"]*1e-5, label = 'pump', color='blue', linestyle='-')
ax.plot(df["time"], df["delp.valve"]*1e-5, label = 'valve', color='green', linestyle='-')
ax.plot(df["time"], df["delp.heater"]*1e-5, label = 'heater', color='red', linestyle='-')
ax.set_xlabel('time')
ax.set_ylabel('Press. Drop (kPa)')
ax.legend()
ax.grid(True) 


ax = axs[3]
ax.plot(df["time"], df["pump.inlet.m_flow"], label = 'pump', color='red', linestyle='-')
ax.set_xlabel('time')
ax.set_ylabel('inflow (kg/s)')
ax.legend()
ax.grid(True) 


ax = axs[4]
ax.plot(df["time"], df["pump.Ws"]/1000, label = '', color='blue', linestyle='-')
ax.set_xlabel('time')
ax.set_ylabel('Power (kW)')
ax.legend()
ax.grid(True) 

plt.show()
