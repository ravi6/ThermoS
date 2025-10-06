import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv('./work/plant_res.csv')

# Create the figure and axes
fig, axs = plt.subplots(nrows=3,ncols=1)
ax = axs[0]
ax.plot(df["time"], df["hx.seg_c.Tf"], label = "cold", color='blue', linestyle='-')
ax.plot(df["time"], df["hx.seg_h.Tf"], label = "hot", color='red', linestyle='-')
ax.plot(df["time"], df["hx.Tw"], label = "Wall", color='green', linestyle='-')
ax.set_xlabel('time')
ax.set_ylabel('Temp. (C)')
ax.set_title('Hx_Seg Test')
ax.legend()
ax.grid(True)

ax = axs[1]
ax.plot(df["time"], df["hx.Qwf_c"], label="cold", color='blue', linestyle='-')
ax.plot(df["time"], df["hx.Qwf_h"], label="hot", color='red', linestyle='-')
ax.set_xlabel('time')
ax.set_ylabel('Heat Loss/Gain kW')
ax.legend()
ax.grid(True) 

ax = axs[2]
ax.plot(df["time"], df["hx.portA_c.m_flow"], label = 'cold', color='blue', linestyle='-')
ax.plot(df["time"], df["hx.portA_h.m_flow"], label = 'hot', color='red', linestyle='-')
ax.set_xlabel('time')
ax.set_ylabel('Mass Flow (kg/s)')
ax.legend()
ax.grid(True) 

plt.show()
