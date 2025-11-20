import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

df = pd.read_csv('./work/plant_res.csv')

# Create the figure and axes
fig, axs = plt.subplots(nrows=2,ncols=1)
ax = axs[0]
ax.plot(df["time"], df["tank.T"]-273, label = "Tank T", color='red', linestyle='-')
ax.set_xlabel('time')
ax.set_ylabel('Temp. (C)')
ax.set_title('OB tank tests')
#ax.legend()
ax.grid(True)

ax = axs[1]
ax.plot(df["time"], df["tank.p"]/1e5, label="", color='blue', linestyle='-')
ax.set_xlabel('time')
ax.set_ylabel('Pressure kPa')
#ax.legend()
ax.grid(True) 
plt.show()
